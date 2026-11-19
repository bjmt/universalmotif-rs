//! Discovery of gapped motifs with [`gapped_motif_finder`]: two parts with a
//! spacer of varying length between them.
//!
//! Pairs of k-mers at every spacing searched are counted with
//! [`spaced_pairs`](crate::finder::spaced_pairs), and the most enriched become seed motifs, which are then
//! refined by rescanning, trimmed, given the spacer sizes their hits support
//! and grown into neighbouring informative columns. A motif is kept only if
//! its hits favour some spacings more than independent placement of its two
//! parts would, which is tested in sequences set aside from the search, so
//! that two motifs which merely occur in the same sequences are not reported
//! as one.
//!
//! The finder is experimental: its options, its results and the motifs it
//! finds may change in later releases.

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::create::{CreateOptions, from_strings, from_typed_matrix};
use crate::error::{Result, arg_err};
use crate::finder::{NegativePairs, spaced_pairs_with};
use crate::mask::{DustOptions, mask_complexity};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};
use crate::ops::{add_gap, normalize};
use crate::packed::EncodedSequences;
use crate::progress;
use crate::pvalue::with_threads;
use crate::scan::{HitTable, NoOverlaps, ScanOptions, scan_encoded};
use crate::shuffle::{ShuffleMethod, ShuffleOptions, shuffled_backgrounds};
use crate::stats::ln_fisher_greater;

/// Options for [`gapped_motif_finder`].
#[derive(Debug, Clone, PartialEq)]
pub struct GappedFinderOptions {
    /// The k-mer sizes to seed motifs from, each from 2 to 6; 4, 5 and 6 by
    /// default.
    pub k: Vec<usize>,
    /// The shortest spacer searched; 0 by default.
    pub min_gap: usize,
    /// The longest spacer searched, at most 31, so that the motifs found
    /// can be compared with their gaps; 31 by default.
    pub max_gap: usize,
    /// The most motifs to find; 5 by default.
    pub nmotifs: usize,
    /// P-value for a site to count as a hit; 1e-4 by default.
    pub hit_pvalue: f64,
    /// A motif is kept if the P-value of its spacing, corrected for the
    /// pairs and seeds tried, is at most this, and its hits are enriched
    /// with at most this P-value; 1e-3 by default.
    pub stop_pvalue: f64,
    /// Letters added on either side of each part of a seed; 3 by default.
    pub flank: usize,
    /// Columns with fewer bits than this are trimmed from the ends of each
    /// part, and parts grow into columns with at least this many; 0.5 by
    /// default.
    pub min_ic: f64,
    /// The most enriched pairs of k-mers built into motifs for each k-mer
    /// size and motif; 6 by default.
    pub n_seed_pairs: usize,
    /// Pairs are tested only if found in at least this fraction of the
    /// sequences, and in at least 5; 0.02 by default. Lower values find rarer
    /// motifs at the cost of more tests to correct for.
    pub min_fraction: f64,
    /// The fraction of the sequences, and of the negatives, set aside to test
    /// the motifs found in the rest; 0.3 by default. The P-values then come
    /// from sequences the motifs were not fitted to and are corrected only
    /// for the motifs tested, which keeps them honest however many seeds
    /// are tried. With 0, the P-values come from the sequences the motifs
    /// were found in, corrected for the pairs and seeds tried, and are
    /// somewhat optimistic.
    pub holdout: f64,
    /// Whether pairs of 6-mers are counted in each sequence holding them
    /// with at most one mismatch, which finds the spacing of rarer motifs
    /// better but takes several times longer; [`Mismatch::Auto`] by
    /// default. Shorter k-mers always match exactly.
    pub mismatch: Mismatch,
    /// Mask low-complexity regions, such as homopolymers and short tandem
    /// repeats, in the sequences and the negatives before the search, with
    /// DUST at its default settings; on by default. Their repeated units
    /// otherwise look like parts with a favoured spacing.
    pub mask: bool,
    /// Search both strands; on by default.
    pub rc: bool,
    /// Shuffling for the negative sequences when none are given; 2 by
    /// default.
    pub shuffle_k: usize,
    /// How many shuffles of the sequences make up the negatives, when none
    /// are given; 3 by default, which steadies the background so that
    /// motifs in few of the sequences are found whatever the seed.
    pub bkg_shuffles: usize,
    /// Seed for shuffling the sequences into negatives.
    pub seed: u64,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for GappedFinderOptions {
    fn default() -> Self {
        GappedFinderOptions {
            k: vec![4, 5, 6],
            min_gap: 0,
            max_gap: 31,
            nmotifs: 5,
            hit_pvalue: 1e-4,
            stop_pvalue: 1e-3,
            flank: 3,
            min_ic: 0.5,
            n_seed_pairs: 6,
            min_fraction: 0.02,
            holdout: 0.3,
            mismatch: Mismatch::Auto,
            mask: true,
            rc: true,
            shuffle_k: 2,
            bkg_shuffles: 3,
            seed: 1,
            nthreads: 1,
        }
    }
}

/// When [`gapped_motif_finder`] lets 6-mers match with one mismatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Mismatch {
    /// Only exact matches.
    Off,
    /// Always allow one mismatch.
    On,
    /// Search with exact matches, and search again allowing one mismatch,
    /// at half of `stop_pvalue`, only if no motif is found, so that the
    /// slower search is run only when the faster one has nothing to show.
    #[default]
    Auto,
}

/// A motif found by [`gapped_motif_finder`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct GappedMotif {
    /// The motif, a gapped PPM named `gapped_1`, `gapped_2` and so on.
    pub motif: Motif,
    /// The k-mer size of the pair it grew from.
    pub k: usize,
    /// The pair of k-mers it grew from, with the spacings they were counted
    /// over, as `ACGT[0-15]TGCA`.
    pub seed_pair: String,
    /// Positive sequences with at least one hit.
    pub seqs_pos: u64,
    /// Negative sequences with at least one hit.
    pub seqs_neg: u64,
    /// The number of positive sequences.
    pub n_pos: usize,
    /// The number of negative sequences.
    pub n_neg: usize,
    /// Pairs of k-mers tested at the size of the pair it grew from.
    pub candidates: u64,
    /// Pairs built into motifs at that size.
    pub seed_pairs: usize,
    /// The one-sided Fisher's exact test P-value of sequences with hits
    /// being more common among the positives: in those held out,
    /// uncorrected, or else in all of them, corrected as the spacing is.
    pub pvalue: f64,
    /// Base 10 log of `pvalue`, which stays finite when `pvalue` is too small
    /// for a double.
    pub log10_pvalue: f64,
    /// The P-value of the motif's spacing against independent placement of
    /// its parts and against the negatives: in the sequences held out,
    /// corrected for the motifs tested there, or else in all of them,
    /// corrected for the pairs and seeds tried at its size and for the
    /// number of sizes; capped at 1.
    pub spacing_pvalue: f64,
    /// Base 10 log of `spacing_pvalue`.
    pub log10_spacing_pvalue: f64,
}

/// A candidate motif and what it was judged on.
struct Candidate {
    motif: Motif,
    hits: Vec<(usize, usize, usize)>,
    held_out_hits: Vec<(usize, usize, usize)>,
    seqs_pos: u64,
    seqs_neg: u64,
    ln_spacing: f64,
    ln_p: f64,
    k: usize,
    seed: String,
    n_tested: usize,
    tried: usize,
}

/// Sequences as text, for cutting out sites, and encoded, for scanning.
struct SeqSet {
    text: Vec<String>,
    encoded: EncodedSequences,
}

impl SeqSet {
    fn new(text: Vec<String>, nthreads: usize) -> Result<SeqSet> {
        let encoded = EncodedSequences::new(&text, &Alphabet::Dna, nthreads)?;
        Ok(SeqSet { text, encoded })
    }
}

/// The hits of a motif's parts by spacing; see `Search::profile`.
struct Profile {
    hits: HitTable,
    spacers: Vec<usize>,
    observed: Vec<f64>,
    expected: Vec<f64>,
}

/// Sequences split at random into those kept and a `fraction` held out.
fn split(seqs: Vec<String>, fraction: f64, seed: u64) -> (Vec<String>, Vec<String>) {
    let (mut kept, mut held) = (Vec::new(), Vec::new());
    for (i, s) in seqs.into_iter().enumerate() {
        // SplitMix64 of the seed and index, as a number in [0, 1).
        let mut z = seed.wrapping_add((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        if ((z >> 11) as f64 / (1u64 << 53) as f64) < fraction {
            held.push(s);
        } else {
            kept.push(s);
        }
    }
    (kept, held)
}

/// A set with the letters of some hits replaced by N.
fn masked(set: &mut SeqSet, hits: &[(usize, usize, usize)], nthreads: usize) -> Result<SeqSet> {
    let mut text = std::mem::take(&mut set.text);
    for &(i, s, e) in hits {
        text[i].replace_range(s..e, &"N".repeat(e - s));
    }
    SeqSet::new(text, nthreads)
}

/// Everything the search steps share.
struct Search<'a> {
    opts: &'a GappedFinderOptions,
    work: SeqSet,
    neg: &'a SeqSet,
    /// Positives and negatives set aside to test the motifs.
    held_out: Option<(SeqSet, SeqSet)>,
    /// The negatives' pair counts for each k-mer size, kept between searches
    /// as the negatives never change.
    neg_pairs: RefCell<HashMap<usize, NegativePairs>>,
}

/// Find motifs made of two parts with a spacer of varying length.
/// Experimental, as the [module docs](self) say.
///
/// A share of the sequences and of the negatives, `holdout`, is first set
/// aside. In the rest, pairs of k-mers at every spacing from `min_gap` to
/// `max_gap` found in at least `min_fraction` of the sequences are counted in
/// the sequences and in the negatives, by default dinucleotide shuffles of
/// the sequences, and tested for enrichment with Fisher's exact test, for
/// each k-mer size in `k`. For each size, each of the `n_seed_pairs` most enriched
/// distinct pairs, over the window of spacings that gives it, seeds a motif:
/// its sites, with `flank` letters on either side of each part, make a gapped
/// motif, refined by rescanning and trimmed of columns below `min_ic` bits at
/// either end. The spacer sizes kept are those whose hits outnumber both the
/// negatives and what the two parts would give if they were placed
/// independently of each other. Each part then grows into informative
/// columns, a motif that ends up with a single part is tried again as a
/// direct and as an inverted repeat of itself, and motifs left without a gap
/// are dropped.
///
/// Each motif is then tested in the sequences held out: the P-value of its
/// spacing, against independent placement of its parts and against the
/// negatives at each spacing, multiplied by the number of motifs tested, must
/// be at most `stop_pvalue`, and sequences with hits must be enriched over
/// the negatives at 0.05. Of those that pass, the one with the strongest
/// spacing is kept, its hits are masked and the search repeated. As the
/// tests use sequences the motifs were not fitted to, the P-values stay
/// honest however widely the search looks. With `holdout` at 0, both
/// P-values come from all the sequences instead, corrected for the pairs and
/// seeds tried, and are somewhat optimistic.
///
/// With `mismatch` at [`Mismatch::Auto`], the default, a search that finds
/// nothing is run again with 6-mers allowed one mismatch, which can recover
/// the spacing of motifs whose parts vary too much to be counted exactly.
/// The second search keeps motifs at half of `stop_pvalue`, so that the two
/// together report no more false motifs than one search would.
///
/// Low-complexity regions are masked first with `mask`. For genomic
/// sequences, negatives drawn from the genome and matched for composition
/// give far better calibrated P-values than the default shuffles, which lose
/// the loose repeats real genomes are full of; their repeated units then
/// look like parts with a favoured spacing. Spacings are also judged against
/// how often the parts meet at each spacing in the negatives, so that parts
/// which run together there too, as A-rich stretches do, are not taken for
/// a gapped motif.
///
/// Sequences are read as DNA, in either case; letters other than A, C, G and
/// T never match.
///
/// # References
///
/// - Morgulis A, Gertz EM, Schäffer AA, Agarwala R (2006). A fast and symmetric
///   DUST implementation to mask low-complexity DNA sequences. *Journal of
///   Computational Biology* 13(5):1028–1040.
///   <https://doi.org/10.1089/cmb.2006.13.1028>
pub fn gapped_motif_finder<S: AsRef<str> + Sync>(
    sequences: &[S],
    negatives: Option<&[S]>,
    opts: &GappedFinderOptions,
) -> Result<Vec<GappedMotif>> {
    if opts.mismatch != Mismatch::Auto {
        return search_once(sequences, negatives, opts);
    }
    let exact = search_once(
        sequences,
        negatives,
        &GappedFinderOptions {
            mismatch: Mismatch::Off,
            ..opts.clone()
        },
    )?;
    if !exact.is_empty() || !opts.k.contains(&6) {
        return Ok(exact);
    }
    search_once(
        sequences,
        negatives,
        &GappedFinderOptions {
            mismatch: Mismatch::On,
            // Halved, so that two searches of sequences with nothing to
            // find report a motif no more often than one would.
            stop_pvalue: opts.stop_pvalue / 2.0,
            ..opts.clone()
        },
    )
}

fn search_once<S: AsRef<str> + Sync>(
    sequences: &[S],
    negatives: Option<&[S]>,
    opts: &GappedFinderOptions,
) -> Result<Vec<GappedMotif>> {
    let progress = progress::current();
    if sequences.is_empty() {
        return arg_err("no sequences given");
    }
    if opts.k.is_empty() || opts.k.iter().any(|k| !(2..=6).contains(k)) {
        return arg_err("each k must be from 2 to 6");
    }
    if opts.min_gap > opts.max_gap || opts.max_gap > 31 {
        return arg_err("gaps must satisfy min_gap <= max_gap <= 31");
    }
    if !(opts.hit_pvalue > 0.0 && opts.hit_pvalue < 1.0)
        || !(opts.stop_pvalue > 0.0 && opts.stop_pvalue <= 1.0)
    {
        return arg_err("hit_pvalue must be in (0, 1) and stop_pvalue in (0, 1]");
    }
    if !(0.0..=1.0).contains(&opts.min_fraction) {
        return arg_err("min_fraction must be in [0, 1]");
    }
    if !(0.0..0.9).contains(&opts.holdout) {
        return arg_err("holdout must be in [0, 0.9)");
    }
    let pos: Vec<String> = sequences.iter().map(|s| dna_text(s.as_ref())).collect();
    let neg: Vec<String> = match negatives {
        Some(n) => n.iter().map(|s| dna_text(s.as_ref())).collect(),
        None => shuffled_backgrounds(
            &pos,
            &ShuffleOptions {
                k: opts.shuffle_k,
                method: ShuffleMethod::Euler,
                seed: opts.seed,
                window: None,
                nthreads: opts.nthreads,
            },
            opts.bkg_shuffles,
        )?,
    };
    if neg.is_empty() {
        return arg_err("no negative sequences");
    }
    // Masked after shuffling, as the masked letters would otherwise be
    // shuffled into the background.
    let (pos, neg) = if opts.mask {
        let dust = DustOptions {
            letter: 'N',
            nthreads: opts.nthreads,
            ..Default::default()
        };
        (mask_complexity(&pos, &dust)?, mask_complexity(&neg, &dust)?)
    } else {
        (pos, neg)
    };
    let (n_pos, n_neg) = (pos.len(), neg.len());
    let (pos, neg, held_out) = if opts.holdout > 0.0 {
        let (pos, test_pos) = split(pos, opts.holdout, opts.seed);
        let (neg, test_neg) = split(neg, opts.holdout, opts.seed ^ 0x9e37_79b9_7f4a_7c15);
        if test_pos.is_empty() || test_neg.is_empty() || pos.len() < 2 || neg.is_empty() {
            return arg_err("too few sequences to hold some out; set holdout to 0");
        }
        let test = (
            SeqSet::new(test_pos, opts.nthreads)?,
            SeqSet::new(test_neg, opts.nthreads)?,
        );
        (pos, neg, Some(test))
    } else {
        (pos, neg, None)
    };
    let min_count = 5.max((opts.min_fraction * pos.len() as f64) as u32);
    let mut ks = opts.k.clone();
    ks.sort_unstable();
    ks.dedup();
    let neg = SeqSet::new(neg, opts.nthreads)?;
    let mut search = Search {
        opts,
        work: SeqSet::new(pos, opts.nthreads)?,
        neg: &neg,
        held_out,
        neg_pairs: RefCell::new(HashMap::new()),
    };
    let ln_stop = opts.stop_pvalue.ln();
    let mut out = Vec::new();
    progress::begin(&progress, "searching", (opts.nmotifs * ks.len()) as u64);
    let mut searched = 0;
    for found in 0..opts.nmotifs {
        let name = format!("gapped_{}", found + 1);
        let mut candidates = Vec::new();
        for &k in &ks {
            progress::check(&progress)?;
            candidates.extend(search.seed_candidates(k, min_count, &name)?);
            searched += 1;
            progress::add(&progress, 1);
        }
        if candidates.is_empty() {
            break;
        }
        // The spacing is the evidence for a gapped motif, as the hits
        // expected at each spacing include those the negatives have there.
        // It is tested in the sequences held out, corrected for the motifs
        // tested, and the hits there need only be enriched over the
        // negatives at 0.05. Without sequences held out, both are tested in
        // those the motifs were found in, corrected for the pairs and seeds
        // tried at each k-mer size, each size taking an equal share of
        // stop_pvalue.
        let n_candidates = candidates.len() as f64;
        let ln_enriched = if search.held_out.is_some() {
            0.05f64.ln()
        } else {
            ln_stop
        };
        for c in &mut candidates {
            if search.held_out.is_some() {
                search.test(c)?;
                c.ln_spacing += n_candidates.ln();
            } else {
                c.ln_spacing +=
                    (c.n_tested as f64).ln() + (c.tried as f64).ln() + (ks.len() as f64).ln();
            }
        }
        // A single part's spacing comes from pairs within it, which says
        // nothing about a gap, so motifs that kept two parts come first.
        let mut best: Option<Candidate> = None;
        for c in candidates {
            if c.ln_spacing > ln_stop || c.ln_p > ln_enriched {
                continue;
            }
            let better = match &best {
                None => true,
                Some(b) => {
                    (!c.motif.gap.isgapped, c.ln_spacing) < (!b.motif.gap.isgapped, b.ln_spacing)
                }
            };
            if better {
                best = Some(c);
            }
        }
        let Some(best) = best else { break };
        // Mask the hits before the next search.
        search.work = masked(&mut search.work, &best.hits, opts.nthreads)?;
        if let Some((test_pos, _)) = &mut search.held_out {
            *test_pos = masked(test_pos, &best.held_out_hits, opts.nthreads)?;
        }
        out.push(GappedMotif {
            motif: best.motif,
            k: best.k,
            seed_pair: best.seed,
            seqs_pos: best.seqs_pos,
            seqs_neg: best.seqs_neg,
            n_pos,
            n_neg,
            candidates: best.n_tested as u64,
            seed_pairs: best.tried,
            pvalue: best.ln_p.min(0.0).exp(),
            log10_pvalue: best.ln_p.min(0.0) / std::f64::consts::LN_10,
            spacing_pvalue: best.ln_spacing.min(0.0).exp(),
            log10_spacing_pvalue: best.ln_spacing.min(0.0) / std::f64::consts::LN_10,
        });
    }
    // A search that stopped early counts the searches it skipped.
    progress::add(
        &progress,
        (opts.nmotifs * ks.len()).saturating_sub(searched) as u64,
    );
    Ok(out)
}

impl Search<'_> {
    /// Motifs built from the most enriched pairs of k-mers.
    fn seed_candidates(&self, k: usize, min_count: u32, name: &str) -> Result<Vec<Candidate>> {
        let o = self.opts;
        // One mismatch in a pair of shorter k-mers matches too much of the
        // background to help.
        let mismatch = o.mismatch == Mismatch::On && k == 6;
        let mut cache = self.neg_pairs.borrow_mut();
        let pairs = spaced_pairs_with(
            &self.work.text,
            &self.neg.text,
            k,
            o.min_gap,
            o.max_gap,
            o.rc,
            min_count,
            mismatch,
            o.nthreads,
            Some(cache.entry(k).or_default()),
        )?;
        drop(cache);
        let mut candidates: Vec<Candidate> = Vec::new();
        let mut covered: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let mut seen = HashSet::new();
        let mut tried = 0;
        for p in &pairs {
            if tried == o.n_seed_pairs {
                break;
            }
            if !seen.insert((p.left, p.right)) {
                continue;
            }
            let (lo, hi) = (p.gap, p.gap + p.width - 1);
            let (left, right) = (kmer_string(p.left, k), kmer_string(p.right, k));
            let occ = occurrences(
                &self.work.text,
                &left,
                &right,
                lo,
                hi,
                o.rc,
                mismatch,
                o.nthreads,
            )?;
            // Seeds lying mostly within the hits of a motif already built are
            // the same motif again.
            let inside = occ
                .iter()
                .filter(|&&(i, start, _, _)| {
                    covered
                        .get(&i)
                        .is_some_and(|v| v.iter().any(|&(s, e)| s <= start && start < e))
                })
                .count();
            if !occ.is_empty() && inside * 2 > occ.len() {
                continue;
            }
            tried += 1;
            let Some(mut cand) = self.build(&occ, k, name)? else {
                continue;
            };
            cand.seed = format!("{left}[{lo}-{hi}]{right}");
            for &(i, s, e) in &cand.hits {
                covered.entry(i).or_default().push((s, e));
            }
            candidates.push(cand);
        }
        for c in &mut candidates {
            c.n_tested = pairs.len();
            c.tried = tried;
        }
        Ok(candidates)
    }

    /// Scan with a motif, its pseudocount re-applied, keeping the best of
    /// overlapping hits.
    fn scan(&self, motif: &Motif, negatives: bool) -> Result<HitTable> {
        self.scan_in(motif, if negatives { self.neg } else { &self.work })
    }

    fn scan_in(&self, motif: &Motif, seqs: &SeqSet) -> Result<HitTable> {
        let seqs = &seqs.encoded;
        let opts = ScanOptions {
            pvalue: self.opts.hit_pvalue,
            rc: self.opts.rc,
            no_overlaps: Some(NoOverlaps::default()),
            nthreads: self.opts.nthreads,
            ..Default::default()
        };
        // A step of the search, not a scan to report on its own.
        progress::without_progress(|| scan_encoded(&[normalize(motif)?], seqs, &opts))
    }

    /// The letters of each hit, read along the motif, leaving out those with
    /// letters other than A, C, G and T.
    fn hit_sites(&self, hits: &HitTable, keep: impl Fn(usize) -> bool) -> Vec<String> {
        let mut sites = Vec::with_capacity(hits.len());
        for i in 0..hits.len() {
            if !keep(i) {
                continue;
            }
            let seq = &self.work.text[hits.sequence[i] as usize];
            let site = &seq[hits.start[i] as usize..hits.end[i] as usize];
            sites.push(if hits.minus[i] {
                rc_string(site)
            } else {
                site.to_string()
            });
        }
        clean(sites)
    }

    /// A gapped motif grown from the occurrences of a seed pair, with its
    /// hits and enrichment.
    fn build(&self, occ: &[Occurrence], k: usize, name: &str) -> Result<Option<Candidate>> {
        let o = self.opts;
        let Some(lo) = occ.iter().map(|x| x.2).min() else {
            return Ok(None);
        };
        // Parts grow by the flanks; inward only as far as the shortest spacer
        // allows, so that the gap stays at least 0.
        let inner = o.flank.min(lo / 2);
        let (mut wl, mut wr) = (o.flank + k + inner, inner + k + o.flank);
        // Each site with its flanks, read along the pattern; the first wl and
        // last wr letters are the parts, and what lies between is the spacer.
        let mut sites = clean(
            occ.iter()
                .filter_map(|&(i, start, g, rev)| {
                    site(&self.work.text[i], start, k + g + k, rev, o.flank, o.flank)
                })
                .collect(),
        );
        if sites.len() < 5 {
            return Ok(None);
        }
        let mut motif = with_pieces(&sites, wl, wr, name)?;
        // Refine: rescan and rebuild from the hits, twice.
        for _ in 0..2 {
            let hits = self.scan(&motif, false)?;
            if hits.len() < 5 {
                break;
            }
            sites = self.hit_sites(&hits, |_| true);
            motif = with_pieces(&sites, wl, wr, name)?;
        }
        // Trim weak columns at both ends of each part: the outer ones are
        // dropped, and the inner ones join the spacer.
        let ic = column_ic(&ppm(&motif)?);
        let (a0, a1) = strong_span(&ic[..wl], o.min_ic);
        let (b0, b1) = strong_span(&ic[wl..], o.min_ic);
        if a1 <= a0 || b1 <= b0 {
            return Ok(None);
        }
        sites = clean(
            sites
                .iter()
                .map(|x| x[a0..x.len() - (wr - b1)].to_string())
                .collect(),
        );
        (wl, wr) = (a1 - a0, b1 - b0);
        if sites.len() < 5 {
            return Ok(None);
        }
        motif = with_pieces(&sites, wl, wr, name)?;
        // The gap: scan with the parts over every spacing searched, and keep
        // the run of spacer sizes whose hits outnumber both the negatives and
        // independent placement of the parts.
        let Some((fit, mut ln_spacing)) = self.fit_gap(&motif, wl, wr)? else {
            return Ok(None);
        };
        if ln_spacing > o.stop_pvalue.ln() {
            return Ok(None);
        }
        // Grow each part into informative columns left in the spacer or
        // beyond its outer end.
        motif = self.grow(fit)?;
        if !motif.gap.isgapped {
            // A single part: perhaps the motif is half of a repeat.
            if let Some((rep, ln)) = self.repeat(&motif)? {
                (motif, ln_spacing) = (rep, ln);
            }
        }
        // Motifs without a gap are left to the ungapped finder.
        if !motif.gap.isgapped {
            return Ok(None);
        }
        let (m, hp, sp, sn, ln_p) = self.judge(normalize(&motif)?)?;
        let hits = (0..hp.len())
            .map(|i| {
                (
                    hp.sequence[i] as usize,
                    hp.start[i] as usize,
                    hp.end[i] as usize,
                )
            })
            .collect();
        Ok(Some(Candidate {
            motif: m,
            hits,
            seqs_pos: sp,
            seqs_neg: sn,
            ln_spacing,
            ln_p,
            held_out_hits: Vec::new(),
            k,
            seed: String::new(),
            n_tested: 0,
            tried: 0,
        }))
    }

    /// A motif with its hits in the positives, the positive and negative
    /// sequences with hits, and the natural log of the enrichment P-value.
    fn judge(&self, m: Motif) -> Result<(Motif, HitTable, u64, u64, f64)> {
        let hp = self.scan(&m, false)?;
        let hn = self.scan(&m, true)?;
        let (sp, sn) = (n_unique(&hp.sequence), n_unique(&hn.sequence));
        let (n1, n2) = (self.work.text.len() as u64, self.neg.text.len() as u64);
        let ln_p = ln_fisher_greater(sp, n1 - sp, sn, n2 - sn);
        Ok((m, hp, sp, sn, ln_p))
    }

    /// The motif with the gap its hits support, built from the sites of those
    /// hits, and the natural log of the P-value of its spacing against
    /// independent placement of the two parts.
    fn fit_gap(&self, motif: &Motif, wl: usize, wr: usize) -> Result<Option<(Motif, f64)>> {
        let lo = self.opts.min_gap;
        let Some(prof) = self.profile(motif, wl, wr, &self.work, self.neg)? else {
            return Ok(None);
        };
        let (hp, sp, cp, expect) = (prof.hits, prof.spacers, prof.observed, prof.expected);
        let n_in = cp.len();
        let mut best = 0;
        for j in 1..n_in {
            if cp[j] - expect[j] > cp[best] - expect[best] {
                best = j;
            }
        }
        if cp[best] - expect[best] <= 0.0 {
            return Ok(None);
        }
        // Clearly more hits than expected, by two standard errors.
        let supported = |j: usize| cp[j] > expect[j] + 2.0 * (expect[j] + 1.0).sqrt();
        let (mut a, mut b) = (best, best);
        while a > 0 && supported(a - 1) {
            a -= 1;
        }
        while b + 1 < n_in && supported(b + 1) {
            b += 1;
        }
        let ln_spacing = ln_poisson_upper(
            cp[a..=b].iter().sum::<f64>() as u64,
            expect[a..=b].iter().sum(),
        );
        let sites = self.hit_sites(&hp, |i| (lo + a..=lo + b).contains(&sp[i]));
        if sites.len() < 5 {
            return Ok(None);
        }
        Ok(Some((
            with_pieces(&sites, wl, wr, &motif.name)?,
            ln_spacing,
        )))
    }

    /// Hits of a motif's two parts at every spacing from min_gap to well
    /// beyond max_gap in a set of positives, and for each spacing searched
    /// the hits seen and those expected from the negatives or from
    /// independent placement of the parts, whichever is more.
    fn profile(
        &self,
        motif: &Motif,
        wl: usize,
        wr: usize,
        pos: &SeqSet,
        neg: &SeqSet,
    ) -> Result<Option<Profile>> {
        let (lo, hi) = (self.opts.min_gap, self.opts.max_gap);
        // Spacings beyond those searched show how often the parts meet by
        // chance when they do not depend on each other.
        let ext = hi + 20.max(hi - lo + 1);
        let wide = normalize(&add_gap(motif, Some(&[wl]), &[lo], &[ext])?)?;
        let hp = self.scan_in(&wide, pos)?;
        let hn = self.scan_in(&wide, neg)?;
        if hp.len() < 5 {
            return Ok(None);
        }
        let spacer = |h: &HitTable, i: usize| h.gap_sizes(i).map_or(0, |g| g[0]);
        let sp: Vec<usize> = (0..hp.len()).map(|i| spacer(&hp, i)).collect();
        let n_sizes = ext - lo + 1;
        let mut cp = vec![0.0; n_sizes];
        let mut cn = vec![0.0; n_sizes];
        for &g in &sp {
            cp[g - lo] += 1.0;
        }
        for i in 0..hn.len() {
            cn[spacer(&hn, i) - lo] += 1.0;
        }
        let scale = pos.text.len() as f64 / neg.text.len().max(1) as f64;
        // Under independent placement, hits at each spacing are proportional
        // to the room the sequences leave for it, which the spacings beyond
        // those searched measure.
        let n_in = hi - lo + 1;
        let flat = |seqs: &[String], counts: &[f64]| -> Option<Vec<f64>> {
            let room: Vec<f64> = (lo..=ext)
                .map(|g| {
                    seqs.iter()
                        .map(|s| (s.len() + 1).saturating_sub(wl + wr + g) as f64)
                        .sum()
                })
                .collect();
            let tail_hits: f64 = counts[n_in..].iter().sum();
            let tail_room: f64 = room[n_in..].iter().sum();
            (tail_room > 0.0).then(|| {
                room.iter()
                    .map(|r| r * (tail_hits + 1.0) / tail_room)
                    .collect()
            })
        };
        // Sequences too short for any spacing beyond those searched leave no
        // way to tell a favoured spacing from chance.
        let Some(flat_pos) = flat(&pos.text, &cp) else {
            return Ok(None);
        };
        // Without that room in the negatives, their clustering is left out.
        let flat_neg = flat(&neg.text, &cn).unwrap_or_else(|| cn.clone());
        // Parts made of letters that run together, as in A-rich stretches,
        // meet at short spacings more often than independent placement
        // allows, in the background as well; how much more, shrunk towards
        // none by a few pseudo-hits, raises the expected count.
        let expected: Vec<f64> = (0..n_in)
            .map(|j| {
                let excess = ((cn[j] + 5.0) / (flat_neg[j] + 5.0)).max(1.0);
                (cn[j] * scale).max(flat_pos[j] * excess)
            })
            .collect();
        cp.truncate(n_in);
        Ok(Some(Profile {
            hits: hp,
            spacers: sp,
            observed: cp,
            expected,
        }))
    }

    /// Test a candidate in the sequences held out: the spacing of its hits
    /// against independent placement of its parts, and the enrichment of
    /// sequences with hits, each as the natural log of a P-value.
    fn test(&self, c: &mut Candidate) -> Result<()> {
        let (test_pos, test_neg) = self.held_out.as_ref().unwrap();
        let m = &c.motif;
        let wl = m.gap.gaploc[0];
        let (lo, hi) = (self.opts.min_gap, self.opts.max_gap);
        let range = m.gap.mingap[0].clamp(lo, hi) - lo..=m.gap.maxgap[0].clamp(lo, hi) - lo;
        c.ln_spacing = match self.profile(m, wl, m.width() - wl, test_pos, test_neg)? {
            Some(p) => ln_poisson_upper(
                p.observed[range.clone()].iter().sum::<f64>() as u64,
                p.expected[range].iter().sum(),
            ),
            None => 0.0,
        };
        let hp = self.scan_in(m, test_pos)?;
        let hn = self.scan_in(m, test_neg)?;
        let (sp, sn) = (n_unique(&hp.sequence), n_unique(&hn.sequence));
        let (n1, n2) = (test_pos.text.len() as u64, test_neg.text.len() as u64);
        c.ln_p = ln_fisher_greater(sp, n1 - sp, sn, n2 - sn);
        c.seqs_pos += sp;
        c.seqs_neg += sn;
        c.held_out_hits = (0..hp.len())
            .map(|i| {
                (
                    hp.sequence[i] as usize,
                    hp.start[i] as usize,
                    hp.end[i] as usize,
                )
            })
            .collect();
        Ok(())
    }

    /// Extend the parts of a gapped motif one column at a time, inwards into
    /// the spacer or outwards, while the new column holds at least min_ic
    /// bits.
    fn grow(&self, mut motif: Motif) -> Result<Motif> {
        let min_ic = self.opts.min_ic;
        for _ in 0..8 {
            if !motif.gap.isgapped {
                return Ok(motif);
            }
            let wl = motif.gap.gaploc[0];
            let wr = motif.width() - wl;
            let hits = self.scan(&motif, false)?;
            if hits.len() < 5 {
                return Ok(motif);
            }
            let mut windows = Vec::with_capacity(hits.len());
            for i in 0..hits.len() {
                // One more letter on either side, read along the motif.
                let seq = &self.work.text[hits.sequence[i] as usize];
                let (s, e) = (hits.start[i] as usize, hits.end[i] as usize);
                if s == 0 || e + 1 > seq.len() {
                    continue;
                }
                let w = &seq[s - 1..e + 1];
                windows.push(if hits.minus[i] {
                    rc_string(w)
                } else {
                    w.to_string()
                });
            }
            let windows = clean(windows);
            if windows.len() < 5 {
                return Ok(motif);
            }
            let spacer = windows.iter().map(|w| w.len() - 2 - wl - wr).min().unwrap();
            let ic_at =
                |f: &dyn Fn(&[u8]) -> u8| letters_ic(windows.iter().map(|w| f(w.as_bytes())));
            let left_out = ic_at(&|w| w[0]);
            let right_out = ic_at(&|w| w[w.len() - 1]);
            let (left_in, right_in) = if spacer >= 1 {
                (ic_at(&|w| w[1 + wl]), ic_at(&|w| w[w.len() - 2 - wr]))
            } else {
                (0.0, 0.0)
            };
            let (lo_t, ro_t) = (left_out >= min_ic, right_out >= min_ic);
            let (mut li_t, mut ri_t) = (left_in >= min_ic, right_in >= min_ic);
            if li_t && ri_t && spacer < 2 {
                if right_in < left_in {
                    ri_t = false;
                } else {
                    li_t = false;
                }
            }
            if !(lo_t || ro_t || li_t || ri_t) {
                return Ok(motif);
            }
            let (a, b) = (usize::from(!lo_t), usize::from(!ro_t));
            let sites: Vec<String> = windows
                .iter()
                .map(|w| w[a..w.len() - b].to_string())
                .collect();
            let nl = wl + usize::from(lo_t) + usize::from(li_t);
            let nr = wr + usize::from(ro_t) + usize::from(ri_t);
            motif = with_pieces(&sites, nl, nr, &motif.name)?;
        }
        Ok(motif)
    }

    /// The better of a motif's direct and inverted repeats with the spacings
    /// its hits favour, as the motif and the natural log of its spacing
    /// P-value, if either spacing beats independent placement.
    fn repeat(&self, motif: &Motif) -> Result<Option<(Motif, f64)>> {
        let m = ppm(motif)?;
        let w = m.ncol();
        let mut best: Option<(Motif, f64)> = None;
        for other in [m.clone(), m.reversed()] {
            let mut data = m.as_slice().to_vec();
            data.extend_from_slice(other.as_slice());
            let opts = CreateOptions {
                name: motif.name.clone(),
                ..Default::default()
            };
            let two = from_typed_matrix(
                Matrix::from_vec(4, 2 * w, data)?,
                None,
                Some(Alphabet::Dna),
                MotifType::Ppm,
                &opts,
            )?;
            if let Some((fit, ln)) = self.fit_gap(&two, w, w)?
                && ln <= self.opts.stop_pvalue.ln()
                && best.as_ref().is_none_or(|b| ln < b.1)
            {
                best = Some((fit, ln));
            }
        }
        Ok(best)
    }
}

/// Where a seed pair occurs: sequence, start, spacer and whether it was read
/// on the reverse strand.
type Occurrence = (usize, usize, usize, bool);

/// Every place where `left`, a spacer of `lo` to `hi` letters and `right`
/// occur, with `mismatch` allowing one letter of either part to differ, read
/// along the pattern; the start is on the forward strand.
#[allow(clippy::too_many_arguments)]
fn occurrences(
    seqs: &[String],
    left: &str,
    right: &str,
    lo: usize,
    hi: usize,
    rc: bool,
    mismatch: bool,
    nthreads: usize,
) -> Result<Vec<Occurrence>> {
    let k = left.len();
    let mut forms = vec![(left.to_string(), right.to_string(), false)];
    if rc {
        forms.push((rc_string(right), rc_string(left), true));
    }
    let allowed = usize::from(mismatch);
    // Each sequence on the thread pool, gathered in order.
    let per_seq: Vec<Vec<Occurrence>> = with_threads(nthreads, || {
        seqs.par_iter()
            .enumerate()
            .with_min_len(256)
            .map(|(i, s)| {
                let s = s.as_bytes();
                let mut out = Vec::new();
                for (a, b, rev) in &forms {
                    let (a, b) = (a.as_bytes(), b.as_bytes());
                    for start in 0..s.len().saturating_sub(k - 1) {
                        // Mismatches in the first part; a letter other than A, C, G
                        // and T never matches.
                        let Some(ma) = distance(&s[start..start + k], a, allowed) else {
                            continue;
                        };
                        for g in lo..=hi {
                            let j = start + k + g;
                            if j + k > s.len() {
                                break;
                            }
                            if distance(&s[j..j + k], b, allowed - ma).is_some() {
                                out.push((i, start, g, *rev));
                            }
                        }
                    }
                }
                out
            })
            .collect()
    })?;
    let mut out: Vec<Occurrence> = per_seq.into_iter().flatten().collect();
    if mismatch {
        out.sort_unstable();
        out.dedup();
    }
    Ok(out)
}

/// The number of letters of `seq` differing from `pattern`, if at most
/// `allowed` and all of them A, C, G or T.
fn distance(seq: &[u8], pattern: &[u8], allowed: usize) -> Option<usize> {
    if allowed == 0 {
        // The first letter rules out most places before the rest is compared.
        return (seq.first() == pattern.first() && seq == pattern).then_some(0);
    }
    let mut d = 0;
    for (&x, &p) in seq.iter().zip(pattern) {
        if x != p {
            if !matches!(x, b'A' | b'C' | b'G' | b'T') {
                return None;
            }
            d += 1;
            if d > allowed {
                return None;
            }
        }
    }
    Some(d)
}

/// The letters of a site with flanks, read along the pattern.
fn site(
    seq: &str,
    start: usize,
    width: usize,
    rev: bool,
    flank_l: usize,
    flank_r: usize,
) -> Option<String> {
    // On the reverse strand the pattern's left flank lies to the right.
    let (fa, fb) = if rev {
        (flank_r, flank_l)
    } else {
        (flank_l, flank_r)
    };
    let a = start.checked_sub(fa)?;
    let b = start + width + fb;
    if b > seq.len() {
        return None;
    }
    let s = &seq[a..b];
    Some(if rev { rc_string(s) } else { s.to_string() })
}

/// A motif from sites whose first `left` and last `right` letters are its
/// two parts, with a gap spanning the shortest to the longest spacer between
/// them, or no gap if every spacer is empty.
fn with_pieces(sites: &[String], left: usize, right: usize, name: &str) -> Result<Motif> {
    if let Some(short) = sites.iter().find(|x| x.len() < left + right) {
        return arg_err(format!(
            "the site {short} is shorter than its two parts, {left} and {right} letters"
        ));
    }
    let parts: Vec<String> = sites
        .iter()
        .map(|x| format!("{}{}", &x[..left], &x[x.len() - right..]))
        .collect();
    let strs: Vec<&str> = parts.iter().map(String::as_str).collect();
    let opts = CreateOptions {
        name: name.to_string(),
        ..Default::default()
    };
    let motif = from_strings(&strs, Some(Alphabet::Dna), &opts)?;
    let spacers = sites.iter().map(|x| x.len() - left - right);
    let (lo, hi) = (
        spacers.clone().min().unwrap_or(0),
        spacers.max().unwrap_or(0),
    );
    if hi == 0 {
        return Ok(motif);
    }
    add_gap(&motif, Some(&[left]), &[lo], &[hi])
}

/// The motif's probabilities, with its pseudocount re-applied.
fn ppm(motif: &Motif) -> Result<Matrix> {
    Ok(
        convert_type(&normalize(motif)?, MotifType::Ppm, None, false)?
            .matrix
            .clone(),
    )
}

/// The information content of each column of a DNA PPM.
fn column_ic(ppm: &Matrix) -> Vec<f64> {
    ppm.columns()
        .map(|c| {
            2.0 + c
                .iter()
                .map(|&p| p * p.clamp(1e-12, 1.0).log2())
                .sum::<f64>()
        })
        .collect()
}

/// The information content of one column given as its letters, with a
/// quarter of a count added to each.
fn letters_ic(letters: impl Iterator<Item = u8>) -> f64 {
    let mut counts = [0.25; 4];
    for l in letters {
        if let Some(i) = b"ACGT".iter().position(|&x| x == l) {
            counts[i] += 1.0;
        }
    }
    let total: f64 = counts.iter().sum();
    2.0 + counts
        .iter()
        .map(|c| {
            let p = c / total;
            p * p.log2()
        })
        .sum::<f64>()
}

/// The columns from the first to the last with at least min_ic bits.
fn strong_span(ic: &[f64], min_ic: f64) -> (usize, usize) {
    let first = ic.iter().position(|&x| x >= min_ic);
    let last = ic.iter().rposition(|&x| x >= min_ic);
    match (first, last) {
        (Some(a), Some(b)) => (a, b + 1),
        _ => (0, 0),
    }
}

/// The natural log of P(X >= x) for X Poisson with mean mu.
fn ln_poisson_upper(x: u64, mu: f64) -> f64 {
    // At or below the mean the tail holds about half the distribution or
    // more, and the terms summed below would miss most of it, so it is taken
    // as 1.
    if x == 0 || x as f64 <= mu {
        return 0.0;
    }
    let mu = mu.max(1e-12);
    let n = 200 + (4.0 * (x as f64 + mu).sqrt()) as u64;
    let terms: Vec<f64> = (x..x + n)
        .map(|i| i as f64 * mu.ln() - mu - libm::lgamma(i as f64 + 1.0))
        .collect();
    let top = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (top + terms.iter().map(|t| (t - top).exp()).sum::<f64>().ln()).min(0.0)
}

/// A sequence in upper case, with letters outside ASCII, which never match,
/// as N, so that sites can be cut out by position.
fn dna_text(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii() {
                c.to_ascii_uppercase()
            } else {
                'N'
            }
        })
        .collect()
}

fn n_unique(seqs: &[u32]) -> u64 {
    seqs.iter().collect::<HashSet<_>>().len() as u64
}

/// Sites made only of A, C, G and T.
fn clean(sites: Vec<String>) -> Vec<String> {
    sites
        .into_iter()
        .filter(|x| x.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')))
        .collect()
}

fn kmer_string(code: u64, k: usize) -> String {
    (0..k)
        .map(|i| b"ACGT"[((code >> (2 * (k - 1 - i))) & 3) as usize] as char)
        .collect()
}

fn rc_string(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'a' => 't',
            b'c' => 'g',
            b'g' => 'c',
            b't' => 'a',
            other => other as char,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn occurrences_allow_one_mismatch() {
        let seqs = vec![
            "TTAGGTCAGGAGGTCATT".to_string(),
            "AGGACAGGAGGTCA".to_string(),
        ];
        let exact = occurrences(&seqs, "AGGTCA", "AGGTCA", 1, 5, false, false, 1).unwrap();
        assert_eq!(exact, vec![(0, 2, 2, false)]);
        let one = occurrences(&seqs, "AGGTCA", "AGGTCA", 1, 5, false, true, 1).unwrap();
        assert_eq!(one, vec![(0, 2, 2, false), (1, 0, 2, false)]);
    }

    #[test]
    fn sites_and_pieces() {
        assert_eq!(site("AACGTTT", 2, 3, false, 1, 1).as_deref(), Some("ACGTT"));
        assert_eq!(site("AACGTTT", 2, 3, true, 1, 1).as_deref(), Some("AACGT"));
        assert_eq!(site("AACGTTT", 0, 3, false, 1, 1), None);
        let sites = ["TGACAAGTCA", "TGACAGTCA"].map(String::from);
        let m = with_pieces(&sites, 4, 4, "m").unwrap();
        assert_eq!(
            (m.gap.gaploc[0], m.gap.mingap[0], m.gap.maxgap[0]),
            (4, 1, 2)
        );
    }

    #[test]
    fn poisson_tail() {
        assert!((ln_poisson_upper(2, 1.0) - (1.0 - 2.0 * (-1.0f64).exp()).ln()).abs() < 1e-9);
        assert_eq!(ln_poisson_upper(0, 3.0), 0.0);
        // Counts at or below the mean are no evidence at all.
        assert_eq!(ln_poisson_upper(50, 500.0), 0.0);
        assert_eq!(ln_poisson_upper(10, 1000.0), 0.0);
    }
}
