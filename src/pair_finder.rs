//! Discovery of pairs of motifs that favour a spacing.
//!
//! [`motif_pair_finder`] looks for two motifs, or two copies of one, whose
//! sites keep a fixed or nearly fixed distance apart more often than in a
//! background. The parts come from motifs given beforehand, such as a
//! factor's known motif, and from [`motif_finder`] run on the same
//! sequences, so that a partner no database holds can still be found. Only
//! the spacing is searched, which finds pairs present in a small share of
//! the sequences that [`crate::gapped::gapped_motif_finder`] misses because
//! their parts are too rare or too variable to count exactly.

use std::collections::HashMap;

use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::create::{CreateOptions, from_typed_matrix};
use crate::error::{Result, arg_err};
use crate::finder::{FinderOptions, motif_finder};
use crate::mask::{DustOptions, mask_complexity};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType, Strand};
use crate::ops::{add_gap, normalize};
use crate::packed::EncodedSequences;
use crate::progress;
use crate::scan::{NoOverlaps, ScanOptions, scan_encoded};
use crate::shuffle::{ShuffleMethod, ShuffleOptions, shuffled_backgrounds};
use crate::stats::ln_pbinom_upper;

/// Options for [`motif_pair_finder`].
#[derive(Debug, Clone, PartialEq)]
pub struct PairFinderOptions {
    /// The most motifs [`motif_finder`] adds to the parts; 8 by default.
    pub nmotifs: usize,
    /// The longest spacer searched, at most 31; 31 by default.
    pub max_gap: usize,
    /// The P-value a site of a part must reach; 1e-4 by default. Parts too
    /// short to reach it are scanned at sixteen times the P-value of a
    /// perfect site of their width, up to 1e-3.
    pub hit_pvalue: f64,
    /// The P-value a spacing must reach in the sequences it is chosen in
    /// to be tested in those held out; 1e-3 by default.
    pub train_pvalue: f64,
    /// The P-value a spacing must reach in the sequences held out,
    /// multiplied by the number of spacings tested there; 1e-3 by default.
    pub stop_pvalue: f64,
    /// The fewest pairs of sites a spacing needs where it is chosen; 5 by
    /// default.
    pub min_pairs: u64,
    /// The fewest times more pairs than expected a spacing needs in the
    /// sequences held out; 2 by default.
    pub min_fold: f64,
    /// The most spacers a spacing may span; 3 by default.
    pub max_window: usize,
    /// The share of the sequences and of the background held out to test
    /// the spacings chosen in the rest; 0.3 by default.
    pub holdout: f64,
    /// Pairs whose sites overlap are merged when this share or more of the
    /// sites of one, the one with fewer, overlap those of the other; 0.3 by
    /// default. Above 1, nothing is merged.
    pub merge_overlap: f64,
    /// Mask low-complexity regions with DUST before the search.
    pub mask: bool,
    /// The size of the k-mers kept by the shuffled background, used without
    /// negative sequences; 2 by default.
    pub shuffle_k: usize,
    /// How many shuffles of the sequences make up the background, without
    /// negative sequences; 3 by default.
    pub bkg_shuffles: usize,
    /// Seed for the background, [`motif_finder`] and the sequences held
    /// out.
    pub seed: u64,
    /// Threads to use, with 0 for all.
    pub nthreads: usize,
}

impl Default for PairFinderOptions {
    fn default() -> Self {
        PairFinderOptions {
            nmotifs: 8,
            max_gap: 31,
            hit_pvalue: 1e-4,
            train_pvalue: 1e-3,
            stop_pvalue: 1e-3,
            min_pairs: 5,
            min_fold: 2.0,
            max_window: 3,
            holdout: 0.3,
            merge_overlap: 0.3,
            mask: true,
            shuffle_k: 2,
            bkg_shuffles: 3,
            seed: 1,
            nthreads: 1,
        }
    }
}

/// A pair found by [`motif_pair_finder`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct MotifPair {
    /// The two parts as one motif, gapped from `min_gap` to `max_gap` (or
    /// to 1 when `min_gap` is 0 and the window wider), or ungapped when the
    /// parts sit side by side.
    pub motif: Motif,
    /// The index of the first part in [`PairFinderResult::parts`].
    pub left: usize,
    /// The index of the second part.
    pub right: usize,
    /// The strand of the first part. Pairs are read so that it is the plus
    /// strand where the arrangement allows.
    pub left_strand: Strand,
    /// The strand of the second part.
    pub right_strand: Strand,
    /// `"tandem"`, `"convergent"` or `"divergent"`.
    pub arrangement: &'static str,
    /// The shortest spacer of the spacing.
    pub min_gap: usize,
    /// The longest spacer of the spacing.
    pub max_gap: usize,
    /// Pairs of sites at the spacing in the sequences held out.
    pub pairs: u64,
    /// Those expected from the background.
    pub expected: f64,
    /// `pairs` over `expected`.
    pub fold: f64,
    /// The P-value of the spacing in the sequences held out, multiplied by
    /// the number of spacings tested there and capped at 1.
    pub pvalue: f64,
    /// Base 10 log of `pvalue`.
    pub log10_pvalue: f64,
    /// Whether the pair is a direct repeat of one part whose next rung, the
    /// part's width plus the spacer further on, is enriched as well, as in
    /// a tandem repeat.
    pub repeat: bool,
}

/// What [`motif_pair_finder`] returns.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PairFinderResult {
    /// The parts searched: the motifs given, then those found by
    /// [`motif_finder`], leaving out any whose sites mostly overlap those of
    /// an earlier part.
    pub parts: Vec<Motif>,
    /// The pairs found, most significant first.
    pub pairs: Vec<MotifPair>,
}

/// A pair of parts in one arrangement: the first part, whether it is on
/// the plus strand, the second part and whether it is.
type Config = (u32, bool, u32, bool);

/// The pairs of sites of each arrangement, by spacer.
#[derive(Default)]
struct Counts {
    train: Vec<u64>,
    test: Vec<u64>,
}

/// A spacing chosen in the training sequences.
struct Chosen {
    config: Config,
    lo: usize,
    width: usize,
    ln_train: f64,
}

/// Find pairs of motifs whose sites favour a spacing.
///
/// The parts are the `motifs` given, followed by up to `nmotifs` motifs
/// from [`motif_finder`] run on the same sequences and background; parts
/// whose sites in the sequences mostly overlap those of an earlier part are
/// dropped. Each part is scanned over both sets on both strands, and every
/// pair of sites up to `max_gap` letters apart is counted by arrangement
/// and spacer, a pair and its reverse complement counting as one.
///
/// A share of the sequences and of the background, `holdout`, is set aside.
/// In the rest, each arrangement's window of 1 to `max_window` spacers with
/// the smallest binomial P-value against the background's spacers is
/// chosen; those reaching `train_pvalue` are tested once in the sequences
/// held out, where the P-value multiplied by the number tested must reach
/// `stop_pvalue` and the pairs must be `min_fold` times those expected.
/// Pairs that pass and share most of their sites are merged, keeping the
/// most significant.
pub fn motif_pair_finder<S: AsRef<str> + Sync>(
    sequences: &[S],
    negatives: Option<&[S]>,
    motifs: &[Motif],
    opts: &PairFinderOptions,
) -> Result<PairFinderResult> {
    if sequences.is_empty() {
        return arg_err("no sequences given");
    }
    if motifs.is_empty() && opts.nmotifs == 0 {
        return arg_err("give motifs or set nmotifs above 0, so that there are parts to pair");
    }
    if opts.max_gap > 31 {
        return arg_err("max_gap must be at most 31");
    }
    if opts.max_window == 0 || opts.max_window > opts.max_gap + 1 {
        return arg_err("max_window must be from 1 to max_gap + 1");
    }
    for (name, p) in [
        ("hit_pvalue", opts.hit_pvalue),
        ("train_pvalue", opts.train_pvalue),
        ("stop_pvalue", opts.stop_pvalue),
    ] {
        if !(p > 0.0 && p <= 1.0) {
            return arg_err(format!("{name} must be in (0, 1]"));
        }
    }
    if !(opts.holdout > 0.0 && opts.holdout < 0.9) {
        return arg_err("holdout must be in (0, 0.9)");
    }
    if opts.merge_overlap <= 0.0 {
        return arg_err("merge_overlap must be above 0");
    }
    if let Some(m) = motifs.iter().find(|m| m.alphabet != Alphabet::Dna) {
        return arg_err(format!("{} is not a DNA motif", m.name));
    }
    let progress = progress::current();
    progress::begin(&progress, "finding pairs", 4);
    let pos: Vec<String> = sequences.iter().map(|s| s.as_ref().to_string()).collect();
    let neg: Vec<String> = match negatives {
        Some(n) => n.iter().map(|s| s.as_ref().to_string()).collect(),
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
    progress::check(&progress)?;
    let mut parts: Vec<Motif> = motifs.iter().map(ppm_motif).collect::<Result<_>>()?;
    if opts.nmotifs > 0 {
        let found = progress::without_progress(|| {
            motif_finder(
                &pos,
                Some(&neg),
                &FinderOptions {
                    nmotifs: opts.nmotifs,
                    seed: opts.seed,
                    nthreads: opts.nthreads,
                    ..Default::default()
                },
            )
        })?;
        for f in found {
            parts.push(ppm_motif(&f.motif)?);
        }
    }
    progress::add(&progress, 1);
    progress::check(&progress)?;
    let enc_pos = EncodedSequences::new(&pos, &Alphabet::Dna, opts.nthreads)?;
    let enc_neg = EncodedSequences::new(&neg, &Alphabet::Dna, opts.nthreads)?;
    let mut hits_pos = Vec::with_capacity(parts.len());
    for m in &parts {
        hits_pos.push(sites(m, &enc_pos, opts)?);
    }
    let keep = distinct(&hits_pos);
    let parts: Vec<Motif> = keep.iter().map(|&i| parts[i].clone()).collect();
    let hits_pos: Vec<Vec<Site>> = keep.iter().map(|&i| hits_pos[i].clone()).collect();
    let mut hits_neg = Vec::with_capacity(parts.len());
    for m in &parts {
        hits_neg.push(sites(m, &enc_neg, opts)?);
    }
    progress::add(&progress, 1);
    progress::check(&progress)?;

    let held_pos: Vec<bool> = (0..pos.len())
        .map(|i| held_out(i, opts.seed, opts.holdout))
        .collect();
    let held_neg: Vec<bool> = (0..neg.len())
        .map(|i| held_out(i, opts.seed ^ 0x9e37_79b9_7f4a_7c15, opts.holdout))
        .collect();
    let n_spacers = opts.max_gap + 1;
    let pairs_pos = pair_sites(&hits_pos, pos.len(), opts.max_gap);
    let pairs_neg = pair_sites(&hits_neg, neg.len(), opts.max_gap);
    let count = |pairs: &HashMap<Config, Vec<PairSite>>, held: &[bool]| {
        let mut out: HashMap<Config, Counts> = HashMap::new();
        for (cfg, list) in pairs {
            let c = out.entry(*cfg).or_insert_with(|| Counts {
                train: vec![0; n_spacers],
                test: vec![0; n_spacers],
            });
            for p in list {
                if held[p.seq] {
                    c.test[p.spacer] += 1;
                } else {
                    c.train[p.spacer] += 1;
                }
            }
        }
        out
    };
    let counts_pos = count(&pairs_pos, &held_pos);
    let counts_neg = count(&pairs_neg, &held_neg);
    let zero = vec![0u64; n_spacers];

    // Choose each arrangement's best window where the spacings are learnt.
    let mut configs: Vec<&Config> = counts_pos.keys().collect();
    configs.sort_unstable();
    let ln_train = opts.train_pvalue.ln();
    let mut chosen = Vec::new();
    for cfg in configs {
        let nt = &counts_pos[cfg].train;
        let nb = counts_neg.get(cfg).map_or(&zero, |c| &c.train);
        let mut best: Option<(f64, usize, usize)> = None;
        for w in 1..=opts.max_window {
            for lo in 0..=(n_spacers - w) {
                if nt[lo..lo + w].iter().sum::<u64>() < opts.min_pairs {
                    continue;
                }
                let (ln_p, _, _) = window_test(nt, nb, lo, w);
                if best.is_none_or(|b| ln_p < b.0) {
                    best = Some((ln_p, lo, w));
                }
            }
        }
        if let Some((ln_p, lo, width)) = best
            && ln_p <= ln_train
        {
            chosen.push(Chosen {
                config: *cfg,
                lo,
                width,
                ln_train: ln_p,
            });
        }
    }
    progress::add(&progress, 1);
    progress::check(&progress)?;

    // Test each once in the sequences held out.
    let ln_stop = opts.stop_pvalue.ln();
    let ln_tested = (chosen.len().max(1) as f64).ln();
    let mut calls = Vec::new();
    for c in &chosen {
        let nt = &counts_pos[&c.config].test;
        let nb = counts_neg.get(&c.config).map_or(&zero, |x| &x.test);
        let (ln_p, k, expected) = window_test(nt, nb, c.lo, c.width);
        let ln_adj = (ln_p + ln_tested).min(0.0);
        let fold = if expected > 0.0 {
            k as f64 / expected
        } else {
            0.0
        };
        if ln_adj <= ln_stop && fold >= opts.min_fold {
            calls.push((c, ln_adj, k, expected, fold));
        }
    }
    calls.sort_by(|a, b| {
        a.1.total_cmp(&b.1)
            .then(a.0.ln_train.total_cmp(&b.0.ln_train))
    });

    // Merge pairs that share most of their sites.
    let own = |c: &Chosen| -> Vec<PairSite> {
        pairs_pos[&c.config]
            .iter()
            .filter(|p| (c.lo..c.lo + c.width).contains(&p.spacer))
            .copied()
            .collect()
    };
    let mut kept: Vec<(&Chosen, f64, u64, f64, f64)> = Vec::new();
    let mut kept_sites: Vec<Vec<PairSite>> = Vec::new();
    for call in calls {
        let s = own(call.0);
        if kept_sites
            .iter()
            .any(|k| overlap_share(&s, k) >= opts.merge_overlap)
        {
            continue;
        }
        kept.push(call);
        kept_sites.push(s);
    }

    // Count every sequence and background sequence for the repeat check.
    let all = |pairs: &HashMap<Config, Vec<PairSite>>, cfg: &Config| -> Vec<u64> {
        let mut v = vec![0u64; n_spacers];
        if let Some(list) = pairs.get(cfg) {
            for p in list {
                v[p.spacer] += 1;
            }
        }
        v
    };
    let mut out = Vec::with_capacity(kept.len());
    for (i, (c, ln_adj, k, expected, fold)) in kept.into_iter().enumerate() {
        let (l, lplus, r, rplus) = c.config;
        let repeat = l == r && lplus == rplus && {
            let next = 2 * c.lo + parts[l as usize].width();
            next + c.width <= n_spacers && {
                let (ln_p, kk, e) = window_test(
                    &all(&pairs_pos, &c.config),
                    &all(&pairs_neg, &c.config),
                    next,
                    c.width,
                );
                ln_p <= ln_train && e > 0.0 && kk as f64 >= 2.0 * e
            }
        };
        let hi = c.lo + c.width - 1;
        out.push(MotifPair {
            motif: joined(
                &parts[l as usize],
                lplus,
                &parts[r as usize],
                rplus,
                c.lo,
                hi,
                &format!("pair_{}", i + 1),
            )?,
            left: l as usize,
            right: r as usize,
            left_strand: if lplus { Strand::Plus } else { Strand::Minus },
            right_strand: if rplus { Strand::Plus } else { Strand::Minus },
            arrangement: match (lplus, rplus) {
                (true, true) | (false, false) => "tandem",
                (true, false) => "convergent",
                (false, true) => "divergent",
            },
            min_gap: c.lo,
            max_gap: hi,
            pairs: k,
            expected,
            fold,
            pvalue: ln_adj.exp(),
            log10_pvalue: ln_adj / std::f64::consts::LN_10,
            repeat,
        });
    }
    progress::add(&progress, 1);
    Ok(PairFinderResult { parts, pairs: out })
}

/// A site of a part: its sequence, start, end and whether it is on the
/// plus strand.
#[derive(Clone, Copy)]
struct Site {
    seq: usize,
    start: usize,
    end: usize,
    plus: bool,
}

/// A pair of sites: its sequence, where it starts and ends, and its spacer.
#[derive(Clone, Copy)]
struct PairSite {
    seq: usize,
    from: usize,
    to: usize,
    spacer: usize,
}

/// The motif as a PPM, its pseudocount applied.
fn ppm_motif(m: &Motif) -> Result<Motif> {
    convert_type(&normalize(m)?, MotifType::Ppm, None, false)
}

/// The sites of a part, without overlaps.
fn sites(m: &Motif, seqs: &EncodedSequences, opts: &PairFinderOptions) -> Result<Vec<Site>> {
    let short = 16.0 * 0.25f64.powi(m.width() as i32);
    let scan = ScanOptions {
        pvalue: short.min(1e-3).max(opts.hit_pvalue),
        rc: true,
        no_overlaps: Some(NoOverlaps::default()),
        nthreads: opts.nthreads,
        ..Default::default()
    };
    let h = scan_encoded(std::slice::from_ref(m), seqs, &scan)?;
    Ok((0..h.len())
        .map(|i| Site {
            seq: h.sequence[i] as usize,
            start: h.start[i] as usize,
            end: h.end[i] as usize,
            plus: !h.minus[i],
        })
        .collect())
}

/// The parts to keep, in order, leaving out any half or more of whose sites
/// overlap those of a part already kept.
fn distinct(hits: &[Vec<Site>]) -> Vec<usize> {
    let mut keep = Vec::new();
    let mut taken: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (i, h) in hits.iter().enumerate() {
        if h.is_empty() {
            continue;
        }
        let overlapping = h
            .iter()
            .filter(|s| {
                taken
                    .get(&s.seq)
                    .is_some_and(|v| v.iter().any(|&(a, b)| s.start < b && a < s.end))
            })
            .count();
        if !keep.is_empty() && 2 * overlapping >= h.len() {
            continue;
        }
        keep.push(i);
        for s in h {
            taken.entry(s.seq).or_default().push((s.start, s.end));
        }
    }
    keep
}

/// Every pair of sites up to max_gap letters apart, by arrangement. A pair
/// read on the other strand is the same pair, so each is filed under one
/// reading: the one with the first part on the plus strand where there is
/// one, and otherwise the one with the lower first part.
fn pair_sites(hits: &[Vec<Site>], n_seqs: usize, max_gap: usize) -> HashMap<Config, Vec<PairSite>> {
    let mut by_seq: Vec<Vec<(u32, Site)>> = vec![Vec::new(); n_seqs];
    for (m, h) in hits.iter().enumerate() {
        for s in h {
            by_seq[s.seq].push((m as u32, *s));
        }
    }
    let mut out: HashMap<Config, Vec<PairSite>> = HashMap::new();
    for list in &mut by_seq {
        list.sort_unstable_by_key(|(_, s)| s.start);
        for i in 0..list.len() {
            let (mi, a) = list[i];
            for &(mj, b) in &list[i + 1..] {
                if b.start > a.end + max_gap {
                    break;
                }
                if b.start < a.end {
                    continue;
                }
                let fwd = (mi, a.plus, mj, b.plus);
                let rev = (mj, !b.plus, mi, !a.plus);
                let rank = |c: &Config| (!c.1, c.0, c.2, !c.3);
                let cfg = if rank(&fwd) <= rank(&rev) { fwd } else { rev };
                out.entry(cfg).or_default().push(PairSite {
                    seq: a.seq,
                    from: a.start,
                    to: b.end,
                    spacer: b.start - a.end,
                });
            }
        }
    }
    out
}

/// The natural log of the binomial P-value of the pairs in a window of
/// spacers, against the background's share of pairs there (each spacer
/// given one more pair), with the pairs seen and those expected.
fn window_test(nt: &[u64], nb: &[u64], lo: usize, w: usize) -> (f64, u64, f64) {
    let n: u64 = nt.iter().sum();
    let k: u64 = nt[lo..lo + w].iter().sum();
    let total: f64 = nb.iter().map(|&x| x as f64 + 1.0).sum();
    let q: f64 = nb[lo..lo + w].iter().map(|&x| x as f64 + 1.0).sum::<f64>() / total;
    if k == 0 {
        return (0.0, 0, n as f64 * q);
    }
    (ln_pbinom_upper(k, n, q), k, n as f64 * q)
}

/// The share of the pairs of the smaller set that overlap a pair of the
/// other.
fn overlap_share(a: &[PairSite], b: &[PairSite]) -> f64 {
    let (a, b) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if a.is_empty() {
        return 0.0;
    }
    let mut by_seq: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for p in b {
        by_seq.entry(p.seq).or_default().push((p.from, p.to));
    }
    let hit = a
        .iter()
        .filter(|p| {
            by_seq
                .get(&p.seq)
                .is_some_and(|v| v.iter().any(|&(f, t)| p.from < t && f < p.to))
        })
        .count();
    hit as f64 / a.len() as f64
}

/// Whether sequence i is held out: SplitMix64 of the seed and index, as a
/// number in [0, 1), below the fraction.
fn held_out(i: usize, seed: u64, fraction: f64) -> bool {
    let mut z = seed.wrapping_add((i as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15));
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    ((z >> 11) as f64 / (1u64 << 53) as f64) < fraction
}

/// Two parts, each on its strand, as one motif with a spacer of lo to hi
/// letters between them.
fn joined(
    left: &Motif,
    lplus: bool,
    right: &Motif,
    rplus: bool,
    lo: usize,
    hi: usize,
    name: &str,
) -> Result<Motif> {
    let columns = |m: &Motif, plus: bool| -> Vec<Vec<f64>> {
        let cols: Vec<Vec<f64>> = m.matrix.columns().map(<[f64]>::to_vec).collect();
        if plus {
            cols
        } else {
            cols.into_iter()
                .rev()
                .map(|c| c.into_iter().rev().collect())
                .collect()
        }
    };
    let mut cols = columns(left, lplus);
    let wl = cols.len();
    cols.extend(columns(right, rplus));
    let motif = from_typed_matrix(
        Matrix::from_columns(&cols)?,
        None,
        Some(Alphabet::Dna),
        MotifType::Ppm,
        &CreateOptions {
            name: name.to_string(),
            ..Default::default()
        },
    )?;
    if hi == 0 {
        return Ok(motif);
    }
    add_gap(&motif, Some(&[wl]), &[lo], &[hi])
}
