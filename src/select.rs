//! Choosing and refining motifs after discovery, with [`select_motifs`] and
//! [`refine_motifs`].
//!
//! Discovery tends to return more motifs than a set of sequences needs: rare
//! motifs, and several spellings of the same site that each catch some of
//! its instances. [`select_motifs`] ranks motifs greedily by how many new
//! targets each reaches, sequences or sites within them, so that the motifs
//! adding little can be left out. [`refine_motifs`] rebuilds each motif from
//! its hits, widened and trimmed, so that its edges and letters follow the
//! sites it actually matches.

use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::motif::{Motif, MotifType};
use crate::ops::{TrimFrom, trim_motif};
use crate::scan::{HitTable, NoOverlaps, ScanOptions, scan_table};

/// A stretch of a sequence that motifs should reach, such as a seqlet or a
/// peak summit. Positions are half-open and count from 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Site {
    /// Index of the sequence.
    pub sequence: usize,
    /// First position.
    pub start: usize,
    /// Position just past the site.
    pub end: usize,
}

/// Options for [`select_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct SelectOptions {
    /// P-value cutoff for hits.
    pub pvalue: f64,
    /// The share of all targets a motif must add, when its turn comes, to be
    /// kept; 0.01 by default.
    pub min_gain: f64,
    /// With sites, how many positions a hit must share with a site to reach
    /// it; 1 by default.
    pub min_overlap: usize,
    /// With a similarity matrix, a motif whose similarity to one already
    /// kept reaches this is left out, and its targets are not counted.
    pub max_similarity: Option<f64>,
    /// Also scan the reverse strand of DNA and RNA.
    pub rc: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for SelectOptions {
    fn default() -> Self {
        SelectOptions {
            pvalue: 1e-4,
            min_gain: 0.01,
            min_overlap: 1,
            max_similarity: None,
            rc: true,
            nthreads: 1,
        }
    }
}

/// The result of [`select_motifs`], one entry per motif in the order given.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Selection {
    /// The number of targets: the sites, or the sequences without them.
    pub n_targets: usize,
    /// The targets each motif reaches on its own.
    pub targets: Vec<usize>,
    /// The motif's place in the greedy order, from 1, or `None` for a motif
    /// left out as similar to one kept.
    pub rank: Vec<Option<usize>>,
    /// The share of all targets the motif added when its turn came; 0 for
    /// motifs left out as similar.
    pub gain: Vec<f64>,
    /// The share of all targets reached by the motifs up to and including
    /// this one in the greedy order.
    pub cumulative: Vec<f64>,
    /// Whether the motif is kept: it added at least the minimum gain and is
    /// not too similar to a motif kept before it.
    pub kept: Vec<bool>,
    /// For motifs left out as similar, the index of the kept motif they
    /// resemble.
    pub similar_to: Vec<Option<usize>>,
}

/// Rank motifs by how many new targets each reaches, choosing at each step
/// the motif that adds the most to those reached so far, and keep those
/// that add at least `min_gain` of all targets.
///
/// Targets are `sites` when given, reached by a hit sharing at least
/// `min_overlap` positions with them, and otherwise the sequences, reached
/// by any hit. The gains of motifs chosen later can only be smaller, so the
/// kept motifs are the first ones in the order. With a `similarity` matrix,
/// one row and column per motif, a motif whose similarity to a kept motif
/// reaches `max_similarity` when its turn comes is left out instead, and the
/// targets it reaches are not counted as reached; this removes, say, motifs
/// whose contribution weight matrices repeat another's.
pub fn select_motifs<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    sites: Option<&[Site]>,
    similarity: Option<&[Vec<f64>]>,
    opts: &SelectOptions,
) -> Result<Selection> {
    if !(0.0..=1.0).contains(&opts.min_gain) {
        return arg_err("min_gain must be between 0 and 1");
    }
    if opts.min_overlap == 0 {
        return arg_err("min_overlap must be at least 1");
    }
    if let Some(sim) = similarity
        && (sim.len() != motifs.len() || sim.iter().any(|r| r.len() != motifs.len()))
    {
        return arg_err("the similarity matrix must have one row and one column per motif");
    }
    let n = motifs.len();
    let index = SiteIndex::new(sites, sequences.len())?;
    let scan = ScanOptions {
        pvalue: opts.pvalue,
        rc: opts.rc,
        nthreads: opts.nthreads,
        // Without sites, any hit reaches its sequence, so the first will do.
        first_hit: sites.is_none(),
        ..Default::default()
    };
    let hits = if n == 0 {
        HitTable::default()
    } else {
        scan_table(motifs, sequences, &scan)?
    };
    let words = index.n_targets.div_ceil(64);
    let mut reached = vec![vec![0u64; words]; n];
    for i in 0..hits.len() {
        let m = hits.motif[i] as usize;
        index.reached_by(
            hits.sequence[i] as usize,
            hits.start[i],
            hits.end[i],
            opts.min_overlap,
            |t| {
                reached[m][t / 64] |= 1 << (t % 64);
            },
        );
    }
    let count = |bits: &[u64]| bits.iter().map(|w| w.count_ones() as usize).sum::<usize>();
    let total = index.n_targets.max(1) as f64;
    let mut out = Selection {
        n_targets: index.n_targets,
        targets: reached.iter().map(|r| count(r)).collect(),
        rank: vec![None; n],
        gain: vec![0.0; n],
        cumulative: vec![0.0; n],
        kept: vec![false; n],
        similar_to: vec![None; n],
    };
    let mut covered = vec![0u64; words];
    let mut left: Vec<usize> = (0..n).collect();
    let mut chosen: Vec<usize> = Vec::new();
    let mut rank = 0;
    let mut so_far = 0usize;
    while !left.is_empty() {
        // The motif adding the most, ties going to the earlier motif.
        let (pos, best, new) = left
            .iter()
            .enumerate()
            .map(|(p, &i)| {
                (
                    p,
                    i,
                    reached[i]
                        .iter()
                        .zip(&covered)
                        .map(|(a, c)| (a & !c).count_ones() as usize)
                        .sum(),
                )
            })
            .fold(None, |acc: Option<(usize, usize, usize)>, x| match acc {
                Some(a) if a.2 >= x.2 => Some(a),
                _ => Some(x),
            })
            .unwrap();
        left.swap_remove(pos);
        left.sort_unstable();
        if let (Some(sim), Some(max)) = (similarity, opts.max_similarity)
            && let Some(&k) = chosen.iter().find(|&&k| out.kept[k] && sim[best][k] >= max)
        {
            out.similar_to[best] = Some(k);
            out.cumulative[best] = so_far as f64 / total;
            continue;
        }
        rank += 1;
        so_far += new;
        for (c, a) in covered.iter_mut().zip(&reached[best]) {
            *c |= a;
        }
        out.rank[best] = Some(rank);
        out.gain[best] = new as f64 / total;
        out.cumulative[best] = so_far as f64 / total;
        out.kept[best] = index.n_targets > 0 && out.gain[best] >= opts.min_gain && new > 0;
        chosen.push(best);
    }
    Ok(out)
}

/// The sites of one sequence as (start, end, target index), sorted by start.
type SeqSites = Vec<(usize, usize, usize)>;

/// Sites grouped by sequence and sorted by start, to find those a hit
/// reaches, or every sequence as a target of its own when there are none.
struct SiteIndex {
    n_targets: usize,
    /// For each sequence, its sites, and the largest site length, which
    /// bounds the search.
    by_seq: Option<(Vec<SeqSites>, usize)>,
}

impl SiteIndex {
    fn new(sites: Option<&[Site]>, nseq: usize) -> Result<SiteIndex> {
        let Some(sites) = sites else {
            return Ok(SiteIndex {
                n_targets: nseq,
                by_seq: None,
            });
        };
        let mut by_seq = vec![Vec::new(); nseq];
        let mut longest = 0;
        for (t, s) in sites.iter().enumerate() {
            if s.sequence >= nseq {
                return arg_err(format!(
                    "site {t} is in sequence {}, but there are {nseq} sequences",
                    s.sequence
                ));
            }
            if s.end <= s.start {
                return arg_err(format!("site {t} ends before it starts"));
            }
            by_seq[s.sequence].push((s.start, s.end, t));
            longest = longest.max(s.end - s.start);
        }
        for v in &mut by_seq {
            v.sort_unstable();
        }
        Ok(SiteIndex {
            n_targets: sites.len(),
            by_seq: Some((by_seq, longest)),
        })
    }

    /// Call `f` with each target the hit `start..end` in sequence `seq`
    /// reaches.
    fn reached_by(
        &self,
        seq: usize,
        start: i64,
        end: i64,
        min_overlap: usize,
        mut f: impl FnMut(usize),
    ) {
        let Some((by_seq, longest)) = &self.by_seq else {
            f(seq);
            return;
        };
        let (start, end) = (start.max(0) as usize, end.max(0) as usize);
        let v = &by_seq[seq];
        // Sites starting before start - longest end before the hit starts.
        let from = v.partition_point(|s| s.0 + longest < start);
        for &(s, e, t) in &v[from..] {
            if s >= end {
                break;
            }
            if end.min(e).saturating_sub(start.max(s)) >= min_overlap {
                f(t);
            }
        }
    }
}

/// Options for [`refine_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct RefineOptions {
    /// P-value cutoff for the hits a motif is rebuilt from.
    pub pvalue: f64,
    /// Positions added on either side of each hit before rebuilding, which
    /// trimming then removes again where they hold no information; 5 by
    /// default.
    pub flank: usize,
    /// Rounds of rescanning and rebuilding; 3 by default.
    pub iterations: usize,
    /// Share hits out between motifs: where hits of several motifs overlap,
    /// only the best is used, so that each site builds one motif. On by
    /// default; without it, a weak motif whose hits fall on a stronger
    /// motif's sites drifts into a copy of it.
    pub competitive: bool,
    /// Edge columns below this information content, in bits, are trimmed;
    /// 0.25 by default.
    pub min_ic: f64,
    /// The widest a refined motif may be; its most informative stretch of
    /// this many columns is kept. 30 by default.
    pub max_width: usize,
    /// Motifs with fewer hits than this are left as they are; 20 by default.
    pub min_sites: usize,
    /// With sites, how many positions a hit must share with a site to be
    /// used; 1 by default.
    pub min_overlap: usize,
    /// Also scan the reverse strand of DNA and RNA.
    pub rc: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for RefineOptions {
    fn default() -> Self {
        RefineOptions {
            pvalue: 1e-4,
            flank: 5,
            iterations: 3,
            competitive: true,
            min_ic: 0.25,
            max_width: 30,
            min_sites: 20,
            min_overlap: 1,
            rc: true,
            nthreads: 1,
        }
    }
}

/// A motif from [`refine_motifs`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Refined {
    /// The refined motif, of the type the motif was given as.
    pub motif: Motif,
    /// The hits it was last rebuilt from, or 0 if it was left as it was.
    pub sites: usize,
}

/// Rebuild each motif from its hits, so that it follows the sites it
/// matches.
///
/// In each round the motifs are scanned, keeping only hits that reach the
/// `sites` when they are given, and each motif is rebuilt from its hits
/// widened by `flank` positions on either side, read along the hit's
/// strand, with a pseudocount of the motif's own (or 1 if it has none).
/// Edge columns below `min_ic` are then trimmed, and a motif wider than
/// `max_width` keeps its most informative stretch of that width. Gapped
/// motifs are rebuilt part by part from the letters their hits put in each
/// part, whatever the spacer, with flanks added and trimmed at their outer
/// ends only, so that no part is lost, and their gaps keep their ranges.
/// Motifs with fewer than `min_sites` hits are left as they are, though they
/// still compete for sites. The refined motifs keep the names and other
/// details of those given.
pub fn refine_motifs<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    sites: Option<&[Site]>,
    opts: &RefineOptions,
) -> Result<Vec<Refined>> {
    if opts.max_width == 0 {
        return arg_err("max_width must be at least 1");
    }
    if opts.min_overlap == 0 {
        return arg_err("min_overlap must be at least 1");
    }
    let index = SiteIndex::new(sites, sequences.len())?;
    let mut current: Vec<Motif> = motifs.to_vec();
    let mut used = vec![0usize; motifs.len()];
    if motifs.is_empty() {
        return Ok(Vec::new());
    }
    let scan = ScanOptions {
        pvalue: opts.pvalue,
        rc: opts.rc,
        nthreads: opts.nthreads,
        use_gaps: true,
        no_overlaps: Some(NoOverlaps {
            ignore_motif: opts.competitive,
            ..Default::default()
        }),
        ..Default::default()
    };
    for _ in 0..opts.iterations {
        let hits = scan_table(&current, sequences, &scan)?;
        let mut next = current.clone();
        let mut windows: Vec<Vec<(usize, usize, bool)>> = vec![Vec::new(); current.len()];
        let mut gapped: Vec<Vec<usize>> = vec![Vec::new(); current.len()];
        for i in 0..hits.len() {
            let (seq, start, end) = (hits.sequence[i] as usize, hits.start[i], hits.end[i]);
            if sites.is_some() {
                let mut hit = false;
                index.reached_by(seq, start, end, opts.min_overlap, |_| hit = true);
                if !hit {
                    continue;
                }
            }
            let (a, b) = (start - opts.flank as i64, end + opts.flank as i64);
            if a < 0 || b as usize > sequences[seq].as_ref().len() {
                continue;
            }
            windows[hits.motif[i] as usize].push((seq, a as usize, hits.minus[i]));
            gapped[hits.motif[i] as usize].push(i);
        }
        for (k, m) in current.iter().enumerate() {
            let w = &windows[k];
            if w.len() < opts.min_sites {
                continue;
            }
            let new = if m.gap.isgapped {
                rebuild_gapped(m, &hits, &gapped[k], sequences, opts)?
            } else {
                rebuild(m, sequences, w, m.width() + 2 * opts.flank, opts)?
            };
            if let Some(new) = new {
                next[k] = new;
                used[k] = w.len();
            }
        }
        current = next;
    }
    Ok(current
        .into_iter()
        .zip(used)
        .map(|(motif, sites)| Refined { motif, sites })
        .collect())
}

/// A motif rebuilt from windows of `width` letters, trimmed and capped, in
/// the type of `m`.
fn rebuild<S: AsRef<str>>(
    m: &Motif,
    sequences: &[S],
    windows: &[(usize, usize, bool)],
    width: usize,
    opts: &RefineOptions,
) -> Result<Option<Motif>> {
    let alphabet = &m.alphabet;
    let k = alphabet.len();
    let fold = !matches!(alphabet, Alphabet::Custom(_));
    let lookup: Vec<Option<usize>> = (0u8..=255)
        .map(|b| {
            let c = b as char;
            alphabet.index_of(c).or_else(|| {
                if fold {
                    alphabet.index_of(c.to_ascii_uppercase())
                } else {
                    None
                }
            })
        })
        .collect();
    let mut counts = vec![0.0; k * width];
    for &(seq, start, minus) in windows {
        let s = &sequences[seq].as_ref().as_bytes()[start..start + width];
        for j in 0..width {
            let (b, col) = if minus {
                (s[width - 1 - j], j)
            } else {
                (s[j], j)
            };
            if let Some(mut r) = lookup[b as usize] {
                if minus {
                    if !alphabet.is_nucleotide() {
                        continue;
                    }
                    // Nucleotide rows are in ACGT (or ACGU) order.
                    r = k - 1 - r;
                }
                counts[col * k + r] += 1.0;
            }
        }
    }
    let matrix = crate::matrix::Matrix::from_vec(k, width, counts)?;
    let pseudocount = if m.pseudocount > 0.0 {
        m.pseudocount
    } else {
        1.0
    };
    let mut p = m.clone().into_parts();
    p.matrix = matrix;
    p.motif_type = MotifType::Pcm;
    p.nsites = Some(windows.len() as f64);
    p.pseudocount = pseudocount;
    p.icscore = Motif::compute_icscore(&p.matrix, &p.bkg, MotifType::Pcm, pseudocount, p.nsites);
    p.consensus = Motif::compute_consensus(&p.matrix, &p.alphabet, MotifType::Pcm, pseudocount);
    let pcm = Motif::from_parts(p)?;
    let Some(mut pcm) = trim_motif(&pcm, opts.min_ic, TrimFrom::Both)? else {
        return Ok(None);
    };
    if pcm.width() > opts.max_width {
        pcm = widest_informative(&pcm, opts.max_width)?;
    }
    convert_type(&pcm, m.motif_type, Some(pseudocount), false).map(Some)
}

/// A gapped motif rebuilt from its hits `idx`: each part from the letters
/// its hits put in it, whatever their spacers, with `flank` columns added at
/// the outer ends only and trimmed back as far as `min_ic` asks without
/// emptying a part. The gaps keep their ranges.
fn rebuild_gapped<S: AsRef<str>>(
    m: &Motif,
    hits: &HitTable,
    idx: &[usize],
    sequences: &[S],
    opts: &RefineOptions,
) -> Result<Option<Motif>> {
    let alphabet = &m.alphabet;
    let k = alphabet.len();
    let fold = !matches!(alphabet, Alphabet::Custom(_));
    let index = |b: u8| {
        let c = b as char;
        alphabet.index_of(c).or_else(|| {
            if fold {
                alphabet.index_of(c.to_ascii_uppercase())
            } else {
                None
            }
        })
    };
    let f = opts.flank;
    let core = m.width();
    let width = core + 2 * f;
    let mut counts = vec![0.0; k * width];
    let mut letters = Vec::new();
    let mut n = 0usize;
    for &i in idx {
        if hits.minus[i] && !alphabet.is_nucleotide() {
            continue;
        }
        let seq = sequences[hits.sequence[i] as usize].as_ref().as_bytes();
        let (a, b) = (hits.start[i] as usize, hits.end[i] as usize);
        hits.matched_into(i, sequences, &mut letters);
        let parts: Vec<u8> = letters.iter().copied().filter(|&c| c != b'.').collect();
        if parts.len() != core {
            continue;
        }
        // The flanks read along the hit: before and after it on the +
        // strand, and the complements of after and before it on the -.
        let (left, right): (Vec<Option<usize>>, Vec<Option<usize>>) = if hits.minus[i] {
            let comp = |c: u8| index(c).map(|r| k - 1 - r);
            (
                seq[b..b + f].iter().rev().map(|&c| comp(c)).collect(),
                seq[a - f..a].iter().rev().map(|&c| comp(c)).collect(),
            )
        } else {
            (
                seq[a - f..a].iter().map(|&c| index(c)).collect(),
                seq[b..b + f].iter().map(|&c| index(c)).collect(),
            )
        };
        let row = left
            .into_iter()
            .chain(parts.iter().map(|&c| index(c)))
            .chain(right);
        for (col, r) in row.enumerate() {
            if let Some(r) = r {
                counts[col * k + r] += 1.0;
            }
        }
        n += 1;
    }
    if n < opts.min_sites {
        return Ok(None);
    }
    let pseudocount = if m.pseudocount > 0.0 {
        m.pseudocount
    } else {
        1.0
    };
    let mut p = m.clone().into_parts();
    p.matrix = crate::matrix::Matrix::from_vec(k, width, counts)?;
    p.motif_type = MotifType::Pcm;
    p.nsites = Some(n as f64);
    p.pseudocount = pseudocount;
    p.gap.gaploc = m.gap.gaploc.iter().map(|g| g + f).collect();
    p.icscore = Motif::compute_icscore(&p.matrix, &p.bkg, MotifType::Pcm, pseudocount, p.nsites);
    p.consensus = Motif::compute_consensus(&p.matrix, &p.alphabet, MotifType::Pcm, pseudocount);
    let pcm = Motif::from_parts(p)?;
    // Trim the outer ends, keeping at least one column of the first and
    // last parts.
    let ic: Vec<f64> = pcm
        .matrix
        .columns()
        .map(|c| {
            crate::convert::position_icscore(
                c,
                &pcm.bkg,
                MotifType::Pcm,
                pseudocount,
                n as f64,
                false,
            )
        })
        .collect();
    let first_gap = pcm.gap.gaploc[0];
    let last_gap = *pcm.gap.gaploc.last().unwrap();
    let mut lo = 0;
    while lo + 1 < first_gap && ic[lo] < opts.min_ic {
        lo += 1;
    }
    let mut hi = width;
    while hi > last_gap + 1 && ic[hi - 1] < opts.min_ic {
        hi -= 1;
    }
    let (matrix, gap) = crate::ops::keep_columns(&pcm, lo..hi);
    let trimmed = crate::ops::refresh(&pcm, matrix, gap, MotifType::Pcm)?;
    convert_type(&trimmed, m.motif_type, Some(pseudocount), false).map(Some)
}

/// The stretch of `width` columns with the most information.
fn widest_informative(m: &Motif, width: usize) -> Result<Motif> {
    let nsites = m.nsites.unwrap_or(100.0);
    let ic: Vec<f64> = m
        .matrix
        .columns()
        .map(|c| {
            crate::convert::position_icscore(c, &m.bkg, m.motif_type, m.pseudocount, nsites, false)
        })
        .collect();
    let best = (0..=ic.len() - width)
        .max_by(|&a, &b| {
            let sa: f64 = ic[a..a + width].iter().sum();
            let sb: f64 = ic[b..b + width].iter().sum();
            sa.partial_cmp(&sb)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(b.cmp(&a))
        })
        .unwrap_or(0);
    let mut p = m.clone().into_parts();
    p.matrix = m.matrix.slice_columns(best..best + width);
    p.icscore = Motif::compute_icscore(&p.matrix, &p.bkg, p.motif_type, p.pseudocount, p.nsites);
    p.consensus = Motif::compute_consensus(&p.matrix, &p.alphabet, p.motif_type, p.pseudocount);
    Motif::from_parts(p)
}
