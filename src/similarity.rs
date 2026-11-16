//! Similarity of contribution weight matrices (CWMs) that keeps their sign
//! and magnitude, with [`compare_cwms`] and [`cwm_similarity_matrix`].
//!
//! [`crate::compare`] compares motifs column by column as probabilities,
//! which for a CWM means taking absolute values and scaling every column to
//! sum to 1. That loses the sign of negative patterns and lets columns
//! carrying almost no contribution count as much as the core. The measures
//! here work on the whole signed matrices instead, as TF-MoDISco and the
//! tools around it compare patterns. As in [`crate::compare`], gapped CWMs
//! are compared with every combination of their gap sizes, the gaps filled
//! by columns of zeros, unless [`CwmCompareOptions::use_gaps`] is off.
//!
//! # References
//!
//! - Shrikumar A, Tian K, Avsec Ž, Shcherbina A, Banerjee A, Sharmin M, Nair S,
//!   Kundaje A (2018). Technical note on transcription factor motif discovery
//!   from importance scores (TF-MoDISco) version 0.5.6.5. arXiv:1811.00416.
//!   <https://doi.org/10.48550/arXiv.1811.00416>

use std::str::FromStr;

use rayon::prelude::*;

use crate::compare::{Comparison, Layout, consensus_range, layout_from_frames, motif_letters};
use crate::error::{Error, Result, arg_err};
use crate::gapfill;
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

/// How two CWMs are scored once aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum CwmMetric {
    /// The cosine of the two matrices, zero-padded where they do not
    /// overlap: the sum of the products of facing values over the overlap,
    /// divided by the product of the two whole matrices' Euclidean norms.
    /// Overhanging columns thus lower the score by their share of either
    /// motif's norm. It ranges from -1 to 1.
    #[default]
    Cosine,
    /// The continuous Jaccard similarity of TF-MoDISco, on matrices scaled
    /// to unit L1 norm and zero-padded: the sum of `min(|a|, |b|)` signed by
    /// whether `a` and `b` agree, divided by the sum of `max(|a|, |b|)`. It
    /// ranges from -1 to 1 and is harsher than the cosine on differences in
    /// magnitude.
    ContinuousJaccard,
}

impl CwmMetric {
    /// The name used for the metric in option strings.
    pub fn as_str(&self) -> &'static str {
        match self {
            CwmMetric::Cosine => "COSINE",
            CwmMetric::ContinuousJaccard => "CJACCARD",
        }
    }

    /// The `min_score` that [`crate::merge::merge_similar`] uses by default:
    /// 0.8 for the cosine and 0.45 for continuous Jaccard, about where
    /// duplicate patterns of one TF-MoDISco run begin.
    pub fn default_min_score(&self) -> f64 {
        match self {
            CwmMetric::Cosine => 0.8,
            CwmMetric::ContinuousJaccard => 0.45,
        }
    }
}

impl std::fmt::Display for CwmMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for CwmMetric {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_uppercase().as_str() {
            "COSINE" => Ok(CwmMetric::Cosine),
            "CJACCARD" => Ok(CwmMetric::ContinuousJaccard),
            _ => arg_err(format!(
                "the CWM metric must be COSINE or CJACCARD, not '{s}'"
            )),
        }
    }
}

/// Options for [`compare_cwms`] and [`cwm_similarity_matrix`].
#[derive(Debug, Clone, PartialEq)]
pub struct CwmCompareOptions {
    /// How aligned matrices are scored; the cosine by default.
    pub metric: CwmMetric,
    /// Compare absolute values, so that a pattern and its opposite, such as
    /// a positive pattern and a negative one with the same letters, count as
    /// alike. Off by default, when they score about 0.
    pub ignore_sign: bool,
    /// Fewest overlapping columns an alignment may have, lowered to the width
    /// of the narrower motif when needed; 5 by default,
    /// where trees use 6.
    pub min_overlap: usize,
    /// Also try the reverse complement of each query, for DNA and RNA.
    pub rc: bool,
    /// Compare gapped CWMs with their gaps, as [`crate::compare`] does; on
    /// by default.
    pub use_gaps: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for CwmCompareOptions {
    fn default() -> Self {
        CwmCompareOptions {
            metric: CwmMetric::Cosine,
            ignore_sign: false,
            min_overlap: 5,
            rc: true,
            use_gaps: true,
            nthreads: 1,
        }
    }
}

/// A motif's columns, made ready for a metric: absolute values if the sign
/// is ignored, and scaled to unit norm, Euclidean for the cosine and L1 for
/// continuous Jaccard.
#[derive(Debug, Clone)]
struct Prepared {
    cols: Vec<Vec<f64>>,
    /// The sum of absolute values of each column, after scaling.
    abs: Vec<f64>,
}

fn prepare(cols: impl Iterator<Item = Vec<f64>>, opts: &CwmCompareOptions) -> Prepared {
    let mut cols: Vec<Vec<f64>> = cols
        .map(|c| {
            if opts.ignore_sign {
                c.iter().map(|v| v.abs()).collect()
            } else {
                c
            }
        })
        .collect();
    let norm = match opts.metric {
        CwmMetric::Cosine => cols.iter().flatten().map(|v| v * v).sum::<f64>().sqrt(),
        CwmMetric::ContinuousJaccard => cols.iter().flatten().map(|v| v.abs()).sum::<f64>(),
    };
    if norm > 0.0 {
        for v in cols.iter_mut().flatten() {
            *v /= norm;
        }
    }
    let abs = cols
        .iter()
        .map(|c| c.iter().map(|v| v.abs()).sum())
        .collect();
    Prepared { cols, abs }
}

/// The best alignment of a target to an oriented query.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CwmAlignment {
    pub score: f64,
    /// Column of the oriented query facing target column 0.
    pub offset: i64,
    pub strand: Strand,
    pub overlap: usize,
}

/// The score of one alignment, with query column `j + offset` facing target
/// column `j` over target columns `lo..hi`.
fn score_at(
    q: &Prepared,
    t: &Prepared,
    offset: i64,
    lo: usize,
    hi: usize,
    metric: CwmMetric,
) -> f64 {
    let qcol = |j: usize| (j as i64 + offset) as usize;
    match metric {
        CwmMetric::Cosine => (lo..hi)
            .map(|j| {
                q.cols[qcol(j)]
                    .iter()
                    .zip(&t.cols[j])
                    .map(|(a, b)| a * b)
                    .sum::<f64>()
            })
            .sum(),
        CwmMetric::ContinuousJaccard => {
            let (mut num, mut max, mut qa, mut ta) = (0.0, 0.0, 0.0, 0.0);
            for j in lo..hi {
                let qc = &q.cols[qcol(j)];
                for (a, b) in qc.iter().zip(&t.cols[j]) {
                    let (aa, bb) = (a.abs(), b.abs());
                    num += aa.min(bb) * a.signum() * b.signum() * f64::from(*a != 0.0 && *b != 0.0);
                    max += aa.max(bb);
                }
                qa += q.abs[qcol(j)];
                ta += t.abs[j];
            }
            let qtotal: f64 = q.abs.iter().sum();
            let ttotal: f64 = t.abs.iter().sum();
            let den = max + (qtotal - qa) + (ttotal - ta);
            if den > 0.0 { num / den } else { 0.0 }
        }
    }
}

/// The best alignment of `t` to `q` and, with `qrc`, to `q`'s reverse
/// complement, over every offset with at least `min_overlap` columns in
/// common. Ties go to the forward strand and then the smaller offset.
fn best_cwm_alignment(
    q: &Prepared,
    qrc: Option<&Prepared>,
    t: &Prepared,
    opts: &CwmCompareOptions,
) -> CwmAlignment {
    let (wq, wt) = (q.cols.len() as i64, t.cols.len() as i64);
    let mo = (opts.min_overlap.max(1) as i64).min(wq).min(wt);
    let mut best = CwmAlignment {
        score: f64::NEG_INFINITY,
        offset: 0,
        strand: Strand::Plus,
        overlap: 0,
    };
    for (strand, oq) in [(Strand::Plus, Some(q)), (Strand::Minus, qrc)] {
        let Some(oq) = oq else { continue };
        for offset in -(wt - mo)..=(wq - mo) {
            let lo = (-offset).max(0) as usize;
            let hi = (wq - offset).min(wt) as usize;
            let s = score_at(oq, t, offset, lo, hi, opts.metric);
            if s > best.score {
                best = CwmAlignment {
                    score: s,
                    offset,
                    strand,
                    overlap: hi - lo,
                };
            }
        }
    }
    best
}

fn check(motifs: &[Motif], opts: &CwmCompareOptions) -> Result<bool> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs given");
    };
    if let Some(m) = motifs.iter().find(|m| m.motif_type != MotifType::Cwm) {
        return arg_err(format!(
            "CWM similarity needs CWMs, but {} is a {}; compare it by PCC instead",
            m.name, m.motif_type
        ));
    }
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if opts.min_overlap == 0 {
        return arg_err("min_overlap must be at least 1");
    }
    Ok(opts.rc && first.alphabet.is_nucleotide())
}

/// Each motif prepared for the metric, and its reverse complement when one
/// is needed.
fn prepare_all(
    motifs: &[Motif],
    opts: &CwmCompareOptions,
    rc: bool,
) -> Vec<(Prepared, Option<Prepared>)> {
    motifs
        .iter()
        .map(|m| {
            let fwd = prepare(m.matrix.columns().map(<[f64]>::to_vec), opts);
            let rev = rc.then(|| prepare(m.matrix.reversed().columns().map(<[f64]>::to_vec), opts));
            (fwd, rev)
        })
        .collect()
}

/// Bins for each query column's products with the pooled columns.
const PVALUE_BINS: usize = 201;

/// A query's per-column score distributions against the pooled columns, on
/// one grid of `PVALUE_BINS` bins from `lo` in steps of `step`.
struct QueryNull {
    lo: f64,
    step: f64,
    /// `pmf[orientation][column]`, orientation 0 forward and 1 reverse
    /// complemented.
    pmf: [Vec<Vec<f64>>; 2],
}

impl QueryNull {
    fn bin(&self, v: f64) -> usize {
        (((v - self.lo) / self.step).round().max(0.0) as usize).min(PVALUE_BINS - 1)
    }
}

/// The distributions of a query's column products with the pooled columns.
/// With `rc`, the pool holds the complement of each of its columns, so the
/// reverse complemented query's products are the forward ones in reverse
/// column order.
fn query_null(q: &Prepared, rc: bool, pool: &[Vec<f64>]) -> QueryNull {
    let products: Vec<Vec<f64>> = q
        .cols
        .iter()
        .map(|col| {
            pool.iter()
                .map(|p| p.iter().zip(col).map(|(a, b)| a * b).sum())
                .collect()
        })
        .collect();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in products.iter().flatten() {
        lo = lo.min(*v);
        hi = hi.max(*v);
    }
    let step = if hi > lo {
        (hi - lo) / (PVALUE_BINS - 1) as f64
    } else {
        1.0
    };
    let mut null = QueryNull {
        lo,
        step,
        pmf: [Vec::new(), Vec::new()],
    };
    let fwd: Vec<Vec<f64>> = products
        .iter()
        .map(|vals| {
            let mut pmf = vec![0.0; PVALUE_BINS];
            for &v in vals {
                pmf[null.bin(v)] += 1.0;
            }
            let inv = 1.0 / vals.len() as f64;
            pmf.iter_mut().for_each(|p| *p *= inv);
            pmf
        })
        .collect();
    let rev = if rc {
        fwd.iter().rev().cloned().collect()
    } else {
        Vec::new()
    };
    null.pmf = [fwd, rev];
    null
}

/// The Bonferroni-corrected P-values of a query's alignments against the
/// pooled null of [`compare_cwms`].
///
/// The null of an alignment is the sum of the distributions of the query
/// columns it covers. Alignments are grouped by orientation and first query
/// column, and from each first column the distributions are convolved one
/// column at a time, so that every alignment's tail is read from the sum it
/// needs without convolving the same columns again. Only one sum is kept at a
/// time, so memory stays small for wide motifs.
fn cosine_pvalues(
    null: &QueryNull,
    q: [&Prepared; 2],
    targets: &[&Prepared],
    aligned: &[CwmAlignment],
    tested: &[usize],
) -> Vec<f64> {
    let mut pvalues = vec![f64::NAN; aligned.len()];
    // Each alignment's orientation, first query column, overlap and
    // observed bin.
    let mut keys: Vec<(usize, usize, usize, usize, usize)> = Vec::new();
    for (t, a) in aligned.iter().enumerate() {
        if a.overlap == 0 {
            continue;
        }
        let o = usize::from(a.strand == Strand::Minus);
        let lo = (-a.offset).max(0) as usize;
        let observed: usize = targets[t]
            .cols
            .iter()
            .enumerate()
            .skip(lo)
            .take(a.overlap)
            .map(|(j, tc)| {
                let qc = &q[o].cols[(j as i64 + a.offset) as usize];
                null.bin(qc.iter().zip(tc).map(|(x, y)| x * y).sum())
            })
            .sum();
        keys.push((o, a.offset.max(0) as usize, a.overlap, observed, t));
    }
    keys.sort_unstable();
    let mut i = 0;
    while i < keys.len() {
        let (o, start, ..) = keys[i];
        let mut pmf = vec![1.0];
        let mut len = 0;
        while i < keys.len() && keys[i].0 == o && keys[i].1 == start {
            let overlap = keys[i].2;
            while len < overlap {
                pmf = crate::compare::convolve(&pmf, &null.pmf[o][start + len]);
                len += 1;
            }
            // The chance of each bin or above, summed from the top so that
            // small tails keep their precision.
            let mut tail = vec![0.0; pmf.len() + 1];
            for b in (0..pmf.len()).rev() {
                tail[b] = tail[b + 1] + pmf[b];
            }
            while i < keys.len() && keys[i].0 == o && keys[i].1 == start && keys[i].2 == overlap {
                let (.., observed, t) = keys[i];
                let p = tail[observed.min(pmf.len())].min(1.0);
                pvalues[t] = (p * tested[t] as f64).min(1.0);
                i += 1;
            }
        }
    }
    pvalues
}

/// The number of alignments tried between two motifs.
fn n_alignments(wq: usize, wt: usize, opts: &CwmCompareOptions, rc: bool) -> usize {
    let mo = opts.min_overlap.max(1).min(wq).min(wt);
    (wq + wt + 1 - 2 * mo) * if rc { 2 } else { 1 }
}

/// Compare each query CWM with every CWM, itself included, keeping their
/// sign and magnitude.
///
/// Returns one [`Comparison`] per query and target pair, ordered by query and
/// then target, as [`crate::compare::compare_motifs`] does. `score` is the
/// similarity of the best alignment. Every motif must be a CWM.
///
/// For the cosine, `pvalue` is the chance that a target of the same width,
/// made of columns drawn at random from all the motifs compared (in either
/// orientation, for DNA and RNA), scores at least as well in the same
/// alignment, Bonferroni-corrected for the number of alignments tried, as
/// [`crate::compare`] does for column correlations. Each column of a
/// target is taken from its motif scaled to unit norm, and the score of an
/// alignment is a sum over its columns, so the chance is found exactly from
/// the distribution of each query column's products with the pooled
/// columns, binned into 201 bins. `qvalue` is its Benjamini-Hochberg
/// adjustment among the comparisons of the same query. Continuous Jaccard is
/// a ratio whose numerator and denominator both change with the columns, so
/// its `pvalue` and `qvalue` are NaN.
///
/// # Example
///
/// ```
/// use universalmotif::Alphabet;
/// use universalmotif::create::{CreateOptions, from_typed_matrix};
/// use universalmotif::matrix::Matrix;
/// use universalmotif::motif::MotifType;
/// use universalmotif::similarity::{CwmCompareOptions, compare_cwms};
///
/// let cwm = |cols: &[[f64; 4]], name: &str| {
///     let opts = CreateOptions { name: name.into(), motif_type: MotifType::Cwm, ..Default::default() };
///     from_typed_matrix(Matrix::from_columns(cols)?, None, Some(Alphabet::Dna), MotifType::Cwm, &opts)
/// };
/// let up = cwm(&[[0.0, 0.0, 0.9, 0.0], [0.8, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.7]], "up")?;
/// let down = cwm(&[[0.0, 0.0, -0.9, 0.0], [-0.8, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, -0.7]], "down")?;
/// let rows = compare_cwms(&[up, down], &[0], &CwmCompareOptions { min_overlap: 3, ..Default::default() })?;
/// assert!((rows[0].score - 1.0).abs() < 1e-12);
/// // The opposite pattern scores no better than 0 in any alignment.
/// assert!(rows[1].score <= 1e-12);
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn compare_cwms(
    motifs: &[Motif],
    queries: &[usize],
    opts: &CwmCompareOptions,
) -> Result<Vec<Comparison>> {
    if opts.use_gaps && gapfill::any_gapped(motifs) {
        if let Some(&q) = queries.iter().find(|&&q| q >= motifs.len()) {
            return arg_err(format!(
                "query index {q} is out of range for {} motifs",
                motifs.len()
            ));
        }
        let exp = gapfill::expand(motifs, &[0.25; 4])?;
        let instances: Vec<usize> = queries
            .iter()
            .flat_map(|&q| exp.of[q].iter().copied())
            .collect();
        let rows = compare_cwms(&exp.motifs, &instances, opts)?;
        // Columns of zeros add nothing to a cosine, so the scores stand.
        return Ok(gapfill::collapse(rows, &exp, queries, false, None));
    }
    let rc = check(motifs, opts)?;
    if let Some(&q) = queries.iter().find(|&&q| q >= motifs.len()) {
        return arg_err(format!(
            "query index {q} is out of range for {} motifs",
            motifs.len()
        ));
    }
    let prep = prepare_all(motifs, opts, rc);
    let letters = motif_letters(motifs)?;
    let n = motifs.len();
    let cosine = opts.metric == CwmMetric::Cosine;
    // The pooled columns, scaled as their motifs are and with their
    // complements when the reverse strand is searched.
    let pool: Vec<Vec<f64>> = if cosine {
        let fwd = prep.iter().flat_map(|p| p.0.cols.iter().cloned());
        let comp = prep
            .iter()
            .flat_map(|p| p.0.cols.iter().map(|c| c.iter().rev().copied().collect()));
        if rc {
            fwd.chain(comp).collect()
        } else {
            fwd.collect()
        }
    } else {
        Vec::new()
    };
    let targets: Vec<&Prepared> = prep.iter().map(|p| &p.0).collect();
    let per_query: Vec<Vec<Comparison>> = with_threads(opts.nthreads, || {
        queries
            .par_iter()
            .map(|&q| {
                let (fwd, rev) = (&prep[q].0, prep[q].1.as_ref());
                let aligned: Vec<CwmAlignment> = (0..n)
                    .into_par_iter()
                    .map(|t| best_cwm_alignment(fwd, rev, &prep[t].0, opts))
                    .collect();
                let mut pvalues = vec![f64::NAN; n];
                let mut qvalues = vec![f64::NAN; n];
                if cosine {
                    let null = query_null(fwd, rc, &pool);
                    let tested: Vec<usize> = (0..n)
                        .map(|t| n_alignments(motifs[q].width(), motifs[t].width(), opts, rc))
                        .collect();
                    let oriented = [fwd, rev.unwrap_or(fwd)];
                    pvalues = cosine_pvalues(&null, oriented, &targets, &aligned, &tested);
                    qvalues = crate::stats::p_adjust_bh_except(&pvalues, Some(q));
                }
                let wq = motifs[q].width() as i64;
                aligned
                    .iter()
                    .enumerate()
                    .map(|(t, a)| {
                        let q_start = a.offset.max(0);
                        // The query's consensus is given in its own orientation.
                        let q_start = if a.strand == Strand::Minus {
                            wq - q_start - a.overlap as i64
                        } else {
                            q_start
                        };
                        Comparison {
                            query: q,
                            target: t,
                            offset: a.offset,
                            strand: a.strand,
                            overlap: a.overlap,
                            score: a.score,
                            pvalue: pvalues[t],
                            qvalue: qvalues[t],
                            query_consensus: consensus_range(&letters[q], q_start, a.overlap),
                            target_consensus: consensus_range(
                                &letters[t],
                                (-a.offset).max(0),
                                a.overlap,
                            ),
                            query_gaps: Vec::new(),
                            target_gaps: Vec::new(),
                        }
                    })
                    .collect()
            })
            .collect()
    })?;
    Ok(per_query.into_iter().flatten().collect())
}

/// The similarity of every pair of CWMs, indexed `[query][target]`, as
/// [`compare_cwms`] scores them. Both measures are symmetric, so the matrix
/// is too.
pub fn cwm_similarity_matrix(motifs: &[Motif], opts: &CwmCompareOptions) -> Result<Vec<Vec<f64>>> {
    let rc = check(motifs, opts)?;
    if motifs.len() < 2 {
        return arg_err("need at least 2 motifs for an all-against-all comparison");
    }
    if opts.use_gaps && gapfill::any_gapped(motifs) {
        let n = motifs.len();
        let queries: Vec<usize> = (0..n).collect();
        let rows = compare_cwms(motifs, &queries, opts)?;
        return Ok(rows
            .chunks(n)
            .map(|r| r.iter().map(|c| c.score).collect())
            .collect());
    }
    let prep = prepare_all(motifs, opts, rc);
    let n = motifs.len();
    with_threads(opts.nthreads, || {
        (0..n)
            .into_par_iter()
            .map(|q| {
                (0..n)
                    .map(|t| {
                        best_cwm_alignment(&prep[q].0, prep[q].1.as_ref(), &prep[t].0, opts).score
                    })
                    .collect()
            })
            .collect()
    })
}

/// The best alignment of motif `t` to the anchor `q`, for merging.
/// Lay CWMs out in a common column frame by aligning each one to
/// `motifs[anchor]` by their similarity, as [`compare_cwms`] does, keeping
/// their signs unless `opts.ignore_sign` is set. Motifs that align better as
/// reverse complements are flipped. Every alignment covers at least
/// `opts.min_overlap` columns, lowered to the width of the narrower motif,
/// so every motif is placed against the anchor, however unlike it.
///
/// # Example
///
/// ```
/// use universalmotif::Alphabet;
/// use universalmotif::create::{CreateOptions, from_typed_matrix};
/// use universalmotif::matrix::Matrix;
/// use universalmotif::motif::MotifType;
/// use universalmotif::similarity::{CwmCompareOptions, cwm_anchor_layout};
///
/// let cwm = |cols: &[[f64; 4]], name: &str| {
///     let opts = CreateOptions { name: name.into(), motif_type: MotifType::Cwm, ..Default::default() };
///     from_typed_matrix(Matrix::from_columns(cols)?, None, Some(Alphabet::Dna), MotifType::Cwm, &opts)
/// };
/// // GAT and, one column further on, ATC: the reverse complement of GAT.
/// let gat = cwm(&[[0.0, 0.0, 0.9, 0.0], [0.8, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.7]], "gat")?;
/// let atc = cwm(&[[0.0, 0.0, 0.0, 0.2], [0.7, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.8], [0.0, 0.9, 0.0, 0.0]], "atc")?;
/// let layout = cwm_anchor_layout(&[gat, atc], 0, &CwmCompareOptions { min_overlap: 3, ..Default::default() })?;
/// assert_eq!(layout.width, 4);
/// assert!(layout.placements[1].rc);
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn cwm_anchor_layout(
    motifs: &[Motif],
    anchor: usize,
    opts: &CwmCompareOptions,
) -> Result<Layout> {
    if anchor >= motifs.len() {
        return arg_err("anchor index out of range");
    }
    if opts.use_gaps && gapfill::any_gapped(motifs) {
        let compare_opts = crate::compare::CompareOptions {
            min_overlap: opts.min_overlap,
            rc: opts.rc,
            nthreads: opts.nthreads,
            ..Default::default()
        };
        return crate::compare::anchor_layout(motifs, anchor, &compare_opts);
    }
    let rc = check(motifs, opts)?;
    let prep = prepare_all(motifs, opts, rc);
    let wa = motifs[anchor].width() as i64;
    let widths: Vec<i64> = motifs.iter().map(|m| m.width() as i64).collect();
    let (fwd, rev) = (&prep[anchor].0, prep[anchor].1.as_ref());
    let aligned: Vec<CwmAlignment> = with_threads(opts.nthreads, || {
        prep.par_iter()
            .map(|(t, _)| best_cwm_alignment(fwd, rev, t, opts))
            .collect()
    })?;
    let mut frame = vec![0i64; motifs.len()];
    let mut flip = vec![false; motifs.len()];
    for (i, a) in aligned.iter().enumerate() {
        if i == anchor || a.overlap == 0 {
            continue;
        }
        // Target column 0 faces column `offset` of the oriented anchor; on the
        // minus strand, the flipped target then starts at the forward
        // anchor's column `wa - wt - offset`.
        if a.strand == Strand::Minus {
            frame[i] = wa - widths[i] - a.offset;
            flip[i] = true;
        } else {
            frame[i] = a.offset;
        }
    }
    Ok(layout_from_frames(&frame, &flip, &widths))
}

pub(crate) fn align_cwm_pair(
    q: &Motif,
    t: &Motif,
    opts: &CwmCompareOptions,
) -> Result<CwmAlignment> {
    let rc = check(&[q.clone(), t.clone()], opts)?;
    let both = prepare_all(&[q.clone(), t.clone()], opts, rc);
    Ok(best_cwm_alignment(
        &both[0].0,
        both[0].1.as_ref(),
        &both[1].0,
        opts,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Alphabet;
    use crate::compare::CompareOptions;
    use crate::create::{CreateOptions, from_typed_matrix};
    use crate::matrix::Matrix;

    fn cwm(cols: &[[f64; 4]], name: &str) -> Motif {
        let opts = CreateOptions {
            name: name.into(),
            motif_type: MotifType::Cwm,
            ..Default::default()
        };
        from_typed_matrix(
            Matrix::from_columns(cols).unwrap(),
            None,
            Some(Alphabet::Dna),
            MotifType::Cwm,
            &opts,
        )
        .unwrap()
    }

    /// Reference values computed by hand for a small pair, the continuous
    /// Jaccard as tfmodisco-lite defines it.
    fn reference(a: &[[f64; 4]], b: &[[f64; 4]], metric: CwmMetric) -> f64 {
        let (la, lb) = (a.len() as i64, b.len() as i64);
        let mut best = f64::NEG_INFINITY;
        let rc: Vec<[f64; 4]> = a.iter().rev().map(|c| [c[3], c[2], c[1], c[0]]).collect();
        for q in [a.to_vec(), rc] {
            for off in -(lb - 1)..la {
                let lo = 0.min(off);
                let hi = la.max(off + lb);
                let w = (hi - lo) as usize;
                let mut pa = vec![[0.0; 4]; w];
                let mut pb = vec![[0.0; 4]; w];
                for (i, c) in q.iter().enumerate() {
                    pa[(i as i64 - lo) as usize] = *c;
                }
                for (j, c) in b.iter().enumerate() {
                    pb[(j as i64 + off - lo) as usize] = *c;
                }
                let s = match metric {
                    CwmMetric::Cosine => {
                        let dot: f64 = pa
                            .iter()
                            .zip(&pb)
                            .map(|(x, y)| (0..4).map(|r| x[r] * y[r]).sum::<f64>())
                            .sum();
                        let na: f64 = pa.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
                        let nb: f64 = pb.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
                        dot / (na * nb)
                    }
                    CwmMetric::ContinuousJaccard => {
                        let l1a: f64 = pa.iter().flatten().map(|v| v.abs()).sum();
                        let l1b: f64 = pb.iter().flatten().map(|v| v.abs()).sum();
                        let (mut num, mut den) = (0.0, 0.0);
                        for (x, y) in pa.iter().flatten().zip(pb.iter().flatten()) {
                            let (x, y) = (x / l1a, y / l1b);
                            num += x.abs().min(y.abs())
                                * x.signum()
                                * y.signum()
                                * f64::from(x != 0.0 && y != 0.0);
                            den += x.abs().max(y.abs());
                        }
                        num / den
                    }
                };
                best = best.max(s);
            }
        }
        best
    }

    #[test]
    fn scores_match_padded_reference() {
        let a = [
            [0.1, -0.2, 0.9, 0.0],
            [0.8, 0.0, -0.1, 0.05],
            [0.0, 0.0, 0.0, 0.7],
            [0.2, 0.1, 0.0, -0.3],
        ];
        let b = [
            [0.0, 0.3, 0.0, 0.0],
            [0.0, 0.0, 0.85, 0.1],
            [0.7, 0.0, 0.0, 0.0],
            [0.0, -0.1, 0.0, 0.6],
            [0.05, 0.0, 0.0, 0.0],
        ];
        let motifs = [cwm(&a, "a"), cwm(&b, "b")];
        for metric in [CwmMetric::Cosine, CwmMetric::ContinuousJaccard] {
            let opts = CwmCompareOptions {
                metric,
                min_overlap: 1,
                ..Default::default()
            };
            let rows = compare_cwms(&motifs, &[0], &opts).unwrap();
            assert!(
                (rows[1].score - reference(&a, &b, metric)).abs() < 1e-12,
                "{metric}"
            );
            assert!((rows[0].score - 1.0).abs() < 1e-12);
            let m = cwm_similarity_matrix(&motifs, &opts).unwrap();
            assert!((m[0][1] - m[1][0]).abs() < 1e-12);
        }
    }

    #[test]
    fn sign_and_strand() {
        let up = cwm(
            &[
                [0.0, 0.0, 0.9, 0.0],
                [0.8, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.7],
            ],
            "up",
        );
        let down = cwm(
            &[
                [0.0, 0.0, -0.9, 0.0],
                [-0.8, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, -0.7],
            ],
            "down",
        );
        // The reverse complement of GAT is ATC.
        let rev = cwm(
            &[
                [0.7, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.8],
                [0.0, 0.9, 0.0, 0.0],
            ],
            "rev",
        );
        let motifs = [up, down, rev];
        let opts = CwmCompareOptions {
            min_overlap: 3,
            ..Default::default()
        };
        let rows = compare_cwms(&motifs, &[0], &opts).unwrap();
        assert!(rows[1].score <= 1e-12);
        assert!((rows[2].score - 1.0).abs() < 1e-12);
        assert_eq!(rows[2].strand, Strand::Minus);
        assert_eq!(rows[2].query_consensus, "GAT");
        let ignore = CwmCompareOptions {
            ignore_sign: true,
            ..opts.clone()
        };
        assert!((compare_cwms(&motifs, &[0], &ignore).unwrap()[1].score - 1.0).abs() < 1e-12);
        let fwd = CwmCompareOptions { rc: false, ..opts };
        assert!(compare_cwms(&motifs, &[0], &fwd).unwrap()[2].score < 0.5);
    }

    #[test]
    fn merging_keeps_cwms_signed() {
        use crate::merge::{
            MergeLinkage, MergeOptions, MergeSimilarOptions, merge_motifs, merge_similar,
        };
        let gat = cwm(
            &[
                [0.0, 0.0, 0.9, 0.0],
                [0.8, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.7],
            ],
            "gat",
        );
        // The same site on the other strand, one column longer.
        let atc = cwm(
            &[
                [0.7, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.8],
                [0.0, 0.9, 0.0, 0.0],
                [0.0, 0.2, 0.0, 0.0],
            ],
            "atc",
        );
        let down = cwm(
            &[
                [0.0, -0.9, 0.0, 0.0],
                [0.0, 0.0, -0.8, 0.0],
                [0.0, -0.7, 0.0, 0.0],
            ],
            "down",
        );
        let opts = MergeOptions {
            min_overlap: 3,
            ..Default::default()
        };
        let merged = merge_motifs(&[gat.clone(), atc.clone()], &opts).unwrap();
        assert_eq!(merged.motif_type, MotifType::Cwm);
        assert_eq!(merged.width(), 4);
        // atc has the larger total contribution, so the merge is in its
        // orientation; the core columns agree and average to themselves, and
        // the extra column comes from atc alone.
        let c = merged.matrix.col(2);
        assert!((c[1] - 0.9).abs() < 1e-12 && c[2].abs() < 1e-12);
        assert!((merged.matrix.col(3)[1] - 0.2).abs() < 1e-12);
        let sopts = MergeSimilarOptions {
            min_overlap: 3,
            ..Default::default()
        };
        let res = merge_similar(&[gat.clone(), atc.clone(), down.clone()], &sopts).unwrap();
        assert_eq!(res.clusters, vec![0, 0, 1]);
        assert_eq!(res.motifs[1].name, "down");
        // A pattern and its opposite join when the sign is ignored.
        let gat_down = cwm(
            &[
                [0.0, 0.0, -0.9, 0.0],
                [-0.8, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, -0.7],
            ],
            "gat_down",
        );
        let signed = merge_similar(&[gat.clone(), gat_down.clone()], &sopts).unwrap();
        assert_eq!(signed.clusters, vec![0, 1]);
        let unsigned = MergeSimilarOptions {
            ignore_sign: true,
            linkage: MergeLinkage::Complete,
            ..sopts
        };
        assert_eq!(
            merge_similar(&[gat, gat_down], &unsigned).unwrap().clusters,
            vec![0, 0]
        );
    }

    #[test]
    fn cosine_pvalues() {
        let a = cwm(
            &[
                [0.0, 0.0, 0.9, 0.0],
                [0.8, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.7],
                [0.1, 0.6, 0.0, 0.0],
                [0.0, 0.0, 0.5, 0.0],
            ],
            "a",
        );
        let b = cwm(
            &[
                [0.0, 0.7, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.8],
                [0.6, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.9, 0.1],
                [0.0, 0.5, 0.0, 0.0],
            ],
            "b",
        );
        let c = cwm(
            &[
                [0.5, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.6],
                [0.0, 0.0, 0.7, 0.0],
                [0.0, 0.8, 0.0, 0.0],
                [0.6, 0.0, 0.0, 0.0],
            ],
            "c",
        );
        let motifs = [a, b, c];
        let rows = compare_cwms(&motifs, &[0], &CwmCompareOptions::default()).unwrap();
        // Itself is the best match, and every P-value is a probability.
        assert!(rows[0].pvalue < rows[1].pvalue && rows[0].pvalue < rows[2].pvalue);
        assert!(
            rows.iter()
                .all(|r| (0.0..=1.0).contains(&r.pvalue) && r.qvalue >= r.pvalue)
        );
        // Continuous Jaccard has scores but no P-values.
        let cj = CwmCompareOptions {
            metric: CwmMetric::ContinuousJaccard,
            ..Default::default()
        };
        let rows = compare_cwms(&motifs, &[0], &cj).unwrap();
        assert!(
            rows.iter()
                .all(|r| r.pvalue.is_nan() && r.qvalue.is_nan() && r.score.is_finite())
        );
    }

    #[test]
    fn needs_cwms() {
        let ppm =
            crate::create::from_strings(&["TGACTCA"], None, &CreateOptions::default()).unwrap();
        let c = cwm(&[[0.0, 0.0, 0.9, 0.0]; 5], "c");
        assert!(compare_cwms(&[c, ppm], &[0], &CwmCompareOptions::default()).is_err());
        assert!("cosine".parse::<CwmMetric>().is_ok() && "PCC".parse::<CwmMetric>().is_err());
    }

    #[test]
    fn cwm_layouts() {
        let one = |l: usize, v: f64| {
            let mut c = [0.0; 4];
            c[l] = v;
            c
        };
        // GATAAG, its reverse complement after two weak columns, and itself
        // after three.
        let fwd: Vec<[f64; 4]> = [(2, 0.9), (0, 0.8), (3, 0.7), (0, 0.6), (0, 0.5), (2, 0.4)]
            .iter()
            .map(|&(l, v)| one(l, v))
            .collect();
        let mut rev = vec![one(0, 0.05); 2];
        rev.extend(fwd.iter().rev().map(|c| [c[3], c[2], c[1], c[0]]));
        let mut shifted = vec![one(3, 0.05); 3];
        shifted.extend(fwd.iter().copied());
        let motifs = [cwm(&fwd, "a"), cwm(&rev, "b"), cwm(&shifted, "c")];
        let layout = cwm_anchor_layout(&motifs, 0, &CwmCompareOptions::default()).unwrap();
        let got: Vec<(usize, bool)> = layout.placements.iter().map(|p| (p.left, p.rc)).collect();
        assert_eq!(got, [(3, false), (3, true), (0, false)]);
        assert_eq!(layout.width, 11);
        // anchor_layout() hands CWMs to it.
        let same = crate::compare::anchor_layout(&motifs, 0, &CompareOptions::default()).unwrap();
        assert_eq!(same.placements, layout.placements);
    }
}
