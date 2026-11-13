//! Similarity of contribution weight matrices (CWMs) that keeps their sign
//! and magnitude, with [`compare_cwms`] and [`cwm_similarity_matrix`].
//!
//! [`crate::compare`] compares motifs column by column as probabilities,
//! which for a CWM means taking absolute values and scaling every column to
//! sum to 1. That loses the sign of negative patterns and lets columns
//! carrying almost no contribution count as much as the core. The measures
//! here work on the whole signed matrices instead, as TF-MoDISco and the
//! tools around it compare patterns.

use std::str::FromStr;

use rayon::prelude::*;

use crate::compare::{Comparison, consensus_range, motif_letters};
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

/// How two CWMs are scored once aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
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
    /// of the narrower motif when needed.
    pub min_overlap: usize,
    /// Also try the reverse complement of each query, for DNA and RNA.
    pub rc: bool,
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

/// A query's per-column score distributions against the pooled columns, in
/// each orientation, on one grid of `PVALUE_BINS` bins from `lo` in steps
/// of `step`.
struct QueryNull {
    lo: f64,
    step: f64,
    /// `pmf[orientation][column]`, orientation 0 forward and 1 reversed.
    pmf: [Vec<Vec<f64>>; 2],
}

fn query_null(q: &Prepared, qrc: Option<&Prepared>, pool: &[Vec<f64>]) -> QueryNull {
    let products = |col: &[f64]| -> Vec<f64> {
        pool.iter()
            .map(|p| p.iter().zip(col).map(|(a, b)| a * b).sum())
            .collect()
    };
    let fwd: Vec<Vec<f64>> = q.cols.iter().map(|c| products(c)).collect();
    let rev: Vec<Vec<f64>> = qrc
        .map(|r| r.cols.iter().map(|c| products(c)).collect())
        .unwrap_or_default();
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in fwd.iter().chain(&rev).flatten() {
        lo = lo.min(*v);
        hi = hi.max(*v);
    }
    let step = if hi > lo {
        (hi - lo) / (PVALUE_BINS - 1) as f64
    } else {
        1.0
    };
    let bin = |v: f64| (((v - lo) / step).round().max(0.0) as usize).min(PVALUE_BINS - 1);
    let to_pmf = |vals: &Vec<f64>| {
        let mut pmf = vec![0.0; PVALUE_BINS];
        for &v in vals {
            pmf[bin(v)] += 1.0 / vals.len() as f64;
        }
        pmf
    };
    QueryNull {
        lo,
        step,
        pmf: [
            fwd.iter().map(to_pmf).collect(),
            rev.iter().map(to_pmf).collect(),
        ],
    }
}

/// The Bonferroni-corrected P-value of a cosine alignment against the pooled
/// null of [`compare_cwms`].
///
/// The column distributions are convolved one at a time, keeping only the
/// partial sums that could still end on either side of the observed score:
/// those too low to reach it even with the highest remaining bins are
/// dropped, and those certain to reach it even with the lowest are added to
/// the tail at once. This gives the same tail as the full convolution.
fn cosine_pvalue(
    null: &QueryNull,
    q: &Prepared,
    t: &Prepared,
    a: &CwmAlignment,
    n_tested: usize,
) -> f64 {
    if a.overlap == 0 {
        return f64::NAN;
    }
    let o = usize::from(a.strand == Strand::Minus);
    let lo = (-a.offset).max(0) as usize;
    let cols: Vec<&Vec<f64>> = (lo..lo + a.overlap)
        .map(|j| &null.pmf[o][(j as i64 + a.offset) as usize])
        .collect();
    let mut observed = 0usize;
    for j in lo..lo + a.overlap {
        let qc = (j as i64 + a.offset) as usize;
        let v: f64 = q.cols[qc].iter().zip(&t.cols[j]).map(|(x, y)| x * y).sum();
        observed += (((v - null.lo) / null.step).round().max(0.0) as usize).min(PVALUE_BINS - 1);
    }
    // The lowest and highest bin of each column, and their sums over the
    // columns still to come.
    let range = |c: &Vec<f64>| {
        let first = c.iter().position(|&v| v > 0.0).unwrap_or(0);
        let last = c.iter().rposition(|&v| v > 0.0).unwrap_or(0);
        (first, last)
    };
    let ranges: Vec<(usize, usize)> = cols.iter().map(|c| range(c)).collect();
    let k = cols.len();
    let mut rest_min = vec![0usize; k + 1];
    let mut rest_max = vec![0usize; k + 1];
    for i in (0..k).rev() {
        rest_min[i] = rest_min[i + 1] + ranges[i].0;
        rest_max[i] = rest_max[i + 1] + ranges[i].1;
    }
    let mut tail = 0.0;
    // Partial sums as (offset of index 0, probabilities).
    let mut base = 0usize;
    let mut pmf = vec![1.0];
    for (i, col) in cols.iter().enumerate() {
        let (cmin, cmax) = ranges[i];
        let mut next = vec![0.0; pmf.len() + cmax - cmin];
        for (x, &p) in pmf.iter().enumerate() {
            if p == 0.0 {
                continue;
            }
            for (y, &c) in col[cmin..=cmax].iter().enumerate() {
                next[x + y] += p * c;
            }
        }
        base += cmin;
        // Keep the band whose fate the remaining columns can still decide.
        let (rmin, rmax) = (rest_min[i + 1], rest_max[i + 1]);
        let mut keep_lo = next.len();
        let mut keep_hi = 0;
        for (x, &p) in next.iter().enumerate() {
            if p == 0.0 {
                continue;
            }
            let idx = base + x;
            if idx + rmax < observed {
                continue;
            }
            if idx + rmin >= observed {
                tail += p;
                continue;
            }
            keep_lo = keep_lo.min(x);
            keep_hi = keep_hi.max(x + 1);
        }
        if keep_lo >= keep_hi {
            pmf = Vec::new();
            break;
        }
        let mut band = next[keep_lo..keep_hi].to_vec();
        for (x, v) in band.iter_mut().enumerate() {
            let idx = base + keep_lo + x;
            if idx + rmax < observed || idx + rmin >= observed {
                *v = 0.0;
            }
        }
        base += keep_lo;
        pmf = band;
    }
    // Whatever is left has been through every column.
    tail += pmf
        .iter()
        .enumerate()
        .filter(|(x, _)| base + x >= observed)
        .map(|(_, p)| p)
        .sum::<f64>();
    (tail * n_tested as f64).min(1.0)
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
    let pairs: Vec<(usize, usize)> = queries
        .iter()
        .flat_map(|&q| (0..n).map(move |t| (q, t)))
        .collect();
    let cosine = opts.metric == CwmMetric::Cosine;
    let nulls: Vec<Option<QueryNull>> = if cosine {
        // The pooled columns, scaled as their motifs are and with their
        // complements when the reverse strand is searched.
        let mut pool: Vec<Vec<f64>> = prep.iter().flat_map(|p| p.0.cols.iter().cloned()).collect();
        if rc {
            let comp: Vec<Vec<f64>> = pool
                .iter()
                .map(|c| c.iter().rev().copied().collect())
                .collect();
            pool.extend(comp);
        }
        let mut wanted = vec![false; n];
        for &q in queries {
            wanted[q] = true;
        }
        with_threads(opts.nthreads, || {
            (0..n)
                .into_par_iter()
                .map(|q| wanted[q].then(|| query_null(&prep[q].0, prep[q].1.as_ref(), &pool)))
                .collect()
        })?
    } else {
        (0..n).map(|_| None).collect()
    };
    let mut rows: Vec<Comparison> = with_threads(opts.nthreads, || {
        pairs
            .par_iter()
            .map(|&(q, t)| {
                let a = best_cwm_alignment(&prep[q].0, prep[q].1.as_ref(), &prep[t].0, opts);
                let pvalue = match &nulls[q] {
                    Some(null) => {
                        let oq = if a.strand == Strand::Minus {
                            prep[q].1.as_ref().unwrap()
                        } else {
                            &prep[q].0
                        };
                        let tested = n_alignments(motifs[q].width(), motifs[t].width(), opts, rc);
                        cosine_pvalue(null, oq, &prep[t].0, &a, tested)
                    }
                    None => f64::NAN,
                };
                let wq = motifs[q].width() as i64;
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
                    pvalue,
                    qvalue: f64::NAN,
                    query_consensus: consensus_range(&letters[q], q_start, a.overlap),
                    target_consensus: consensus_range(&letters[t], (-a.offset).max(0), a.overlap),
                }
            })
            .collect()
    })?;
    if cosine {
        for chunk in rows.chunks_mut(n) {
            let adj =
                crate::stats::p_adjust_bh(&chunk.iter().map(|c| c.pvalue).collect::<Vec<_>>());
            for (c, v) in chunk.iter_mut().zip(adj) {
                c.qvalue = v;
            }
        }
    }
    Ok(rows)
}

/// The similarity of every pair of CWMs, indexed `[query][target]`, as
/// [`compare_cwms`] scores them. Both measures are symmetric, so the matrix
/// is too.
pub fn cwm_similarity_matrix(motifs: &[Motif], opts: &CwmCompareOptions) -> Result<Vec<Vec<f64>>> {
    let rc = check(motifs, opts)?;
    if motifs.len() < 2 {
        return arg_err("need at least 2 motifs for an all-against-all comparison");
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

    /// The scores of the numpy versions, which match tfmodisco-lite's
    /// continuous Jaccard, written out by hand for a small pair.
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
}
