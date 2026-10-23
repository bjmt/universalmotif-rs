//! Motif comparison with [`compare_motifs`].
//!
//! Two motifs are compared by sliding one along the other and summing, over the
//! overlapping columns, the Pearson correlation of each column pair centred on
//! a uniform background; for DNA and RNA the reverse complement is tried as
//! well. The best alignment is kept. Its significance comes from the
//! distribution of that sum under a null model of target columns, computed
//! exactly by convolution, Bonferroni-corrected for the number of alignments
//! tried.
//!
//! Motifs of any alphabet can be compared. Columns are centred on 1/n for an
//! n-letter alphabet, which makes the score the Pearson correlation of
//! probability columns. Reverse complements are only tried, and the parametric
//! null is only available, for DNA and RNA.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::convert::{consensus_letter, consensus_letter_aa, convert_type};
use crate::error::{Error, Result, arg_err};
use crate::matrix::which_max;
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

const PCC_BINS: usize = 51;
const ENUM_GRID_K: usize = 5;
const ENUM_DIRICHLET_N: f64 = 4.0;
const TIE_TOLERANCE: f64 = 1e-10;

/// A motif column of probabilities. Fixed-size arrays of 4 letters for DNA and
/// RNA and of 20 for amino acids are used so that code is compiled for those
/// sizes; other alphabets use vectors.
trait Col: Clone + Send + Sync {
    fn from_slice(v: &[f64]) -> Self;
    fn values(&self) -> &[f64];
}

impl<const N: usize> Col for [f64; N] {
    fn from_slice(v: &[f64]) -> Self {
        let mut a = [0.0; N];
        a.copy_from_slice(v);
        a
    }
    fn values(&self) -> &[f64] {
        self
    }
}

impl Col for Vec<f64> {
    fn from_slice(v: &[f64]) -> Self {
        v.to_vec()
    }
    fn values(&self) -> &[f64] {
        self
    }
}

/// Pearson correlation of two probability columns, which is centred on 1/n, or
/// 0.25 for DNA and RNA.
#[inline(always)]
fn pcc_column<C: Col>(q: &C, t: &C) -> f64 {
    let (q, t) = (q.values(), t.values());
    let centre = 1.0 / q.len() as f64;
    let (mut sqq, mut stt, mut sqt) = (0.0, 0.0, 0.0);
    for i in 0..q.len() {
        let (a, b) = (q[i] - centre, t[i] - centre);
        sqq += a * a;
        stt += b * b;
        sqt += a * b;
    }
    if sqq < 1e-12 || stt < 1e-12 {
        return 0.0;
    }
    sqt / (sqq * stt).sqrt()
}

fn score_to_bin(s: f64) -> usize {
    let t = ((s + 1.0) / 2.0).clamp(0.0, 1.0);
    ((t * (PCC_BINS - 1) as f64 + 0.5).floor() as i64).clamp(0, PCC_BINS as i64 - 1) as usize
}

fn sum_score_to_bin(s: f64, len: usize) -> usize {
    let per_col = ((s / len as f64 + 1.0) / 2.0).clamp(0.0, 1.0);
    let bin = (per_col * (PCC_BINS - 1) as f64 * len as f64 + 0.5).floor() as i64;
    bin.clamp(0, (len * (PCC_BINS - 1)) as i64) as usize
}

/// Reverse complement of DNA or RNA columns. In A, C, G, T order the complement
/// reverses each column.
fn rc_columns<C: Col>(m: &[C]) -> Vec<C> {
    m.iter()
        .rev()
        .map(|c| {
            let mut v = c.values().to_vec();
            v.reverse();
            C::from_slice(&v)
        })
        .collect()
}

/// The 56 columns of the K = 5 simplex grid with Dirichlet-multinomial
/// weights, for the parametric null.
fn enum_columns(bkg: &[f64; 4]) -> Vec<([f64; 4], f64)> {
    let lgamma = libm::lgamma;
    let alpha: Vec<f64> = bkg
        .iter()
        .map(|b| (ENUM_DIRICHLET_N * b).max(0.1))
        .collect();
    let alpha_sum: f64 = alpha.iter().sum();
    let k = ENUM_GRID_K as f64;
    let mut cols = Vec::with_capacity(56);
    for a in 0..=ENUM_GRID_K {
        for b in 0..=ENUM_GRID_K - a {
            for c in 0..=ENUM_GRID_K - a - b {
                let d = ENUM_GRID_K - a - b - c;
                let counts = [a, b, c, d].map(|x| x as f64);
                let mut lw = lgamma(k + 1.0);
                for i in 0..4 {
                    lw -= lgamma(counts[i] + 1.0);
                    lw += lgamma(counts[i] + alpha[i]) - lgamma(alpha[i]);
                }
                lw -= lgamma(k + alpha_sum) - lgamma(alpha_sum);
                cols.push((counts.map(|x| x / k), lw));
            }
        }
    }
    let max = cols.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);
    let mut z = 0.0;
    for c in cols.iter_mut() {
        c.1 = (c.1 - max).exp();
        z += c.1;
    }
    for c in cols.iter_mut() {
        c.1 /= z;
    }
    cols
}

/// The null model column correlations are compared against.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NullModel {
    /// Columns of all the motifs being compared.
    #[default]
    Empirical,
    /// A Dirichlet-multinomial grid of columns around the background.
    Parametric,
}

/// Options for comparisons.
#[derive(Debug, Clone, PartialEq)]
pub struct CompareOptions {
    /// Fewest overlapping columns an alignment may have, lowered to the width
    /// of the narrower motif when needed.
    pub min_overlap: usize,
    /// Also try the reverse complement of each query.
    pub rc: bool,
    /// The null model for P-values; the empirical one by default.
    pub null: NullModel,
    /// Background for the parametric null, in A, C, G, T order.
    pub bkg: [f64; 4],
    /// Score whole queries rather than overlaps, as TOMTOM does by default.
    /// Each query column's correlations are measured against the median of
    /// its correlations with the null columns, so that a typical column
    /// scores zero, and query columns hanging over the ends of the target
    /// count as typical columns. The best alignment is then the one with the
    /// highest sum, and its P-value is that of the best of all the alignments
    /// tried, assuming they are independent, rather than a Bonferroni
    /// correction. This favours alignments covering more of the query and
    /// makes comparisons depend on which motif is the query. Off by default,
    /// when only the overlap counts.
    pub complete: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for CompareOptions {
    fn default() -> Self {
        CompareOptions {
            min_overlap: 5,
            rc: true,
            null: NullModel::Empirical,
            bkg: [0.25; 4],
            complete: false,
            nthreads: 1,
        }
    }
}

/// The best alignment of a target motif to a query motif.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Alignment {
    /// Index of the query motif.
    pub query: usize,
    /// Index of the target motif.
    pub target: usize,
    /// Sum of column correlations over the overlap.
    pub score: f64,
    /// The score the alignment was chosen by: `score`, or with
    /// [`CompareOptions::complete`] the sum over the overlap of each query
    /// column's correlation less its median.
    pub complete_score: f64,
    /// Number of overlapping columns, or 0 if no alignment was possible.
    pub overlap: usize,
    /// Column of the oriented query facing target column 0; can be negative.
    pub offset: i64,
    /// `Minus` if the query's reverse complement aligned best.
    pub strand: Strand,
    /// First overlapping column of the oriented query.
    pub query_start: usize,
    /// Number of alignments tried.
    pub n_tested: usize,
}

/// The alphabet shared by all the motifs.
fn shared_alphabet(motifs: &[Motif]) -> Result<Alphabet> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs to compare");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if first.alphabet.len() < 2 {
        return arg_err("motifs need at least two letters to be compared");
    }
    Ok(first.alphabet.clone())
}

/// Run `$f::<C>` with the column type that suits the alphabet.
macro_rules! by_alphabet {
    ($alphabet:expr, $f:ident($($arg:expr),*)) => {
        match $alphabet.len() {
            4 => $f::<[f64; 4]>($($arg),*),
            20 => $f::<[f64; 20]>($($arg),*),
            _ => $f::<Vec<f64>>($($arg),*),
        }
    };
}

/// PPM columns of the motifs.
fn columns<C: Col>(motifs: &[Motif]) -> Result<Vec<Vec<C>>> {
    motifs
        .iter()
        .map(|m| {
            let ppm = convert_type(m, MotifType::Ppm, None, false)?;
            Ok(ppm.matrix.columns().map(C::from_slice).collect())
        })
        .collect()
}

/// The median correlation of each column of each motif with the null
/// columns, which complete scores subtract so that a typical column scores
/// zero. The null is the columns of all the motifs, in both orientations when
/// reverse complements are tried, or the parametric grid with its weights.
fn column_medians<C: Col>(cols: &[Vec<C>], opts: &CompareOptions, rc_ok: bool) -> Vec<Vec<f64>> {
    let null: Vec<(C, f64)> = match opts.null {
        NullModel::Parametric => enum_columns(&opts.bkg)
            .into_iter()
            .map(|(c, w)| (C::from_slice(&c), w))
            .collect(),
        NullModel::Empirical => {
            let mut v: Vec<(C, f64)> = cols.iter().flatten().map(|c| (c.clone(), 1.0)).collect();
            if rc_ok {
                v.extend(cols.iter().flat_map(|m| rc_columns(m)).map(|c| (c, 1.0)));
            }
            v
        }
    };
    cols.par_iter()
        .map(|m| {
            m.iter()
                .map(|q| {
                    let mut v: Vec<(f64, f64)> =
                        null.iter().map(|(t, w)| (pcc_column(q, t), *w)).collect();
                    v.sort_by(|a, b| a.0.total_cmp(&b.0));
                    let half = v.iter().map(|x| x.1).sum::<f64>() / 2.0;
                    let mut cum = 0.0;
                    v.iter()
                        .find(|x| {
                            cum += x.1;
                            cum >= half
                        })
                        .map_or(0.0, |x| x.0)
                })
                .collect()
        })
        .collect()
}

/// Running sums of some per-column values, so that the sum over any run of
/// columns is a difference of two entries.
fn running_sums(v: impl Iterator<Item = f64>) -> Vec<f64> {
    let mut out = vec![0.0];
    for x in v {
        out.push(out.last().unwrap() + x);
    }
    out
}

/// The best alignment of `t` to the query, as (score it was chosen by,
/// overlap, offset, strand, first query column, alignments tried, sum over
/// the overlap). With `medians`, the running sums of the query's column
/// medians in each orientation, alignments are chosen by complete scores.
#[allow(clippy::type_complexity)]
fn best_alignment<C: Col>(
    q: &[C],
    q_rc: Option<&[C]>,
    t: &[C],
    min_overlap: usize,
    medians: Option<(&[f64], &[f64])>,
) -> (f64, usize, i64, Strand, usize, usize, f64) {
    let (wq, wt) = (q.len() as i64, t.len() as i64);
    let mo = (min_overlap as i64).min(wq).min(wt);
    let mut best = (
        f64::NEG_INFINITY,
        0usize,
        0i64,
        Strand::Plus,
        0usize,
        0usize,
        f64::NEG_INFINITY,
    );
    let orientations: Vec<(&[C], Strand)> = match q_rc {
        Some(r) => vec![(q, Strand::Plus), (r, Strand::Minus)],
        None => vec![(q, Strand::Plus)],
    };
    for (qq, strand) in orientations {
        for d in mo - wt..=wq - mo {
            let t_start = if d >= 0 { 0 } else { -d };
            let t_end = (wq - d).min(wt);
            let len = t_end - t_start;
            if len < mo || len <= 0 {
                continue;
            }
            let q_start = d + t_start;
            let overlap_sum = (0..len).fold(0.0, |acc, j| {
                acc + pcc_column(&qq[(q_start + j) as usize], &t[(t_start + j) as usize])
            });
            let s = match medians {
                Some((plus, minus)) => {
                    let b = if strand == Strand::Plus { plus } else { minus };
                    overlap_sum - (b[(q_start + len) as usize] - b[q_start as usize])
                }
                None => overlap_sum,
            };
            best.5 += 1;
            // Ties (such as both strands of a palindrome) go to the longer
            // overlap, so that comparing the two motifs either way round gives
            // the same mean score, and then to the first alignment, rather than
            // to whichever rounds higher.
            let tied_longer = (s - best.0).abs() <= TIE_TOLERANCE && len as usize > best.1;
            if s > best.0 + TIE_TOLERANCE || tied_longer || best.5 == 1 {
                best = (
                    s,
                    len as usize,
                    d,
                    strand,
                    q_start as usize,
                    best.5,
                    overlap_sum,
                );
            }
        }
    }
    best
}

/// Best alignments for `(query, target)` pairs of motif indices. Reverse
/// complements are only tried for DNA and RNA motifs.
#[doc(hidden)]
pub fn align(
    motifs: &[Motif],
    pairs: &[(usize, usize)],
    opts: &CompareOptions,
) -> Result<Vec<Alignment>> {
    let alphabet = shared_alphabet(motifs)?;
    by_alphabet!(
        alphabet,
        align_with(motifs, pairs, opts, alphabet.is_nucleotide())
    )
}

fn align_with<C: Col>(
    motifs: &[Motif],
    pairs: &[(usize, usize)],
    opts: &CompareOptions,
    nucleotide: bool,
) -> Result<Vec<Alignment>> {
    if opts.min_overlap < 1 {
        return arg_err("min_overlap must be at least 1");
    }
    let cols: Vec<Vec<C>> = columns(motifs)?;
    let rc_ok = opts.rc && nucleotide;
    let rc: Vec<Vec<C>> = if rc_ok {
        cols.iter().map(|c| rc_columns(c)).collect()
    } else {
        Vec::new()
    };
    if pairs
        .iter()
        .any(|&(q, t)| q >= cols.len() || t >= cols.len())
    {
        return arg_err("motif index out of range");
    }
    with_threads(opts.nthreads, || {
        let sums = median_sums(&cols, opts, rc_ok);
        pairs
            .par_iter()
            .map(|&(q, t)| {
                let medians = sums.get(q).map(|(p, m)| (p.as_slice(), m.as_slice()));
                let (complete_score, overlap, offset, strand, query_start, n_tested, score) =
                    best_alignment(
                        &cols[q],
                        rc_ok.then(|| rc[q].as_slice()),
                        &cols[t],
                        opts.min_overlap,
                        medians,
                    );
                Alignment {
                    query: q,
                    target: t,
                    score,
                    complete_score,
                    overlap,
                    offset,
                    strand,
                    query_start,
                    n_tested,
                }
            })
            .collect()
    })
}

/// With complete scores, the running sums of each motif's column medians in
/// its own orientation and reversed, as its reverse complement has them;
/// otherwise nothing.
fn median_sums<C: Col>(
    cols: &[Vec<C>],
    opts: &CompareOptions,
    rc_ok: bool,
) -> Vec<(Vec<f64>, Vec<f64>)> {
    if !opts.complete {
        return Vec::new();
    }
    column_medians(cols, opts, rc_ok)
        .into_iter()
        .map(|m| {
            let plus = running_sums(m.iter().copied());
            let minus = running_sums(m.iter().rev().copied());
            (plus, minus)
        })
        .collect()
}

fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0.0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            out[i + j] += ai * bj;
        }
    }
    out
}

/// The probability of a score in `bin` or above. A score beyond what the null
/// can reach gets the probability of the null's highest score rather than 0,
/// since a null built from a coarse grid of columns may not reach scores that
/// real columns do.
fn tail_prob(pmf: &[f64], bin: usize) -> f64 {
    let Some(top) = pmf.iter().rposition(|&p| p > 0.0) else {
        return 0.0;
    };
    pmf[bin.min(top)..]
        .iter()
        .fold(0.0, |a, &b| a + b)
        .clamp(0.0, 1.0)
}

/// Where a motif sits in a shared column frame, as laid out by
/// [`anchor_layout`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct Placement {
    /// Empty columns to the left of the motif.
    pub left: usize,
    /// Whether the motif is shown as its reverse complement.
    pub rc: bool,
}

/// Motifs laid out in a common column frame by [`anchor_layout`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Layout {
    /// Where each motif sits, in input order.
    pub placements: Vec<Placement>,
    /// The width of the frame, in columns.
    pub width: usize,
}

/// Lay motifs out in a common column frame by aligning each one to
/// `motifs[anchor]`. Motifs that do not overlap the anchor by at least
/// `opts.min_overlap` columns are placed at the anchor's left edge.
pub fn anchor_layout(motifs: &[Motif], anchor: usize, opts: &CompareOptions) -> Result<Layout> {
    if anchor >= motifs.len() {
        return arg_err("anchor index out of range");
    }
    let pairs: Vec<(usize, usize)> = (0..motifs.len()).map(|t| (anchor, t)).collect();
    let al = align(motifs, &pairs, opts)?;
    let widths: Vec<i64> = motifs.iter().map(|m| m.width() as i64).collect();
    let wa = widths[anchor];
    let mut frame = vec![0i64; motifs.len()];
    let mut rc = vec![false; motifs.len()];
    for (i, a) in al.iter().enumerate() {
        if i == anchor || a.overlap < opts.min_overlap {
            continue;
        }
        let l = a.overlap as i64;
        let q_start = a.query_start as i64;
        let t_start = (-a.offset).max(0);
        if a.strand == Strand::Minus {
            frame[i] = (wa - q_start - l) - (widths[i] - t_start - l);
            rc[i] = true;
        } else {
            frame[i] = q_start - t_start;
        }
    }
    let lo = *frame.iter().min().unwrap();
    let hi = frame.iter().zip(&widths).map(|(f, w)| f + w).max().unwrap();
    let placements = frame
        .iter()
        .zip(rc)
        .map(|(&f, rc)| Placement {
            left: (f - lo) as usize,
            rc,
        })
        .collect();
    Ok(Layout {
        placements,
        width: (hi - lo) as usize,
    })
}

/// Bonferroni-corrected P-values of alignments. The parametric null is only
/// available for DNA and RNA motifs.
pub(crate) fn alignment_pvalues(
    motifs: &[Motif],
    alignments: &[Alignment],
    opts: &CompareOptions,
) -> Result<Vec<f64>> {
    let alphabet = shared_alphabet(motifs)?;
    if opts.null == NullModel::Parametric && !alphabet.is_nucleotide() {
        return Err(Error::Unsupported(format!(
            "the parametric null is only available for DNA and RNA motifs; use the empirical null for {alphabet} motifs"
        )));
    }
    by_alphabet!(alphabet, alignment_pvalues_with(motifs, alignments, opts))
}

fn alignment_pvalues_with<C: Col>(
    motifs: &[Motif],
    alignments: &[Alignment],
    opts: &CompareOptions,
) -> Result<Vec<f64>> {
    let cols: Vec<Vec<C>> = columns(motifs)?;
    let flat: Vec<C> = cols.iter().flatten().cloned().collect();
    let grid: Option<Vec<(C, f64)>> = (opts.null == NullModel::Parametric).then(|| {
        enum_columns(&opts.bkg)
            .into_iter()
            .map(|(c, w)| (C::from_slice(&c), w))
            .collect()
    });
    let col_pmf = |q: &C| -> Vec<f64> {
        let mut pmf = vec![0.0; PCC_BINS];
        match &grid {
            None => {
                if flat.is_empty() {
                    return pmf;
                }
                for t in &flat {
                    pmf[score_to_bin(pcc_column(q, t))] += 1.0;
                }
                let inv = 1.0 / flat.len() as f64;
                pmf.iter_mut().for_each(|p| *p *= inv);
            }
            Some(g) => {
                for (c, w) in g {
                    pmf[score_to_bin(pcc_column(q, c))] += w;
                }
            }
        }
        pmf
    };
    // Null distributions for every (start, length) of an oriented query.
    let cache = |q: &[C]| -> Vec<Vec<Vec<f64>>> {
        let pmfs: Vec<Vec<f64>> = q.iter().map(col_pmf).collect();
        (0..q.len())
            .map(|start| {
                let mut row: Vec<Vec<f64>> = vec![Vec::new(), pmfs[start].clone()];
                for len in 2..=q.len() - start {
                    let next = convolve(&row[len - 1], &pmfs[start + len - 1]);
                    row.push(next);
                }
                row
            })
            .collect()
    };

    let mut by_query: Vec<Vec<usize>> = vec![Vec::new(); cols.len()];
    for (k, a) in alignments.iter().enumerate() {
        if a.query >= cols.len() {
            return arg_err("motif index out of range");
        }
        by_query[a.query].push(k);
    }
    let rc_ok = opts.rc && motifs.first().is_some_and(|m| m.alphabet.is_nucleotide());
    let per_query: Vec<Vec<(usize, f64)>> = with_threads(opts.nthreads, || {
        let sums = median_sums(&cols, opts, rc_ok);
        by_query
            .par_iter()
            .enumerate()
            .map(|(q, ks)| {
                if ks.is_empty() {
                    return Vec::new();
                }
                let fwd = cache(&cols[q]);
                let need_rc = (opts.complete && rc_ok)
                    || ks.iter().any(|&k| alignments[k].strand == Strand::Minus);
                let rev = need_rc.then(|| cache(&rc_columns(&cols[q])));
                ks.iter()
                    .map(|&k| {
                        let a = &alignments[k];
                        let mut p = 1.0;
                        if a.overlap > 0 && a.query_start < cols[q].len() && a.n_tested > 0 {
                            if opts.complete {
                                let (plus, minus) = &sums[q];
                                let mut orientations = vec![(&fwd, plus.as_slice())];
                                if let Some(r) = rev.as_ref().filter(|_| rc_ok) {
                                    orientations.push((r, minus.as_slice()));
                                }
                                // An observed score beyond what the null reaches
                                // keeps the probability of the null's highest.
                                let own = if a.strand == Strand::Plus {
                                    Some(&fwd)
                                } else {
                                    rev.as_ref()
                                };
                                let floor = own
                                    .and_then(|c| c.get(a.query_start))
                                    .and_then(|r| r.get(a.overlap))
                                    .map_or(0.0, |pmf| {
                                        tail_prob(pmf, sum_score_to_bin(a.score, a.overlap))
                                    });
                                p = best_of_all(
                                    a.complete_score,
                                    cols[q].len(),
                                    cols[a.target].len(),
                                    opts.min_overlap,
                                    &orientations,
                                )
                                .max(floor);
                            } else {
                                let c = if a.strand == Strand::Plus {
                                    Some(&fwd)
                                } else {
                                    rev.as_ref()
                                };
                                if let Some(pmf) = c
                                    .and_then(|c| c.get(a.query_start))
                                    .and_then(|r| r.get(a.overlap))
                                {
                                    let bin = sum_score_to_bin(a.score, a.overlap);
                                    p = (tail_prob(pmf, bin) * a.n_tested as f64).clamp(0.0, 1.0);
                                }
                            }
                        }
                        (k, p)
                    })
                    .collect()
            })
            .collect()
    })?;
    let mut out = vec![1.0; alignments.len()];
    for (k, p) in per_query.into_iter().flatten() {
        out[k] = p;
    }
    Ok(out)
}

/// Null distributions of an oriented query's column sums, by first column and
/// number of columns.
type NullCache = Vec<Vec<Vec<f64>>>;

/// The probability that the best of all the alignments of a target of width
/// `wt` to a query of width `wq` reaches a complete score of `score`,
/// treating the alignments as independent. `orientations` holds, for each
/// orientation of the query tried, its null distributions by first column
/// and overlap and the running sums of its column medians.
fn best_of_all(
    score: f64,
    wq: usize,
    wt: usize,
    min_overlap: usize,
    orientations: &[(&NullCache, &[f64])],
) -> f64 {
    let (wq, wt) = (wq as i64, wt as i64);
    let mo = (min_overlap as i64).min(wq).min(wt);
    // The log of the probability that no alignment reaches the score.
    let mut log_none = 0.0;
    for (cache, medians) in orientations {
        for d in mo - wt..=wq - mo {
            let t_start = if d >= 0 { 0 } else { -d };
            let len = (wq - d).min(wt) - t_start;
            if len < mo || len <= 0 {
                continue;
            }
            let (start, len) = ((d + t_start) as usize, len as usize);
            // The complete score is the overlap's sum less its medians, so
            // reaching it takes an overlap sum of the score plus the medians.
            let raw = score + medians[start + len] - medians[start];
            // A sum of correlations cannot exceed the number of columns.
            if raw > len as f64 + 1e-9 {
                continue;
            }
            let pmf = &cache[start][len];
            let bin = sum_score_to_bin(raw, len);
            let tail = pmf
                .get(bin..)
                .map_or(0.0, |t| t.iter().sum::<f64>().min(1.0));
            if tail >= 1.0 {
                return 1.0;
            }
            log_none += (-tail).ln_1p();
        }
    }
    (-log_none.exp_m1()).clamp(0.0, 1.0)
}

/// One row of a comparison: a target motif's best alignment to a query.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Comparison {
    /// Index of the query motif.
    pub query: usize,
    /// Index of the target motif.
    pub target: usize,
    /// Column of the oriented query facing the target's first column, which can
    /// be negative.
    pub offset: i64,
    /// `Minus` if the query's reverse complement aligned best.
    pub strand: Strand,
    /// Number of overlapping columns.
    pub overlap: usize,
    /// Mean column correlation over the overlap.
    pub score: f64,
    /// The P-value of the alignment's score.
    pub pvalue: f64,
    /// Benjamini-Hochberg adjusted P-value, among the comparisons of the
    /// same query.
    pub qvalue: f64,
    /// Consensus of the aligned query columns, in the query's own orientation.
    pub query_consensus: String,
    /// Consensus of the aligned target columns.
    pub target_consensus: String,
}

/// Consensus letters of some columns, by the same rules as a motif's own
/// consensus for DNA, RNA and proteins, and the most probable letter for other
/// alphabets.
fn consensus_range<C: Col>(cols: &[C], start: i64, len: usize, alphabet: &Alphabet) -> String {
    let start = start.max(0) as usize;
    let len = len.min(cols.len().saturating_sub(start));
    let letters = alphabet.letters();
    cols[start..start + len]
        .iter()
        .map(|c| {
            let v = c.values();
            match alphabet {
                Alphabet::Dna => consensus_letter(v, false, MotifType::Ppm, 0.0),
                Alphabet::Rna => consensus_letter(v, true, MotifType::Ppm, 0.0),
                Alphabet::Aa => consensus_letter_aa(v, MotifType::Ppm, 0.0),
                _ => letters[which_max(v)],
            }
        })
        .collect()
}

/// Compare each query motif with every motif, itself included.
///
/// Returns one [`Comparison`] per query and target pair, ordered by query and
/// then target, without any filtering.
///
/// # Example
///
/// ```
/// use universalmotif::compare::{CompareOptions, compare_motifs};
/// use universalmotif::create::{CreateOptions, from_strings};
///
/// let opts = CreateOptions { nsites: Some(20.0), pseudocount: 1.0, ..Default::default() };
/// let motifs = [
///     from_strings(&["TGACTCA", "TGAGTCA"], None, &opts)?,
///     from_strings(&["TGACGTCA"], None, &opts)?,
///     from_strings(&["CACGTG"], None, &opts)?,
/// ];
/// // The first motif against all three, best matches first.
/// for c in compare_motifs(&motifs, &[0], &CompareOptions::default())? {
///     println!("{} vs {}: {:.3} (P = {:.2e})", c.query, c.target, c.score, c.pvalue);
/// }
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn compare_motifs(
    motifs: &[Motif],
    queries: &[usize],
    opts: &CompareOptions,
) -> Result<Vec<Comparison>> {
    let n = motifs.len();
    let pairs: Vec<(usize, usize)> = queries
        .iter()
        .flat_map(|&q| (0..n).map(move |t| (q, t)))
        .collect();
    let al = align(motifs, &pairs, opts)?;
    let pvals = alignment_pvalues(motifs, &al, opts)?;
    let alphabet = shared_alphabet(motifs)?;
    by_alphabet!(
        alphabet,
        finish_comparisons(motifs, queries, &al, &pvals, &alphabet)
    )
}

fn finish_comparisons<C: Col>(
    motifs: &[Motif],
    queries: &[usize],
    al: &[Alignment],
    pvals: &[f64],
    alphabet: &Alphabet,
) -> Result<Vec<Comparison>> {
    let cols: Vec<Vec<C>> = columns(motifs)?;

    let mut qvals = vec![f64::NAN; al.len()];
    for &q in queries {
        let rows: Vec<usize> = (0..al.len()).filter(|&k| al[k].query == q).collect();
        let adj = crate::stats::p_adjust_bh(&rows.iter().map(|&k| pvals[k]).collect::<Vec<_>>());
        for (k, v) in rows.into_iter().zip(adj) {
            qvals[k] = v;
        }
    }

    Ok(al
        .iter()
        .enumerate()
        .map(|(k, a)| {
            let wq = cols[a.query].len() as i64;
            let len = a.overlap as i64;
            let q_start = if a.strand == Strand::Minus {
                wq - 1 - (a.query_start as i64 + len - 1)
            } else {
                a.query_start as i64
            };
            let t_start = (-a.offset).max(0);
            Comparison {
                query: a.query,
                target: a.target,
                offset: a.offset,
                strand: a.strand,
                overlap: a.overlap,
                score: if a.overlap > 0 {
                    a.score / a.overlap as f64
                } else {
                    f64::NAN
                },
                pvalue: pvals[k],
                qvalue: qvals[k],
                query_consensus: consensus_range(&cols[a.query], q_start, a.overlap, alphabet),
                target_consensus: consensus_range(&cols[a.target], t_start, a.overlap, alphabet),
            }
        })
        .collect())
}

/// All-against-all comparison results, indexed `[query][target]`.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ComparisonMatrices {
    /// Mean column correlations.
    pub scores: Vec<Vec<f64>>,
    /// P-values of the best alignments.
    pub pvalues: Vec<Vec<f64>>,
    /// Benjamini-Hochberg q-values, adjusted within each query.
    pub qvalues: Vec<Vec<f64>>,
}

/// All-against-all comparison. With `need_pvalues`, P-values and q-values
/// are filled in; otherwise they are NaN.
pub fn compare_matrix(
    motifs: &[Motif],
    opts: &CompareOptions,
    need_pvalues: bool,
) -> Result<ComparisonMatrices> {
    let n = motifs.len();
    if n < 2 {
        return arg_err("need at least 2 motifs for an all-against-all comparison");
    }
    let queries: Vec<usize> = (0..n).collect();
    let (scores, pvals, qvals) = if need_pvalues {
        let rows = compare_motifs(motifs, &queries, opts)?;
        let get = |f: &dyn Fn(&Comparison) -> f64| -> Vec<Vec<f64>> {
            rows.chunks(n).map(|r| r.iter().map(f).collect()).collect()
        };
        (get(&|c| c.score), get(&|c| c.pvalue), get(&|c| c.qvalue))
    } else {
        let pairs: Vec<(usize, usize)> = queries
            .iter()
            .flat_map(|&q| (0..n).map(move |t| (q, t)))
            .collect();
        let al = align(motifs, &pairs, opts)?;
        let s: Vec<Vec<f64>> = al
            .chunks(n)
            .map(|r| {
                r.iter()
                    .map(|a| {
                        if a.overlap > 0 {
                            a.score / a.overlap as f64
                        } else {
                            f64::NAN
                        }
                    })
                    .collect()
            })
            .collect();
        let nan = vec![vec![f64::NAN; n]; n];
        (s, nan.clone(), nan)
    };
    Ok(ComparisonMatrices {
        scores,
        pvalues: pvals,
        qvalues: qvals,
    })
}
