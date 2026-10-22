//! Motif comparison (`compare_motifs_lite()`).
//!
//! Two DNA or RNA motifs are compared by sliding one along the other (and
//! its reverse complement) and summing, over the overlapping columns, the
//! Pearson correlation of each column pair centred on a uniform background.
//! The best alignment is kept. Its significance comes from the distribution
//! of that sum under a null model of target columns, computed exactly by
//! convolution, Bonferroni-corrected for the number of alignments tried.
//! This follows universalmotif's `src/compare_motifs_lite.cpp`, itself a port
//! of `yamcmp` from yamtk.

use rayon::prelude::*;

use crate::convert::convert_type;
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

const PCC_BINS: usize = 51;
const ENUM_GRID_K: usize = 5;
const ENUM_DIRICHLET_N: f64 = 4.0;
const CONSENSUS_THRESHOLD: f64 = 0.75;

type Column = [f64; 4];

/// Pearson correlation of two columns, centred on 0.25.
fn pcc_column(q: &Column, t: &Column) -> f64 {
    let (mut sqq, mut stt, mut sqt) = (0.0, 0.0, 0.0);
    for i in 0..4 {
        let (a, b) = (q[i] - 0.25, t[i] - 0.25);
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

fn iupac(mask: u8) -> char {
    match mask {
        0x1 => 'A',
        0x2 => 'C',
        0x4 => 'G',
        0x8 => 'T',
        0x3 => 'M',
        0x5 => 'R',
        0x9 => 'W',
        0x6 => 'S',
        0xA => 'Y',
        0xC => 'K',
        0x7 => 'V',
        0xB => 'H',
        0xD => 'D',
        0xE => 'B',
        _ => 'N',
    }
}

/// IUPAC letter covering the most probable letters until 75% is reached.
fn consensus_char(p: &Column) -> char {
    let mut order = [0usize, 1, 2, 3];
    for i in 1..4 {
        let mut j = i;
        while j > 0 && p[order[j]] > p[order[j - 1]] {
            order.swap(j, j - 1);
            j -= 1;
        }
    }
    let (mut cum, mut mask) = (0.0, 0u8);
    for &k in &order {
        cum += p[k];
        mask |= 1 << k;
        if cum >= CONSENSUS_THRESHOLD {
            break;
        }
    }
    iupac(mask)
}

fn rc_columns(m: &[Column]) -> Vec<Column> {
    m.iter().rev().map(|c| [c[3], c[2], c[1], c[0]]).collect()
}

/// The 56 columns of the K = 5 simplex grid with Dirichlet-multinomial
/// weights, for the parametric null.
fn enum_columns(bkg: &Column) -> Vec<(Column, f64)> {
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
    /// Fewest overlapping columns an alignment may have (lowered to the
    /// width of the narrower motif when needed).
    pub min_overlap: usize,
    /// Also try the reverse complement of each query.
    pub rc: bool,
    pub null: NullModel,
    /// Background (A, C, G, T) for the parametric null.
    pub bkg: [f64; 4],
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
            nthreads: 1,
        }
    }
}

/// The best alignment of a target motif to a query motif.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    pub query: usize,
    pub target: usize,
    /// Sum of column correlations over the overlap.
    pub score: f64,
    /// Number of overlapping columns (0 if no alignment was possible).
    pub overlap: usize,
    /// Query column (in the oriented query) facing target column 0; can be
    /// negative.
    pub offset: i64,
    /// `Minus` if the query's reverse complement aligned best.
    pub strand: Strand,
    /// First overlapping column of the oriented query.
    pub query_start: usize,
    /// Number of alignments tried.
    pub n_tested: usize,
}

/// PPM columns of DNA/RNA motifs.
fn columns(motifs: &[Motif]) -> Result<Vec<Vec<Column>>> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs to compare");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if !first.alphabet.is_nucleotide() {
        return Err(Error::Unsupported(format!(
            "only DNA and RNA motifs can be compared; got {}",
            first.alphabet
        )));
    }
    motifs
        .iter()
        .map(|m| {
            let ppm = convert_type(m, MotifType::Ppm, None, false)?;
            Ok(ppm
                .matrix
                .columns()
                .map(|c| [c[0], c[1], c[2], c[3]])
                .collect())
        })
        .collect()
}

fn best_alignment(
    q: &[Column],
    q_rc: Option<&[Column]>,
    t: &[Column],
    min_overlap: usize,
) -> (f64, usize, i64, Strand, usize, usize) {
    let (wq, wt) = (q.len() as i64, t.len() as i64);
    let mo = (min_overlap as i64).min(wq).min(wt);
    let mut best = (
        f64::NEG_INFINITY,
        0usize,
        0i64,
        Strand::Plus,
        0usize,
        0usize,
    );
    let orientations: Vec<(&[Column], Strand)> = match q_rc {
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
            let s = (0..len).fold(0.0, |acc, j| {
                acc + pcc_column(&qq[(q_start + j) as usize], &t[(t_start + j) as usize])
            });
            best.5 += 1;
            if s > best.0 || best.5 == 1 {
                best = (s, len as usize, d, strand, q_start as usize, best.5);
            }
        }
    }
    best
}

/// Best alignments for `(query, target)` pairs of motif indices.
pub fn align(
    motifs: &[Motif],
    pairs: &[(usize, usize)],
    opts: &CompareOptions,
) -> Result<Vec<Alignment>> {
    if opts.min_overlap < 1 {
        return arg_err("min_overlap must be at least 1");
    }
    let cols = columns(motifs)?;
    let rc: Vec<Vec<Column>> = cols.iter().map(|c| rc_columns(c)).collect();
    if pairs
        .iter()
        .any(|&(q, t)| q >= cols.len() || t >= cols.len())
    {
        return arg_err("motif index out of range");
    }
    with_threads(opts.nthreads, || {
        pairs
            .par_iter()
            .map(|&(q, t)| {
                let (score, overlap, offset, strand, query_start, n_tested) = best_alignment(
                    &cols[q],
                    opts.rc.then_some(rc[q].as_slice()),
                    &cols[t],
                    opts.min_overlap,
                );
                Alignment {
                    query: q,
                    target: t,
                    score,
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

fn tail_prob(pmf: &[f64], bin: usize) -> f64 {
    if bin >= pmf.len() {
        return 0.0;
    }
    pmf[bin..].iter().fold(0.0, |a, &b| a + b).clamp(0.0, 1.0)
}

/// Bonferroni-corrected P-values of alignments.
pub fn alignment_pvalues(
    motifs: &[Motif],
    alignments: &[Alignment],
    opts: &CompareOptions,
) -> Result<Vec<f64>> {
    let cols = columns(motifs)?;
    let flat: Vec<Column> = cols.iter().flatten().copied().collect();
    let grid = (opts.null == NullModel::Parametric).then(|| enum_columns(&opts.bkg));
    let col_pmf = |q: &Column| -> Vec<f64> {
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
    let cache = |q: &[Column]| -> Vec<Vec<Vec<f64>>> {
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
    let per_query: Vec<Vec<(usize, f64)>> = with_threads(opts.nthreads, || {
        by_query
            .par_iter()
            .enumerate()
            .map(|(q, ks)| {
                if ks.is_empty() {
                    return Vec::new();
                }
                let fwd = cache(&cols[q]);
                let need_rc = ks.iter().any(|&k| alignments[k].strand == Strand::Minus);
                let rev = need_rc.then(|| cache(&rc_columns(&cols[q])));
                ks.iter()
                    .map(|&k| {
                        let a = &alignments[k];
                        let mut p = 1.0;
                        if a.overlap > 0 && a.query_start < cols[q].len() && a.n_tested > 0 {
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

/// One row of a comparison: a target motif's best alignment to a query.
#[derive(Debug, Clone, PartialEq)]
pub struct Comparison {
    pub query: usize,
    pub target: usize,
    pub offset: i64,
    pub strand: Strand,
    pub overlap: usize,
    /// Mean column correlation over the overlap.
    pub score: f64,
    pub pvalue: f64,
    /// Benjamini-Hochberg adjusted P-value, among the comparisons of the
    /// same query.
    pub qvalue: f64,
    /// Consensus of the aligned query columns (in the query's own
    /// orientation).
    pub query_consensus: String,
    /// Consensus of the aligned target columns.
    pub target_consensus: String,
}

fn consensus_range(cols: &[Column], start: i64, len: usize) -> String {
    let start = start.max(0) as usize;
    let len = len.min(cols.len().saturating_sub(start));
    cols[start..start + len]
        .iter()
        .map(consensus_char)
        .collect()
}

/// Compare each query motif with every motif (itself included).
///
/// Returns one [`Comparison`] per (query, target) pair, ordered by query and
/// then target, without any filtering.
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
    let cols = columns(motifs)?;

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
                query_consensus: consensus_range(&cols[a.query], q_start, a.overlap),
                target_consensus: consensus_range(&cols[a.target], t_start, a.overlap),
            }
        })
        .collect())
}

/// All-against-all comparison results, indexed `[query][target]`.
#[derive(Debug, Clone, PartialEq)]
pub struct ComparisonMatrices {
    /// Mean column correlations.
    pub scores: Vec<Vec<f64>>,
    pub pvalues: Vec<Vec<f64>>,
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
