//! Gapped motifs as the ungapped motifs they stand for, one for each
//! combination of gap sizes, with the gaps filled by blank columns, so that
//! comparisons, layouts and merges can take the gaps into account.

use crate::compare::Comparison;
use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{GapInfo, Motif, MotifType};

/// Instances beyond which a gapped motif is refused.
const MAX_INSTANCES: usize = 1 << 12;

/// Recomputes the score of a comparison kept between gapped motifs, given
/// the comparison and its pair of instances.
pub(crate) type Rescore<'a> = &'a dyn Fn(&Comparison, usize, usize) -> f64;

/// How an instance aligns to the anchor's.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Fit {
    /// The alignment's score, higher being better.
    pub score: f64,
    /// Its overlap, in columns.
    pub overlap: usize,
    /// Whether the instance is reverse complemented.
    pub minus: bool,
}

/// A set of motifs with each gapped one replaced by its instances.
pub(crate) struct Expansion {
    /// The instances, those of each motif together and in its order.
    pub motifs: Vec<Motif>,
    /// The motif each instance stands for.
    pub origin: Vec<usize>,
    /// The size of each gap of each instance, empty for ungapped motifs.
    pub sizes: Vec<Vec<usize>>,
    /// The instances of each motif.
    pub of: Vec<Vec<usize>>,
}

/// Whether any motif has gaps.
pub(crate) fn any_gapped(motifs: &[Motif]) -> bool {
    motifs
        .iter()
        .any(|m| m.gap.isgapped && !m.gap.gaploc.is_empty())
}

/// Every combination of gap sizes of a motif, smallest first.
fn combos(m: &Motif) -> Result<Vec<Vec<usize>>> {
    let g = &m.gap;
    let n: usize = g
        .mingap
        .iter()
        .zip(&g.maxgap)
        .map(|(lo, hi)| hi.saturating_sub(*lo) + 1)
        .product();
    if n > MAX_INSTANCES {
        return arg_err(format!(
            "{} has {n} combinations of gap sizes, more than the {MAX_INSTANCES} that can be compared",
            m.name
        ));
    }
    let mut out = vec![Vec::new()];
    for (&lo, &hi) in g.mingap.iter().zip(&g.maxgap) {
        out = out
            .into_iter()
            .flat_map(|c| {
                (lo..=hi).map(move |s| {
                    let mut c = c.clone();
                    c.push(s);
                    c
                })
            })
            .collect();
    }
    Ok(out)
}

/// The motif with gaps of the given sizes filled by blank columns: zeros for
/// a CWM, and for other motifs, which become PPMs, the background
/// frequencies, `bkg` for DNA and RNA and uniform for other alphabets.
pub(crate) fn instance(m: &Motif, sizes: &[usize], bkg: &[f64; 4]) -> Result<Motif> {
    let base = if m.motif_type == MotifType::Cwm {
        m.clone()
    } else {
        convert_type(m, MotifType::Ppm, None, false)?
    };
    let n = base.alphabet.len();
    let blank: Vec<f64> = match base.motif_type {
        MotifType::Cwm => vec![0.0; n],
        _ if n == 4 && base.alphabet.is_nucleotide() => {
            let total: f64 = bkg.iter().sum();
            bkg.iter().map(|b| b / total).collect()
        }
        _ => vec![1.0 / n as f64; n],
    };
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(base.width() + sizes.iter().sum::<usize>());
    let mut next = 0;
    for (k, col) in base.matrix.columns().enumerate() {
        if next < sizes.len() && m.gap.gaploc[next] == k {
            for _ in 0..sizes[next] {
                cols.push(blank.clone());
            }
            next += 1;
        }
        cols.push(col.to_vec());
    }
    let mut parts = base.into_parts();
    parts.matrix = Matrix::from_columns(&cols)?;
    parts.gap = GapInfo::default();
    parts.recompute_derived();
    Motif::from_parts(parts)
}

/// Each motif as its instances, or as itself when it has no gaps.
pub(crate) fn expand(motifs: &[Motif], bkg: &[f64; 4]) -> Result<Expansion> {
    let mut out = Expansion {
        motifs: Vec::new(),
        origin: Vec::new(),
        sizes: Vec::new(),
        of: Vec::with_capacity(motifs.len()),
    };
    for (i, m) in motifs.iter().enumerate() {
        let mut mine = Vec::new();
        if m.gap.isgapped && !m.gap.gaploc.is_empty() {
            for sizes in combos(m)? {
                mine.push(out.motifs.len());
                out.motifs.push(instance(m, &sizes, bkg)?);
                out.origin.push(i);
                out.sizes.push(sizes);
            }
        } else {
            mine.push(out.motifs.len());
            out.motifs.push(m.clone());
            out.origin.push(i);
            out.sizes.push(Vec::new());
        }
        out.of.push(mine);
    }
    Ok(out)
}

/// Comparisons of instances reduced to those of the motifs they stand for.
/// `rows` compares every instance of each of `queries` with every instance,
/// query by query. Each pair of motifs keeps the comparison of its instances
/// with the lowest P-value (or, without P-values, the best score: the
/// highest, or with `lower_better` the lowest), its P-value multiplied by
/// the number of pairs of instances tried, and q-values are worked out
/// again within each query, leaving out its comparison with itself. With
/// `rescore`, the score of each comparison kept between gapped motifs is
/// replaced by what it gives for the comparison and its pair of instances.
pub(crate) fn collapse(
    rows: Vec<Comparison>,
    exp: &Expansion,
    queries: &[usize],
    lower_better: bool,
    rescore: Option<Rescore>,
) -> Vec<Comparison> {
    let n = exp.of.len();
    let has_p = rows.iter().any(|r| !r.pvalue.is_nan());
    let better = |a: &Comparison, b: &Comparison| -> bool {
        if has_p && a.pvalue != b.pvalue {
            return a.pvalue < b.pvalue;
        }
        if lower_better {
            a.score < b.score
        } else {
            a.score > b.score
        }
    };
    let mut best: Vec<Vec<Option<(Comparison, usize, usize)>>> = vec![vec![None; n]; n];
    for r in rows {
        let (q, t) = (exp.origin[r.query], exp.origin[r.target]);
        let (iq, it) = (r.query, r.target);
        let keep = match &best[q][t] {
            None => true,
            Some((b, _, _)) => {
                // Scores of comparisons that could not be made are NaN.
                !r.score.is_nan() && (b.score.is_nan() || better(&r, b))
            }
        };
        if keep {
            best[q][t] = Some((r, iq, it));
        }
    }
    let mut out = Vec::with_capacity(queries.len() * n);
    for &q in queries {
        let mut row: Vec<Comparison> = (0..n)
            .map(|t| {
                // Cloned rather than taken, since a query may be asked for more than once.
                let (mut r, iq, it) = best[q][t].clone().expect("every pair was compared");
                if let Some(f) = rescore
                    && !(exp.sizes[iq].is_empty() && exp.sizes[it].is_empty())
                {
                    r.score = f(&r, iq, it);
                }
                r.query = q;
                r.target = t;
                r.query_gaps = exp.sizes[iq].clone();
                r.target_gaps = exp.sizes[it].clone();
                r
            })
            .collect();
        if has_p {
            for (t, r) in row.iter_mut().enumerate() {
                let tried = (exp.of[q].len() * exp.of[t].len()) as f64;
                r.pvalue = (r.pvalue * tried).min(1.0);
            }
            let p: Vec<f64> = row.iter().map(|r| r.pvalue).collect();
            let q_adj = crate::stats::p_adjust_bh_except(&p, Some(q));
            for (r, v) in row.iter_mut().zip(q_adj) {
                r.qvalue = v;
            }
        }
        out.extend(row);
    }
    out
}

/// The ungapped parts of a gapped motif, split at its gaps.
pub(crate) fn split(m: &Motif) -> Result<Vec<Motif>> {
    let mut bounds = vec![0];
    bounds.extend(&m.gap.gaploc);
    bounds.push(m.width());
    let cols: Vec<Vec<f64>> = m.matrix.columns().map(<[f64]>::to_vec).collect();
    bounds
        .windows(2)
        .map(|b| {
            let mut parts = m.clone().into_parts();
            parts.matrix = Matrix::from_columns(&cols[b[0]..b[1]])?;
            parts.gap = GapInfo::default();
            parts.recompute_derived();
            Motif::from_parts(parts)
        })
        .collect()
}

/// The columns of an instance that fill its gaps.
pub(crate) fn blank_columns(m: &Motif, sizes: &[usize]) -> Vec<usize> {
    let mut out = Vec::new();
    let mut before = 0;
    for (k, &s) in sizes.iter().enumerate() {
        let start = m.gap.gaploc[k] + before;
        out.extend(start..start + s);
        before += s;
    }
    out
}
