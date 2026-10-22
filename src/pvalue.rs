//! Motif P-values and score thresholds (`motif_pvalue()`).
//!
//! A motif's P-value for a score is the probability that a random site, drawn
//! from the background, scores at least that much. It is computed exactly
//! from the motif's score distribution (see [`crate::cdf`]), which is
//! universalmotif's `method = "dynamic"`. The exhaustive branch-and-bound
//! method and its random score search are not ported.

use rayon::prelude::*;

use crate::cdf::ScoreDistribution;
use crate::convert::convert_type;
use crate::error::{Error, Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};
use crate::ops::normalize;

/// Options for [`motif_pvalue`] and [`motif_score`].
#[derive(Debug, Clone, PartialEq)]
pub struct PvalueOptions {
    /// Background frequencies: none to use each motif's own, one set for all
    /// motifs, or one set per motif.
    pub bkg: Vec<Vec<f64>>,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for PvalueOptions {
    fn default() -> Self {
        PvalueOptions {
            bkg: Vec::new(),
            nthreads: 1,
        }
    }
}

/// Run `f` on a thread pool of `nthreads` threads (0 for all available).
/// Parallel iterators inside `f` use this pool, so they never spread to
/// more threads than asked for.
pub(crate) fn with_threads<T: Send>(nthreads: usize, f: impl FnOnce() -> T + Send) -> Result<T> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .map_err(|e| Error::Argument(format!("could not start threads: {e}")))?;
    Ok(pool.install(f))
}

/// The PWM and background a motif is scored with.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoringMatrix {
    pub pwm: Matrix,
    pub bkg: Vec<f64>,
    /// Whether the pseudocount had to be reapplied because the PWM had
    /// infinite values.
    pub pseudocount_added: bool,
}

/// Prepare the PWM a motif is scored with, as `motif_pvalue()` does.
///
/// A finite PWM is used as it is when the background is its own, because
/// converting back and forth can move entries across integer score
/// boundaries. Otherwise the motif is converted to a PPM and then to a PWM
/// against `bkg` (its own background by default). If that PWM would have
/// infinite values the pseudocount is reapplied first (see [`normalize`]).
pub fn scoring_matrix(motif: &Motif, bkg: Option<&[f64]>) -> Result<ScoringMatrix> {
    let bkg = match bkg {
        Some(b) => {
            if b.len() != motif.alphabet.len() {
                return arg_err(format!(
                    "the background must have {} values, one per letter",
                    motif.alphabet.len()
                ));
            }
            b.to_vec()
        }
        None => motif.bkg.clone(),
    };
    if motif.motif_type == MotifType::Pwm
        && motif.matrix.as_slice().iter().all(|v| v.is_finite())
        && bkg == motif.bkg
    {
        return Ok(ScoringMatrix {
            pwm: motif.matrix.clone(),
            bkg,
            pseudocount_added: false,
        });
    }
    let mut ppm = convert_type(motif, MotifType::Ppm, None, false)?;
    let pwm = convert_type(&ppm, MotifType::Pwm, None, false)?;
    let mut pseudocount_added = false;
    if pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
        ppm = normalize(&ppm)?;
        pseudocount_added = true;
    }
    let mut parts = ppm.into_parts();
    parts.bkg = bkg.clone();
    let ppm = Motif::from_parts(parts)?;
    let pwm = convert_type(&ppm, MotifType::Pwm, None, false)?;
    Ok(ScoringMatrix {
        pwm: pwm.into_parts().matrix,
        bkg,
        pseudocount_added,
    })
}

/// The score distribution of each motif, built in parallel.
pub fn score_distributions(
    motifs: &[Motif],
    opts: &PvalueOptions,
) -> Result<Vec<ScoreDistribution>> {
    if opts.bkg.len() > 1 && opts.bkg.len() != motifs.len() {
        return arg_err("give one background for all motifs or one per motif");
    }
    let mats = motifs
        .iter()
        .enumerate()
        .map(|(i, m)| {
            let bkg = match opts.bkg.len() {
                0 => None,
                1 => Some(opts.bkg[0].as_slice()),
                _ => Some(opts.bkg[i].as_slice()),
            };
            scoring_matrix(m, bkg)
        })
        .collect::<Result<Vec<_>>>()?;
    with_threads(opts.nthreads, || {
        mats.par_iter()
            .map(|m| ScoreDistribution::new(&m.pwm, &m.bkg))
            .collect()
    })
}

fn check_queries(n: usize, queries: &[(usize, f64)]) -> Result<()> {
    for &(i, v) in queries {
        if i >= n {
            return arg_err(format!("query refers to motif {i}, but there are {n}"));
        }
        if !v.is_finite() {
            return arg_err("scores and P-values must be finite numbers");
        }
    }
    Ok(())
}

/// Answer `(motif index, value)` queries one motif at a time, so only the
/// distributions in flight are held in memory.
fn per_motif(
    motifs: &[Motif],
    queries: &[(usize, f64)],
    opts: &PvalueOptions,
    answer: impl Fn(&ScoreDistribution, f64) -> f64 + Sync,
) -> Result<Vec<f64>> {
    if opts.bkg.len() > 1 && opts.bkg.len() != motifs.len() {
        return arg_err("give one background for all motifs or one per motif");
    }
    let mut by_motif: Vec<Vec<usize>> = vec![Vec::new(); motifs.len()];
    for (q, &(i, _)) in queries.iter().enumerate() {
        by_motif[i].push(q);
    }
    let answers: Vec<Vec<(usize, f64)>> = with_threads(opts.nthreads, || {
        by_motif
            .par_iter()
            .enumerate()
            .filter(|(_, qs)| !qs.is_empty())
            .map(|(i, qs)| {
                let bkg = match opts.bkg.len() {
                    0 => None,
                    1 => Some(opts.bkg[0].as_slice()),
                    _ => Some(opts.bkg[i].as_slice()),
                };
                let m = scoring_matrix(&motifs[i], bkg)?;
                let dist = ScoreDistribution::new(&m.pwm, &m.bkg);
                Ok(qs
                    .iter()
                    .map(|&q| (q, answer(&dist, queries[q].1)))
                    .collect())
            })
            .collect::<Result<Vec<_>>>()
    })??;
    let mut out = vec![0.0; queries.len()];
    for (q, v) in answers.into_iter().flatten() {
        out[q] = v;
    }
    Ok(out)
}

/// P-values for `(motif index, score)` pairs.
pub fn motif_pvalue(
    motifs: &[Motif],
    queries: &[(usize, f64)],
    opts: &PvalueOptions,
) -> Result<Vec<f64>> {
    check_queries(motifs.len(), queries)?;
    per_motif(motifs, queries, opts, |d, s| d.pvalue(s))
}

/// Score thresholds for `(motif index, P-value)` pairs: the lowest score a
/// random site reaches with at most that probability, or infinity if no
/// score is that unlikely.
pub fn motif_score(
    motifs: &[Motif],
    queries: &[(usize, f64)],
    opts: &PvalueOptions,
) -> Result<Vec<f64>> {
    check_queries(motifs.len(), queries)?;
    per_motif(motifs, queries, opts, |d, p| d.score(p))
}
