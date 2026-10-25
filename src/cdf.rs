//! Score distributions of motifs.
//!
//! Motif scores are computed on an integer grid: each PWM entry is multiplied
//! by 1000 and truncated toward zero, and a site's score is the sum of the
//! entries along it. [`ScoreDistribution`] holds the exact distribution of
//! that score for random sequences drawn from a background, computed by
//! dynamic programming over the motif's positions.

use crate::matrix::{Matrix, sum};

/// A PWM value on the integer score grid: `x * 1000`, truncated toward zero.
pub(crate) fn grid_value(x: f64) -> i32 {
    (x * 1000.0) as i32
}

/// Round a score up to the integer grid, so that a threshold means "at
/// least this score". Values within floating-point round-trip error of a
/// grid point are snapped to it first, so that for example `1.001` maps to
/// 1001 rather than 1002.
pub(crate) fn score_grid_ceiling(score: f64) -> f64 {
    let scaled = score * 1000.0;
    if !scaled.is_finite() {
        return scaled;
    }
    let nearest = scaled.round();
    let tolerance = 4.0 * f64::EPSILON * f64::max(1.0, scaled.abs());
    if (scaled - nearest).abs() <= tolerance {
        nearest.ceil()
    } else {
        scaled.ceil()
    }
}

/// The distribution of scores a motif gives random sequences.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoreDistribution {
    /// `cdf[i]` is the probability of a score of at least `i + offset` grid
    /// units.
    cdf: Vec<f64>,
    offset: f64,
    score_min: i64,
    score_max: i64,
}

/// Scores filled at a time by [`ScoreDistribution::new`], 32 KB worth.
const DP_BLOCK: usize = 4096;

impl ScoreDistribution {
    /// Build the distribution for a PWM and background frequencies, one per
    /// row. The PWM must be finite.
    pub fn new(pwm: &Matrix, bkg: &[f64]) -> ScoreDistribution {
        let width = pwm.ncol();
        let ints: Vec<Vec<i32>> = pwm
            .columns()
            .map(|c| c.iter().map(|&v| grid_value(v)).collect())
            .collect();

        let motif_min = ints.iter().flatten().copied().fold(0, i32::min);
        let raw_min = pwm.as_slice().iter().copied().fold(0.0, f64::min);
        let offset = (raw_min * 1000.0).trunc() * width as f64;

        let (mut score_min, mut score_max) = (0i64, 0i64);
        for col in &ints {
            score_min += *col.iter().min().unwrap() as i64;
            score_max += *col.iter().max().unwrap() as i64;
        }

        let shifted: Vec<Vec<usize>> = ints
            .iter()
            .map(|c| c.iter().map(|&v| (v - motif_min) as usize).collect())
            .collect();
        let max_score = shifted.iter().flatten().copied().max().unwrap_or(0);

        // Dynamic programming over positions.
        // Only the reachable scores [lo, hi] are visited; the rest are zero.
        let len = width * max_score + 1;
        let mut pdf_new = vec![0.0; len];
        let mut pdf_old = vec![0.0; len];
        pdf_new[0] = 1.0;
        let (mut lo, mut hi) = (0usize, 0usize);
        // The scores holding values in the buffer that becomes pdf_new next,
        // which are cleared before it is filled again.
        let mut stale = (0usize, 0usize);
        for (i, col) in shifted.iter().enumerate() {
            let max_step = i * max_score;
            std::mem::swap(&mut pdf_old, &mut pdf_new);
            let (lo_k, hi_k) = (lo, hi.min(max_step));
            let (col_lo, col_hi) = (*col.iter().min().unwrap(), *col.iter().max().unwrap());
            let (new_lo, new_hi) = (lo_k + col_lo, hi_k + col_hi);
            pdf_new[stale.0.min(new_lo)..=stale.1.max(new_hi)].fill(0.0);
            // Scores are filled a block at a time, each block taking every
            // letter in turn, so that a block and the scores it reads stay in
            // cache. Every score still sums its terms in letter order. MEME
            // skips zero entries, but adding their zero products changes
            // nothing and lets the loop vectorize.
            let mut start = new_lo;
            while start <= new_hi {
                let end = (start + DP_BLOCK).min(new_hi + 1);
                for (&s, &b) in col.iter().zip(bkg) {
                    let from = start.max(lo_k + s);
                    let to = end.min(hi_k + s + 1);
                    if from >= to {
                        continue;
                    }
                    let src = &pdf_old[from - s..to - s];
                    let dst = &mut pdf_new[from..to];
                    for (d, &o) in dst.iter_mut().zip(src) {
                        *d += o * b;
                    }
                }
                start = end;
            }
            stale = (lo_k, hi_k);
            lo += col_lo;
            hi += col_hi;
        }
        let total = sum(&pdf_new);
        if total > 0.0 {
            pdf_new.iter_mut().for_each(|v| *v /= total);
        }
        for i in (0..len.saturating_sub(1)).rev() {
            pdf_new[i] += pdf_new[i + 1];
        }

        ScoreDistribution {
            cdf: pdf_new,
            offset,
            score_min,
            score_max,
        }
    }

    /// The probability of a random site scoring at least `score`.
    pub fn pvalue(&self, score: f64) -> f64 {
        let idx = score_grid_ceiling(score) - self.offset;
        if idx < 0.0 {
            1.0
        } else if idx < self.cdf.len() as f64 {
            self.cdf[idx as usize].min(1.0)
        } else {
            0.0
        }
    }

    /// The lowest score whose P-value is at most `pvalue`, or infinity if no
    /// score is that unlikely.
    pub fn score(&self, pvalue: f64) -> f64 {
        if pvalue >= 1.0 {
            return self.score_min as f64 / 1000.0;
        }
        for score in self.score_min..=self.score_max {
            let index = (score as f64 - self.offset) as usize;
            if self.cdf[index] <= pvalue {
                return score as f64 / 1000.0;
            }
        }
        f64::INFINITY
    }

    /// The lowest and highest possible scores.
    pub fn score_range(&self) -> (f64, f64) {
        (
            self.score_min as f64 / 1000.0,
            self.score_max as f64 / 1000.0,
        )
    }

    /// The part of the distribution at or above `from`, which answers
    /// [`ScoreDistribution::pvalue`] for those scores in a fraction of the
    /// memory.
    pub fn tail(&self, from: f64) -> ScoreTail {
        // Nothing scores above score_max, so the table stops there.
        let stop = ((self.score_max as f64 - self.offset) as usize + 1).min(self.cdf.len());
        let start = ((score_grid_ceiling(from) - self.offset).max(0.0) as usize).min(stop);
        ScoreTail {
            cdf: self.cdf[start..stop].to_vec(),
            offset: self.offset + start as f64,
        }
    }
}

/// The upper tail of a [`ScoreDistribution`], as returned by
/// [`ScoreDistribution::tail`].
#[derive(Debug, Clone, PartialEq)]
pub struct ScoreTail {
    cdf: Vec<f64>,
    offset: f64,
}

impl ScoreTail {
    /// The number of grid points held.
    pub fn len(&self) -> usize {
        self.cdf.len()
    }

    /// Whether no grid points are held.
    pub fn is_empty(&self) -> bool {
        self.cdf.is_empty()
    }

    /// The probability of a random site scoring at least `score`, which
    /// must be at or above the tail's start.
    pub fn pvalue(&self, score: f64) -> f64 {
        let idx = score_grid_ceiling(score) - self.offset;
        debug_assert!(idx >= 0.0, "score below the tail");
        if idx < self.cdf.len() as f64 {
            self.cdf[idx.max(0.0) as usize].min(1.0)
        } else {
            0.0
        }
    }
}
