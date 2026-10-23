//! Scores of individual matches: score ranges, the score of a word, every
//! word above a score and the probability of a word under the background.
//!
//! Scores use the motif's PWM on the same integer grid as scanning: each entry
//! is multiplied by 1000 and truncated. Unless non-finite values are allowed, a
//! PWM with `-Inf` entries first gets a pseudocount from
//! [`crate::ops::normalize`]; when they are allowed, a word using such an entry
//! scores `-Inf`.

use crate::alphabet::Alphabet;
use crate::cdf::grid_value;
use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::motif::{Motif, MotifType};
use crate::ops::normalize;

/// The motif as a PWM for scoring, with a pseudocount added to remove
/// `-Inf` entries unless `allow_nonfinite`.
fn scoring_pwm(motif: &Motif, allow_nonfinite: bool) -> Result<Motif> {
    let pwm = convert_type(motif, MotifType::Pwm, None, false)?;
    if !allow_nonfinite && pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
        normalize(&pwm)
    } else {
        Ok(pwm)
    }
}

/// Grid values per column, `None` for non-finite entries.
fn grid_columns(pwm: &Motif) -> Vec<Vec<Option<i64>>> {
    pwm.matrix
        .columns()
        .map(|c| {
            c.iter()
                .map(|&v| v.is_finite().then(|| grid_value(v) as i64))
                .collect()
        })
        .collect()
}

/// The lowest and highest possible scores. With non-finite values allowed,
/// the lowest score can be `-Inf`.
pub fn motif_range(motif: &Motif, allow_nonfinite: bool) -> Result<(f64, f64)> {
    let cols = grid_columns(&scoring_pwm(motif, allow_nonfinite)?);
    let mut lo = 0.0;
    let mut hi = 0;
    for c in &cols {
        lo += match c.iter().copied().collect::<Option<Vec<i64>>>() {
            Some(v) => *v.iter().min().unwrap() as f64,
            None => f64::NEG_INFINITY,
        };
        hi += c.iter().flatten().max().copied().unwrap_or(0);
    }
    Ok((lo / 1000.0, hi as f64 / 1000.0))
}

/// Where a threshold of zero lies for [`motif_score`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ThresholdType {
    /// Zero is the lowest possible score and one the highest.
    #[default]
    Total,
    /// Zero is a score of zero and one the highest.
    FromZero,
}

impl ThresholdType {
    /// The name [`str::parse`] reads, such as `"total"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            ThresholdType::Total => "total",
            ThresholdType::FromZero => "fromzero",
        }
    }
}

impl std::str::FromStr for ThresholdType {
    type Err = crate::error::Error;

    /// Parse `"total"` or `"fromzero"`.
    fn from_str(s: &str) -> Result<ThresholdType> {
        match s {
            "total" => Ok(ThresholdType::Total),
            "fromzero" => Ok(ThresholdType::FromZero),
            _ => arg_err("threshold type must be 'total' or 'fromzero'"),
        }
    }
}

impl std::fmt::Display for ThresholdType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Scores at fractions of the motif's score range, from 0 to 1. Non-finite
/// entries are left out of the lowest score, so the `Total` type is only
/// meaningful without them.
pub fn motif_score(
    motif: &Motif,
    thresholds: &[f64],
    threshold_type: ThresholdType,
    allow_nonfinite: bool,
) -> Result<Vec<f64>> {
    if thresholds.iter().any(|t| !(0.0..=1.0).contains(t)) {
        return arg_err("thresholds must be between 0 and 1");
    }
    let cols = grid_columns(&scoring_pwm(motif, allow_nonfinite)?);
    let s_max: i64 = cols
        .iter()
        .map(|c| c.iter().flatten().max().copied().unwrap_or(0))
        .sum();
    let s_min: i64 = cols
        .iter()
        .map(|c| c.iter().flatten().min().copied().unwrap_or(0))
        .sum();
    let (s_max, s_min) = (s_max as f64, s_min as f64);
    Ok(thresholds
        .iter()
        .map(|&t| {
            let out = match threshold_type {
                ThresholdType::Total => (s_max.abs() + s_min.abs()) * t - s_min.abs(),
                ThresholdType::FromZero => s_max * t,
            };
            out / 1000.0
        })
        .collect())
}

/// The row of each letter of `word`, checking it against the motif's width
/// and alphabet. DNA, RNA and amino acid letters match in either case.
fn word_rows(motif: &Motif, word: &str) -> Result<Vec<usize>> {
    let letters = motif.alphabet.letters();
    let exact = matches!(motif.alphabet, Alphabet::Custom(_));
    let rows: Option<Vec<usize>> = word
        .chars()
        .map(|c| {
            letters.iter().position(|&l| {
                if exact {
                    l == c
                } else {
                    l.eq_ignore_ascii_case(&c)
                }
            })
        })
        .collect();
    let n = word.chars().count();
    if n != motif.width() {
        return arg_err(format!(
            "the match has {n} letters but the motif has {} positions",
            motif.width()
        ));
    }
    rows.ok_or_else(|| {
        crate::error::Error::Argument(format!(
            "found letters in '{word}' that are not in the motif's alphabet"
        ))
    })
}

/// The score of each word: each PWM entry is truncated to three decimals,
/// and so is the total.
pub fn score_match<S: AsRef<str>>(
    motif: &Motif,
    words: &[S],
    allow_nonfinite: bool,
) -> Result<Vec<f64>> {
    let pwm = scoring_pwm(motif, allow_nonfinite)?;
    let cols: Vec<Vec<f64>> = pwm.matrix.columns().map(|c| c.to_vec()).collect();
    // Scores are the sum of the entries on the integer grid, as in scanning,
    // so that a word scores exactly what scan_sequences() reports for it.
    words
        .iter()
        .map(|w| {
            let rows = word_rows(&pwm, w.as_ref())?;
            let mut grid: i64 = 0;
            let mut other = 0.0;
            for (&r, c) in rows.iter().zip(&cols) {
                if c[r].is_finite() {
                    grid += grid_value(c[r]) as i64;
                } else {
                    other += c[r];
                }
            }
            Ok(grid as f64 / 1000.0 + other)
        })
        .collect()
}

/// Every possible score of the motif, highest first. There are as many as the alphabet size to the power of
/// the width, so this is only for small motifs.
pub fn get_scores(motif: &Motif, allow_nonfinite: bool) -> Result<Vec<f64>> {
    let cols = grid_columns(&scoring_pwm(motif, allow_nonfinite)?);
    let n = (motif.alphabet.len() as f64).powi(cols.len() as i32);
    if n > i32::MAX as f64 {
        return arg_err(format!(
            "the motif has {n:.0} possible scores, too many to list; use a narrower motif"
        ));
    }
    let mut scores: Vec<Option<i64>> = vec![Some(0)];
    for c in &cols {
        scores = c
            .iter()
            .flat_map(|&v| scores.iter().map(move |&s| s.zip(v).map(|(a, b)| a + b)))
            .collect();
    }
    let mut out: Vec<f64> = scores
        .into_iter()
        .map(|s| s.map_or(f64::NEG_INFINITY, |s| s as f64 / 1000.0))
        .collect();
    out.sort_unstable_by(|a, b| b.total_cmp(a));
    Ok(out)
}

/// Every word scoring at least `score`, found by branch and bound. Within
/// each column letters are ranked by score, and words come out ordered by
/// those ranks, with the column holding the highest-scoring letter varying
/// fastest. A score of `+Inf` has no matches; a score outside the motif's
/// range is an error.
///
/// # Example
///
/// ```
/// use universalmotif::create::{CreateOptions, from_strings};
/// use universalmotif::scoring::{get_matches, motif_range};
///
/// let opts = CreateOptions { nsites: Some(20.0), pseudocount: 1.0, ..Default::default() };
/// let m = from_strings(&["TATAAA"], None, &opts)?;
/// let (_, best) = motif_range(&m, false)?;
/// assert_eq!(get_matches(&m, best, false)?, ["TATAAA"]);
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn get_matches(motif: &Motif, score: f64, allow_nonfinite: bool) -> Result<Vec<String>> {
    let pwm = scoring_pwm(motif, allow_nonfinite)?;
    if score == f64::INFINITY {
        return Ok(Vec::new());
    }
    let (lo, hi) = if pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
        let hi = motif_score(&pwm, &[1.0], ThresholdType::FromZero, true)?[0];
        (f64::NEG_INFINITY, hi)
    } else {
        let s = motif_score(&pwm, &[0.0, 1.0], ThresholdType::Total, true)?;
        (s[0], s[1])
    };
    if score > hi {
        return arg_err(format!(
            "the score is above the highest possible score, {hi}"
        ));
    }
    if score < lo {
        return arg_err(format!(
            "the score is below the lowest possible score, {lo}"
        ));
    }
    let target = (score * 1000.0).trunc() as i64;
    let width = pwm.width();
    // Non-finite entries take a stand-in for minus infinity, low enough to
    // rule out any word using one but small enough that a whole word of
    // them still fits in an i32.
    let floor = -(i32::MAX as i64) / width as i64;
    let cols: Vec<Vec<i64>> = grid_columns(&pwm)
        .into_iter()
        .map(|c| c.into_iter().map(|v| v.unwrap_or(floor)).collect())
        .collect();

    // Letters sorted by score within each column and columns by their best
    // score, both highest first and keeping ties in order.
    let order: Vec<Vec<usize>> = cols
        .iter()
        .map(|c| {
            let mut o: Vec<usize> = (0..c.len()).collect();
            o.sort_by_key(|&i| std::cmp::Reverse(c[i]));
            o
        })
        .collect();
    let mut col_order: Vec<usize> = (0..width).collect();
    col_order.sort_by_key(|&j| std::cmp::Reverse(cols[j][order[j][0]]));
    let sorted: Vec<Vec<i64>> = col_order
        .iter()
        .map(|&j| order[j].iter().map(|&i| cols[j][i]).collect())
        .collect();
    // The best score still possible after each column.
    let mut rest = vec![0i64; width + 1];
    for j in (0..width).rev() {
        rest[j] = rest[j + 1] + sorted[j][0];
    }

    let n = sorted[0].len();
    let mut paths: Vec<Vec<usize>> = (0..n)
        .filter(|&i| sorted[0][i] + rest[1] >= target)
        .map(|i| vec![i])
        .collect();
    let mut scores: Vec<i64> = paths.iter().map(|p| sorted[0][p[0]]).collect();
    for (j, col) in sorted.iter().enumerate().skip(1) {
        let mut next = Vec::new();
        let mut next_scores = Vec::new();
        for (i, &v) in col.iter().enumerate() {
            for (p, &s) in paths.iter().zip(&scores) {
                if s + v + rest[j + 1] >= target {
                    let mut q = p.clone();
                    q.push(i);
                    next.push(q);
                    next_scores.push(s + v);
                }
            }
        }
        paths = next;
        scores = next_scores;
    }

    let letters = pwm.alphabet.letters();
    let mut out = Vec::with_capacity(paths.len());
    for p in paths {
        let mut word = vec![' '; width];
        for (k, &i) in p.iter().enumerate() {
            let j = col_order[k];
            word[j] = letters[order[j][i]];
        }
        out.push(word.into_iter().collect());
    }
    Ok(out)
}

/// The probability of each word under the motif's background. Unless
/// `allow_zero`, a background with zeros gets the motif's pseudocount spread
/// over the letters and scaled by its `nsites`, taking 1 and 100 respectively
/// when the motif has none.
pub fn prob_match<S: AsRef<str>>(motif: &Motif, words: &[S], allow_zero: bool) -> Result<Vec<f64>> {
    let mut bkg = motif.bkg.clone();
    if !allow_zero && bkg.contains(&0.0) {
        let pseudo = if motif.pseudocount == 0.0 {
            1.0
        } else {
            motif.pseudocount
        };
        let nsites = motif.nsites.unwrap_or(100.0);
        let n = bkg.len() as f64;
        bkg.iter_mut().for_each(|b| *b += (pseudo / nsites) / n);
    }
    words
        .iter()
        .map(|w| {
            let rows = word_rows(motif, w.as_ref())?;
            Ok(rows.iter().fold(1.0, |p, &r| p * bkg[r]))
        })
        .collect()
}

/// The probability of each word under a background given as letter
/// frequencies.
pub fn prob_match_bkg<S: AsRef<str>>(bkg: &[(char, f64)], words: &[S]) -> Result<Vec<f64>> {
    words
        .iter()
        .map(|w| {
            w.as_ref()
                .chars()
                .try_fold(1.0, |p, c| match bkg.iter().find(|(l, _)| *l == c) {
                    Some((_, f)) => Ok(p * f),
                    None => arg_err(format!("the letter '{c}' is not in the background")),
                })
        })
        .collect()
}
