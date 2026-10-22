//! Column similarity metrics and averages (universalmotif's
//! `compare_columns()` and `average_ic()`).
//!
//! [`crate::compare`] uses Pearson correlation only, as `compare_motifs()`
//! lite does; the other metrics are here for comparing single columns. The
//! arithmetic follows `src/compare_motifs.cpp`.

use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::motif::{Motif, MotifType};

/// A column similarity or distance metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnMetric {
    /// Pearson correlation coefficient.
    Pcc,
    /// Euclidean distance.
    Eucl,
    /// Sandelin-Wasserman similarity.
    Sw,
    /// Kullback-Leibler divergence (symmetrised).
    Kl,
    /// Background-weighted Euclidean distance.
    Weucl,
    /// Average log-likelihood ratio.
    Allr,
    /// Bhattacharyya coefficient.
    Bhat,
    /// Hellinger distance.
    Hell,
    /// Background-weighted Pearson correlation.
    Wpcc,
    /// Squared Euclidean distance.
    Seucl,
    /// Manhattan distance.
    Man,
    /// Average log-likelihood ratio with a lower limit of -2.
    AllrLl,
}

impl ColumnMetric {
    pub fn parse(s: &str) -> Result<ColumnMetric> {
        Ok(match s {
            "PCC" => ColumnMetric::Pcc,
            "EUCL" => ColumnMetric::Eucl,
            "SW" => ColumnMetric::Sw,
            "KL" => ColumnMetric::Kl,
            "WEUCL" => ColumnMetric::Weucl,
            "ALLR" => ColumnMetric::Allr,
            "BHAT" => ColumnMetric::Bhat,
            "HELL" => ColumnMetric::Hell,
            "WPCC" => ColumnMetric::Wpcc,
            "SEUCL" => ColumnMetric::Seucl,
            "MAN" => ColumnMetric::Man,
            "ALLR_LL" => ColumnMetric::AllrLl,
            _ => {
                return arg_err(format!(
                    "unknown metric '{s}'; use one of PCC, EUCL, SW, KL, WEUCL, ALLR, BHAT, \
                     HELL, WPCC, SEUCL, MAN, ALLR_LL"
                ));
            }
        })
    }

    fn needs_bkg(self) -> bool {
        matches!(
            self,
            ColumnMetric::Weucl | ColumnMetric::Wpcc | ColumnMetric::Allr | ColumnMetric::AllrLl
        )
    }
}

fn total(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(0.0, |a, b| a + b)
}

/// Compare two columns (probabilities over the same letters) with `metric`.
/// The backgrounds are used by the weighted metrics and ALLR, which also
/// uses the site counts. As in the R package, a column whose first value is
/// negative scores 0.
pub fn compare_columns(
    x: &[f64],
    y: &[f64],
    metric: ColumnMetric,
    bkg1: &[f64],
    bkg2: &[f64],
    nsites1: f64,
    nsites2: f64,
) -> Result<f64> {
    if x.len() != y.len() {
        return arg_err(format!(
            "length(x) does not match length(y) [{}, {}]",
            x.len(),
            y.len()
        ));
    }
    if x.len() < 2 {
        return arg_err("columns should have at least 2 entries");
    }
    if metric.needs_bkg() && (bkg1.len() != x.len() || bkg2.len() != x.len()) {
        return arg_err("incorrect background vector length");
    }
    if matches!(metric, ColumnMetric::Allr | ColumnMetric::AllrLl)
        && (nsites1 <= 1.0 || nsites2 <= 1.0)
    {
        return arg_err("nsites1/nsites2 should be greater than 1");
    }
    if !(x[0] >= 0.0 && y[0] >= 0.0) {
        return Ok(0.0);
    }
    let pairs = || x.iter().zip(y);
    let n = x.len() as f64;
    Ok(match metric {
        ColumnMetric::Eucl => total(pairs().map(|(a, b)| (a - b) * (a - b))).sqrt(),
        ColumnMetric::Seucl => total(pairs().map(|(a, b)| (a - b) * (a - b))),
        ColumnMetric::Sw => 2.0 - total(pairs().map(|(a, b)| (a - b) * (a - b))),
        ColumnMetric::Man => total(pairs().map(|(a, b)| (a - b).abs())),
        ColumnMetric::Bhat => total(pairs().map(|(a, b)| (a * b).sqrt())),
        ColumnMetric::Hell => {
            total(pairs().map(|(a, b)| (a.sqrt() - b.sqrt()) * (a.sqrt() - b.sqrt()))).sqrt()
                / 2f64.sqrt()
        }
        // The two terms are added one at a time, as the C++ does.
        ColumnMetric::Kl => {
            0.5 * pairs().fold(0.0, |s, (a, b)| s + a * (a / b).ln() + b * (b / a).ln())
        }
        ColumnMetric::Weucl => total(
            pairs()
                .zip(bkg1.iter().zip(bkg2))
                .map(|((a, b), (w1, w2))| (a - b) * (a - b) * (w1 + w2)),
        )
        .sqrt(),
        ColumnMetric::Pcc => {
            let (sx, sy) = (total(x.iter().copied()), total(y.iter().copied()));
            let top = n * total(pairs().map(|(a, b)| a * b)) - sx * sy;
            let bot = ((n * total(x.iter().map(|a| a * a)) - sx * sx)
                * (n * total(y.iter().map(|b| b * b)) - sy * sy))
                .sqrt();
            if bot == 0.0 { 1.0 } else { top / bot }
        }
        ColumnMetric::Wpcc => {
            let w: Vec<f64> = bkg1.iter().zip(bkg2).map(|(a, b)| a + b).collect();
            let m1 = total(w.iter().zip(x).map(|(w, a)| w * a));
            let m2 = total(w.iter().zip(y).map(|(w, b)| w * b));
            let v1 = total(w.iter().zip(x).map(|(w, a)| w * (a - m1) * (a - m1)));
            let v2 = total(w.iter().zip(y).map(|(w, b)| w * (b - m2) * (b - m2)));
            let cov = total(
                w.iter()
                    .zip(pairs())
                    .map(|(w, (a, b))| w * (a - m1) * (b - m2)),
            );
            let bot = (v1 * v2).sqrt();
            if bot == 0.0 { 1.0 } else { cov / bot }
        }
        ColumnMetric::Allr | ColumnMetric::AllrLl => {
            let left = total(
                pairs()
                    .zip(bkg1)
                    .map(|((a, b), g)| (b * nsites2) * (a / g).ln()),
            );
            let right = total(
                pairs()
                    .zip(bkg2)
                    .map(|((a, b), g)| (a * nsites1) * (b / g).ln()),
            );
            let v = (0.0 + left + right) / (nsites1 + nsites2);
            if metric == ColumnMetric::AllrLl && v < -2.0 {
                -2.0
            } else {
                v
            }
        }
    })
}

/// How [`average`] summarises a set of values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Average {
    /// Arithmetic mean.
    #[default]
    ArithmeticMean,
    /// Geometric mean of the positive values (over all of them).
    GeometricMean,
    Median,
    /// Fisher z-transform: the inverse hyperbolic tangent of the mean of
    /// the hyperbolic tangents.
    FisherZ,
}

impl Average {
    pub fn parse(s: &str) -> Result<Average> {
        match s {
            "a.mean" => Ok(Average::ArithmeticMean),
            "g.mean" => Ok(Average::GeometricMean),
            "median" => Ok(Average::Median),
            "fzt" => Ok(Average::FisherZ),
            _ => arg_err("average must be one of a.mean, g.mean, median, fzt"),
        }
    }
}

/// Average `values` as universalmotif's `average_cpp()` does.
pub fn average(values: &[f64], how: Average) -> f64 {
    let n = values.len() as f64;
    match how {
        Average::ArithmeticMean => total(values.iter().copied()) / n,
        Average::GeometricMean => {
            let s = total(values.iter().filter(|&&v| v > 0.0).map(|v| v.ln()));
            if s == 0.0 { 0.0 } else { (s / n).exp() }
        }
        Average::Median => {
            if values.len() == 1 {
                return values[0];
            }
            let mut v = values.to_vec();
            v.sort_by(f64::total_cmp);
            let h = v.len() / 2;
            if v.len().is_multiple_of(2) {
                (v[h - 1] + v[h]) / 2.0
            } else {
                v[h]
            }
        }
        Average::FisherZ => (total(values.iter().map(|v| v.tanh())) / n).atanh(),
    }
}

/// The average information content of each motif's positions, as
/// universalmotif's `average_ic()`.
pub fn average_ic(motifs: &[Motif], how: Average) -> Result<Vec<f64>> {
    motifs
        .iter()
        .map(|m| {
            let icm = convert_type(m, MotifType::Icm, None, false)?;
            let ic: Vec<f64> = icm
                .matrix
                .columns()
                .map(|c| total(c.iter().copied()))
                .collect();
            Ok(average(&ic, how))
        })
        .collect()
}
