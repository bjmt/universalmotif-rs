//! Column similarity metrics with [`compare_columns`], and averages with
//! [`average`] and [`average_ic`].
//!
//! [`crate::compare`] uses Pearson correlation only; the other metrics are
//! here for comparing single columns.

use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::motif::{Motif, MotifType};

/// A column similarity or distance metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ColumnMetric {
    /// Pearson correlation coefficient.
    Pcc,
    /// Euclidean distance.
    Eucl,
    /// Sandelin-Wasserman similarity.
    Sw,
    /// Symmetrized Kullback-Leibler divergence.
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
    fn needs_bkg(self) -> bool {
        matches!(
            self,
            ColumnMetric::Weucl | ColumnMetric::Wpcc | ColumnMetric::Allr | ColumnMetric::AllrLl
        )
    }

    /// The name [`str::parse`] reads, such as `"PCC"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            ColumnMetric::Pcc => "PCC",
            ColumnMetric::Eucl => "EUCL",
            ColumnMetric::Sw => "SW",
            ColumnMetric::Kl => "KL",
            ColumnMetric::Weucl => "WEUCL",
            ColumnMetric::Allr => "ALLR",
            ColumnMetric::Bhat => "BHAT",
            ColumnMetric::Hell => "HELL",
            ColumnMetric::Wpcc => "WPCC",
            ColumnMetric::Seucl => "SEUCL",
            ColumnMetric::Man => "MAN",
            ColumnMetric::AllrLl => "ALLR_LL",
        }
    }
}

impl std::str::FromStr for ColumnMetric {
    type Err = crate::error::Error;

    /// Parse a metric's short name, such as `"PCC"` or `"EUCL"`.
    fn from_str(s: &str) -> Result<ColumnMetric> {
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
}

impl std::fmt::Display for ColumnMetric {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// `w * ln(r)`, taking a zero weight to give zero even where the log is
/// infinite.
fn xlog(w: f64, r: f64) -> f64 {
    if w == 0.0 { 0.0 } else { w * r.ln() }
}

fn total(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(0.0, |a, b| a + b)
}

/// A motif column to compare with [`compare_columns`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Column<'a> {
    /// The probability of each letter.
    pub probs: &'a [f64],
    /// Background frequencies of the letters, needed by the weighted metrics
    /// and ALLR.
    pub bkg: &'a [f64],
    /// The number of sites the column was built from, needed by ALLR.
    pub nsites: f64,
}

impl<'a> Column<'a> {
    /// A column of probabilities, with no background and 100 sites.
    pub fn new(probs: &'a [f64]) -> Column<'a> {
        Column {
            probs,
            bkg: &[],
            nsites: 100.0,
        }
    }

    /// The column with background frequencies.
    pub fn with_bkg(self, bkg: &'a [f64]) -> Column<'a> {
        Column { bkg, ..self }
    }

    /// The column with a number of sites.
    pub fn with_nsites(self, nsites: f64) -> Column<'a> {
        Column { nsites, ..self }
    }
}

/// Compare two columns of probabilities over the same letters with `metric`.
/// The backgrounds are used by the weighted metrics and ALLR, which also uses
/// the site counts. A column whose first value is negative scores 0.
pub fn compare_columns(x: &Column, y: &Column, metric: ColumnMetric) -> Result<f64> {
    let (bkg1, bkg2, nsites1, nsites2) = (x.bkg, y.bkg, x.nsites, y.nsites);
    let (x, y) = (x.probs, y.probs);
    if x.len() != y.len() {
        return arg_err(format!(
            "x and y must have the same length, not {} and {}",
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
        // Each pair's two terms are added to the running sum one after the
        // other rather than summed first.
        ColumnMetric::Kl => {
            0.5 * pairs().fold(0.0, |s, (a, b)| s + xlog(*a, a / b) + xlog(*b, b / a))
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
            let wsum = total(bkg1.iter().zip(bkg2).map(|(a, b)| a + b));
            let w: Vec<f64> = bkg1.iter().zip(bkg2).map(|(a, b)| (a + b) / wsum).collect();
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
                    .map(|((a, b), g)| xlog(b * nsites2, a / g)),
            );
            let right = total(
                pairs()
                    .zip(bkg2)
                    .map(|((a, b), g)| xlog(a * nsites1, b / g)),
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

/// How [`average`] summarizes a set of values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Average {
    /// Arithmetic mean.
    #[default]
    ArithmeticMean,
    /// Geometric mean of the positive values, with the sum of their logarithms
    /// divided by the number of all values.
    GeometricMean,
    /// The median.
    Median,
    /// Fisher z-transform: the inverse hyperbolic tangent of the mean of
    /// the hyperbolic tangents.
    FisherZ,
}

impl Average {
    /// The name [`str::parse`] reads, such as `"a.mean"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            Average::ArithmeticMean => "a.mean",
            Average::GeometricMean => "g.mean",
            Average::Median => "median",
            Average::FisherZ => "fzt",
        }
    }
}

impl std::str::FromStr for Average {
    type Err = crate::error::Error;

    /// Parse `"a.mean"`, `"g.mean"`, `"median"` or `"fzt"`.
    fn from_str(s: &str) -> Result<Average> {
        match s {
            "a.mean" => Ok(Average::ArithmeticMean),
            "g.mean" => Ok(Average::GeometricMean),
            "median" => Ok(Average::Median),
            "fzt" => Ok(Average::FisherZ),
            _ => arg_err("average must be one of a.mean, g.mean, median, fzt"),
        }
    }
}

impl std::fmt::Display for Average {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Average `values` as chosen by `how`. The geometric mean is 0 when no value
/// is positive or the logarithms of the positive ones sum to 0.
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

/// The average information content of each motif's positions.
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
