//! HOMER motif files (`read_homer()` / `write_homer()`).

use super::{content_lines, fmt_err, format_fixed, format_r_number, join_lines, num, parse_num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

fn fields(line: &str) -> Vec<&str> {
    if line.contains('\t') {
        line.split('\t').map(str::trim).collect()
    } else {
        line.split_whitespace().collect()
    }
}

/// Read HOMER motifs.
///
/// The header line gives the consensus, the name, the log-odds detection
/// threshold (kept in `extrainfo` as `logodds.threshold`) and, in files written by
/// HOMER itself, the target and background site counts and the P-value. A
/// name such as `RUNX(Runt)/HPC7-Runx1-ChIP-Seq/Homer` is split into the name
/// `RUNX` and the family `Runt`. Rows are DNA probabilities.
pub fn read(text: &str) -> Result<Vec<Motif>> {
    let lines = content_lines(text);
    let headers: Vec<usize> = (0..lines.len())
        .filter(|&i| lines[i].starts_with('>'))
        .collect();
    if headers.is_empty() {
        return fmt_err("no motifs found");
    }
    let mut motifs = Vec::with_capacity(headers.len());
    for (k, &h) in headers.iter().enumerate() {
        let end = headers.get(k + 1).copied().unwrap_or(lines.len());
        let f = fields(lines[h]);

        let mut name = f.get(1).copied().unwrap_or("").to_string();
        let mut family = None;
        if name.contains('/') {
            let first = name.split('/').next().unwrap_or("").to_string();
            let mut parts = first.splitn(2, '(');
            name = parts.next().unwrap_or("").to_string();
            family = parts
                .next()
                .map(|rest| rest.split(')').next().unwrap_or("").to_string());
        }

        let (mut nsites, mut bkgsites, mut pval) = (None, None, None);
        if let Some(stats) = f.get(5) {
            for part in stats.split(',') {
                let before_paren = |s: &str| s.split('(').next().unwrap_or("").to_string();
                if let Some(v) = part.strip_prefix("T:") {
                    nsites = parse_num(&before_paren(v));
                } else if let Some(v) = part.strip_prefix("B:") {
                    bkgsites = parse_num(&before_paren(v));
                } else if let Some(v) = part.strip_prefix("P:") {
                    pval = parse_num(v);
                }
            }
        }
        let extrainfo = f
            .get(2)
            .map(|t| vec![("logodds.threshold".to_string(), t.to_string())])
            .unwrap_or_default();

        let mut data = Vec::new();
        for row in &lines[h + 1..end] {
            let vals: Vec<f64> = row
                .split_whitespace()
                .map(|v| num(v, "a probability"))
                .collect::<Result<_>>()?;
            if vals.len() != 4 {
                return fmt_err(format!(
                    "motif '{name}': expected 4 values per row, found {}",
                    vals.len()
                ));
            }
            data.extend(vals);
        }
        if data.is_empty() {
            return fmt_err(format!("empty motif matrix for motif '{name}'"));
        }
        let matrix = Matrix::from_vec(4, data.len() / 4, data)?;
        motifs.push(
            Motif::builder(matrix, Alphabet::Dna)
                .name(name)
                .family(family)
                .motif_type(Some(MotifType::Ppm))
                .nsites(nsites)
                .bkgsites(bkgsites)
                .pval(pval)
                .extrainfo(extrainfo)
                .build()?,
        );
    }
    Ok(motifs)
}

/// How the log-odds detection threshold written into HOMER headers is set.
#[derive(Debug, Clone, PartialEq)]
pub enum Threshold {
    /// A fraction of each motif's maximum possible log-odds score, either one
    /// value for all motifs or one per motif.
    Logodds(Vec<f64>),
    /// Absolute log-odds scores, either one value or one per motif.
    LogoddsAbs(Vec<f64>),
}

impl Default for Threshold {
    fn default() -> Self {
        Threshold::Logodds(vec![0.8])
    }
}

/// Write DNA or RNA motifs in HOMER format.
pub fn write(motifs: &[Motif], threshold: &Threshold) -> Result<String> {
    let values = match threshold {
        Threshold::Logodds(v) | Threshold::LogoddsAbs(v) => v,
    };
    if values.len() != 1 && values.len() != motifs.len() {
        return fmt_err("a single threshold or one per motif should be provided");
    }
    let mut lines = Vec::new();
    for (i, m) in motifs.iter().enumerate() {
        if !m.alphabet.is_nucleotide() {
            return Err(Error::Unsupported(
                "HOMER format only supports DNA and RNA motifs".into(),
            ));
        }
        let t = values[if values.len() == 1 { 0 } else { i }];
        let thr = match threshold {
            Threshold::LogoddsAbs(_) => t,
            Threshold::Logodds(_) => {
                let pwm = convert_type(m, MotifType::Pwm, None, false)?;
                let max: f64 = pwm
                    .matrix
                    .columns()
                    .map(|c| c.iter().copied().fold(f64::NEG_INFINITY, f64::max))
                    .sum();
                max * t
            }
        };
        let ppm = convert_type(m, MotifType::Ppm, None, false)?;
        lines.push(format!(
            ">{}\t{}\t{}",
            m.consensus.as_deref().unwrap_or(""),
            m.name,
            format_r_number(thr)
        ));
        for c in ppm.matrix.columns() {
            let vals: Vec<String> = c.iter().map(|&v| format_fixed(v, 3)).collect();
            lines.push(vals.join("\t"));
        }
    }
    Ok(join_lines(&lines))
}
