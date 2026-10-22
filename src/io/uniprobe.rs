//! UniPROBE motif files.

use super::{content_lines, fmt_err, format_shortest, join_lines, num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

fn is_motif_header(line: &str) -> bool {
    line.find('#').is_some_and(|i| {
        let rest = &line[i + 1..];
        rest.starts_with(char::is_whitespace) && rest.trim_start().starts_with("Motif:")
    })
}

/// Four rows labelled `A:`, `C:`, `G:` and `T:` starting at `lines[i]`.
fn read_rows(lines: &[&str], i: usize, name: &str) -> Result<Matrix> {
    let mut rows = Vec::with_capacity(4);
    for (k, label) in ["A:", "C:", "G:", "T:"].iter().enumerate() {
        let line = lines.get(i + k).copied().unwrap_or("");
        let mut t = line.split_whitespace();
        if t.next() != Some(label) {
            return fmt_err(format!(
                "motif '{name}': expected a row starting with '{label}'"
            ));
        }
        rows.push(
            t.map(|v| num(v, "a probability"))
                .collect::<Result<Vec<f64>>>()?,
        );
    }
    let width = rows[0].len();
    if width == 0 || rows.iter().any(|r| r.len() != width) {
        return fmt_err(format!("motif '{name}' has rows of different lengths"));
    }
    let cols: Vec<Vec<f64>> = (0..width)
        .map(|j| rows.iter().map(|r| r[j]).collect())
        .collect();
    Matrix::from_columns(&cols)
}

/// Read UniPROBE motifs, either from full output files with `# Motif:` headers,
/// where the probability matrix is used and the enrichment score is kept in
/// `extrainfo`, or from minimal files of a name line followed by `A:`, `C:`,
/// `G:` and `T:` rows.
pub fn read(text: &str) -> Result<Vec<Motif>> {
    let lines = content_lines(text);
    let mut motifs = Vec::new();
    if lines.iter().any(|l| is_motif_header(l)) {
        let metas: Vec<Vec<&str>> = (0..lines.len())
            .filter(|&i| is_motif_header(lines[i]))
            .map(|i| {
                lines
                    .get(i + 1)
                    .map(|l| l.split_whitespace().skip(1).collect())
                    .unwrap_or_default()
            })
            .collect();
        let starts: Vec<usize> = (0..lines.len())
            .filter(|&i| lines[i].starts_with("Probability matrix"))
            .map(|i| i + 1)
            .collect();
        if metas.len() != starts.len() {
            return fmt_err("every UniPROBE motif needs a probability matrix");
        }
        for (meta, start) in metas.iter().zip(starts) {
            let name = meta.first().copied().unwrap_or("");
            let matrix = read_rows(&lines, start, name)?;
            let extrainfo = meta
                .get(1)
                .map(|s| vec![("enrichment.score".to_string(), s.to_string())])
                .unwrap_or_default();
            motifs.push(
                Motif::builder(matrix, Alphabet::Dna)
                    .name(name)
                    .motif_type(Some(MotifType::Ppm))
                    .extrainfo(extrainfo)
                    .build()?,
            );
        }
    } else {
        for i in 0..lines.len() {
            if !lines[i].starts_with("A:") {
                continue;
            }
            let name = if i > 0 { lines[i - 1].trim() } else { "" };
            let matrix = read_rows(&lines, i, name)?;
            motifs.push(
                Motif::builder(matrix, Alphabet::Dna)
                    .name(name)
                    .motif_type(Some(MotifType::Ppm))
                    .build()?,
            );
        }
    }
    if motifs.is_empty() {
        return fmt_err("no motifs found");
    }
    Ok(motifs)
}

/// Write DNA motifs in the minimal UniPROBE format: a name line followed by
/// `A:`, `C:`, `G:` and `T:` rows of probabilities, written with full
/// precision.
pub fn write(motifs: &[Motif]) -> Result<String> {
    let mut lines = Vec::new();
    for m in motifs {
        if m.alphabet != Alphabet::Dna {
            return Err(Error::Unsupported(
                "UniPROBE format only supports DNA motifs".into(),
            ));
        }
        let m = convert_type(m, MotifType::Ppm, None, false)?;
        lines.push(m.name.clone());
        for (i, l) in m.alphabet.letters().iter().enumerate() {
            let vals: Vec<String> = (0..m.width())
                .map(|j| format_shortest(m.matrix.get(i, j)))
                .collect();
            lines.push(format!("{l}:\t{}", vals.join("\t")));
        }
        lines.push(String::new());
    }
    Ok(join_lines(&lines))
}
