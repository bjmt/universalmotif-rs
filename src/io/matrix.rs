//! Plain numeric matrices.

use super::{fmt_err, format_fixed, join_lines, num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::Result;
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

/// Whether motif positions are the columns or the rows of the matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Positions {
    #[default]
    Columns,
    Rows,
}

/// How motif names are given.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Headers {
    /// No names.
    None,
    /// The first line of each motif is its name.
    #[default]
    FirstLine,
    /// Lines containing this text are headers; the name is the line with the
    /// text removed. Each header starts a new motif.
    Prefix(String),
}

/// Options for [`read()`].
#[derive(Debug, Clone, PartialEq)]
pub struct ReadOptions {
    pub positions: Positions,
    pub alphabet: Alphabet,
    /// The matrix type; detected from the values when `None`.
    pub motif_type: Option<MotifType>,
    /// Line separating motifs, where an empty string means blank lines. `None`
    /// means motifs are only separated by headers.
    pub sep: Option<String>,
    pub headers: Headers,
    /// Whether each row holds a letter label, or with `Positions::Rows` whether
    /// the first line holds them. The labels are skipped; rows are read in
    /// alphabet order.
    pub rownames: bool,
    /// Lines starting with this text are ignored.
    pub comment: Option<String>,
}

impl Default for ReadOptions {
    fn default() -> Self {
        ReadOptions {
            positions: Positions::Columns,
            alphabet: Alphabet::Dna,
            motif_type: None,
            sep: Some(String::new()),
            headers: Headers::FirstLine,
            rownames: false,
            comment: None,
        }
    }
}

fn is_sep(line: &str, sep: &Option<String>) -> bool {
    match sep {
        Some(s) if s.is_empty() => line.trim().is_empty(),
        Some(s) => line == s,
        None => false,
    }
}

/// Read motifs stored as plain matrices of numbers.
pub fn read(text: &str, opts: &ReadOptions) -> Result<Vec<Motif>> {
    let lines: Vec<&str> = text
        .lines()
        .filter(|l| {
            opts.comment
                .as_ref()
                .is_none_or(|c| !l.starts_with(c.as_str()))
        })
        .collect();

    // Split into (name, rows) blocks.
    let mut blocks: Vec<(Option<String>, Vec<&str>)> = Vec::new();
    match &opts.headers {
        Headers::Prefix(p) => {
            for l in &lines {
                if let Some(i) = l.find(p.as_str()) {
                    let name = format!("{}{}", &l[..i], &l[i + p.len()..]);
                    blocks.push((Some(name), Vec::new()));
                } else if !l.trim().is_empty()
                    && !is_sep(l, &opts.sep)
                    && let Some(b) = blocks.last_mut()
                {
                    b.1.push(l);
                }
            }
        }
        Headers::FirstLine | Headers::None => {
            let mut cur: Vec<&str> = Vec::new();
            for l in lines.iter().chain(std::iter::once(&"")) {
                let end = is_sep(l, &opts.sep) || (opts.sep.is_none() && l.is_empty());
                if end {
                    if !cur.is_empty() {
                        let block = std::mem::take(&mut cur);
                        if opts.headers == Headers::FirstLine {
                            let (name, rows) = block.split_first().unwrap();
                            blocks.push((Some(name.to_string()), rows.to_vec()));
                        } else {
                            blocks.push((None, block));
                        }
                    }
                } else if !l.trim().is_empty() {
                    cur.push(l);
                }
            }
        }
    }
    if blocks.is_empty() {
        return fmt_err("no motifs found");
    }

    let mut motifs = Vec::with_capacity(blocks.len());
    for (name, rows) in blocks {
        let mut rows: &[&str] = &rows;
        if opts.rownames && opts.positions == Positions::Rows {
            rows = rows.get(1..).unwrap_or(&[]);
        }
        let skip = usize::from(opts.rownames && opts.positions == Positions::Columns);
        let table: Vec<Vec<f64>> = rows
            .iter()
            .map(|r| {
                r.split_whitespace()
                    .skip(skip)
                    .map(|v| num(v, "a matrix value"))
                    .collect()
            })
            .collect::<Result<_>>()?;
        let ncol = table.first().map(Vec::len).unwrap_or(0);
        if table.is_empty() || ncol == 0 || table.iter().any(|r| r.len() != ncol) {
            return fmt_err(format!(
                "motif '{}' has an empty or ragged matrix",
                name.as_deref().unwrap_or("")
            ));
        }
        let matrix = match opts.positions {
            Positions::Rows => Matrix::from_columns(&table)?,
            Positions::Columns => {
                let cols: Vec<Vec<f64>> = (0..ncol)
                    .map(|j| table.iter().map(|r| r[j]).collect())
                    .collect();
                Matrix::from_columns(&cols)?
            }
        };
        let mut b = Motif::builder(matrix, opts.alphabet.clone()).motif_type(opts.motif_type);
        if let Some(n) = name {
            b = b.name(n);
        }
        motifs.push(b.build()?);
    }
    Ok(motifs)
}

/// Options for [`write()`].
#[derive(Debug, Clone, PartialEq)]
pub struct WriteOptions {
    pub positions: Positions,
    /// Write letter labels.
    pub rownames: bool,
    /// Convert motifs to this type first.
    pub motif_type: Option<MotifType>,
    /// Line written after each motif; `None` writes nothing.
    pub sep: Option<String>,
    pub headers: Headers,
    /// Decimal places.
    pub digits: usize,
}

impl Default for WriteOptions {
    fn default() -> Self {
        WriteOptions {
            positions: Positions::Columns,
            rownames: false,
            motif_type: None,
            sep: Some(String::new()),
            headers: Headers::FirstLine,
            digits: 6,
        }
    }
}

/// Write motifs as plain matrices of numbers.
pub fn write(motifs: &[Motif], opts: &WriteOptions) -> Result<String> {
    let mut lines = Vec::new();
    for m in motifs {
        let m = match opts.motif_type {
            Some(t) => convert_type(m, t, None, false)?,
            None => m.clone(),
        };
        match &opts.headers {
            Headers::None => {}
            Headers::FirstLine => lines.push(m.name.clone()),
            Headers::Prefix(p) => lines.push(format!("{p}{}", m.name)),
        }
        let letters = m.alphabet.letters();
        match opts.positions {
            Positions::Columns => {
                for (i, l) in letters.iter().enumerate() {
                    let mut vals: Vec<String> = (0..m.width())
                        .map(|j| format_fixed(m.matrix.get(i, j), opts.digits))
                        .collect();
                    if opts.rownames {
                        vals.insert(0, l.to_string());
                    }
                    lines.push(vals.join("\t"));
                }
            }
            Positions::Rows => {
                if opts.rownames {
                    let l: Vec<String> = letters.iter().map(|c| c.to_string()).collect();
                    lines.push(l.join("\t"));
                }
                for c in m.matrix.columns() {
                    let vals: Vec<String> =
                        c.iter().map(|&v| format_fixed(v, opts.digits)).collect();
                    lines.push(vals.join("\t"));
                }
            }
        }
        if let Some(s) = &opts.sep {
            lines.push(s.clone());
        }
    }
    Ok(join_lines(&lines))
}
