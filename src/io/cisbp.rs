//! CIS-BP motif files (`read_cisbp()`, plus a writer universalmotif lacks).

use super::{content_lines, fmt_err, format_shortest, join_lines, num, parse_num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

#[derive(Default, Clone)]
struct Meta {
    name: Option<String>,
    altname: Option<String>,
    family: Option<String>,
    organism: Option<String>,
}

struct Pending {
    meta: Meta,
    letters: Vec<char>,
    rows: Vec<Vec<f64>>,
}

/// The value after a key such as `TF` or `TF Name`, if `line` starts with
/// those words.
fn value_after<'a>(line: &'a str, key: &[&str]) -> Option<&'a str> {
    let mut rest = line;
    for k in key {
        rest = rest.trim_start().strip_prefix(k)?;
        if !rest.is_empty() && !rest.starts_with(char::is_whitespace) {
            return None;
        }
    }
    Some(rest.trim())
}

/// Read CIS-BP motifs.
///
/// Each matrix starts with a `Pos` header naming the letters (DNA or RNA)
/// followed by one row per position. The `TF`, `TF Name`, `Family` and
/// `Species` lines before it give the name, alternate name, family and
/// organism.
pub fn read(text: &str) -> Result<Vec<Motif>> {
    let lines = content_lines(text);
    if !lines.iter().any(|l| l.starts_with("Pos")) {
        return fmt_err("no motifs found");
    }
    let mut motifs = Vec::new();
    let mut meta = Meta::default();
    let mut current: Option<Pending> = None;

    for l in lines {
        let first = l.split_whitespace().next().unwrap_or("");
        if l.starts_with("Pos") {
            if let Some(p) = current.take() {
                finish(p, &mut motifs)?;
            }
            let letters = l
                .split_whitespace()
                .skip(1)
                .map(|t| {
                    let mut c = t.chars();
                    match (c.next(), c.next()) {
                        (Some(x), None) => Ok(x),
                        _ => fmt_err(format!("unexpected letter '{t}' in CIS-BP header")),
                    }
                })
                .collect::<Result<_>>()?;
            current = Some(Pending {
                meta: std::mem::take(&mut meta),
                letters,
                rows: Vec::new(),
            });
        } else if parse_num(first).is_some() {
            if let Some(p) = current.as_mut() {
                let vals: Vec<f64> = l
                    .split_whitespace()
                    .skip(1)
                    .map(|v| num(v, "a probability"))
                    .collect::<Result<_>>()?;
                if vals.len() != p.letters.len() {
                    return fmt_err(format!("CIS-BP row '{l}' does not match its header"));
                }
                p.rows.push(vals);
            }
        } else {
            if let Some(p) = current.take() {
                finish(p, &mut motifs)?;
            }
            if let Some(v) = value_after(l, &["TF", "Name"]) {
                meta.altname = Some(v.to_string());
            } else if let Some(v) = value_after(l, &["TF"]) {
                meta.name = Some(v.to_string());
            } else if let Some(v) = value_after(l, &["Family"]) {
                meta.family = Some(v.to_string());
            } else if let Some(v) = value_after(l, &["Species"]) {
                meta.organism = Some(v.to_string());
            }
        }
    }
    if let Some(p) = current.take() {
        finish(p, &mut motifs)?;
    }
    Ok(motifs)
}

fn finish(p: Pending, motifs: &mut Vec<Motif>) -> Result<()> {
    if p.rows.is_empty() {
        return Ok(());
    }
    let mut order: Vec<usize> = (0..p.letters.len()).collect();
    order.sort_by_key(|&i| p.letters[i]);
    let sorted: Vec<char> = order.iter().map(|&i| p.letters[i]).collect();
    let alphabet = match sorted.as_slice() {
        ['A', 'C', 'G', 'T'] => Alphabet::Dna,
        ['A', 'C', 'G', 'U'] => Alphabet::Rna,
        _ => return fmt_err("CIS-BP motifs must use the letters ACGT or ACGU"),
    };
    let mut data = Vec::with_capacity(p.rows.len() * 4);
    for r in &p.rows {
        data.extend(order.iter().map(|&i| r[i]));
    }
    let matrix = Matrix::from_vec(4, p.rows.len(), data)?;
    let mut b = Motif::builder(matrix, alphabet)
        .motif_type(Some(MotifType::Ppm))
        .altname(p.meta.altname)
        .family(p.meta.family)
        .organism(p.meta.organism);
    if let Some(n) = p.meta.name {
        b = b.name(n);
    }
    motifs.push(b.build()?);
    Ok(())
}

/// Write DNA or RNA motifs in CIS-BP format, as probability matrices. The
/// name, alternate name, family and organism go in the `TF`, `TF Name`,
/// `Family` and `Species` fields. Values are written with full precision.
pub fn write(motifs: &[Motif]) -> Result<String> {
    let mut lines = Vec::new();
    for m in motifs {
        if !m.alphabet.is_nucleotide() {
            return Err(Error::Unsupported(
                "CIS-BP format only supports DNA and RNA motifs".into(),
            ));
        }
        let m = convert_type(m, MotifType::Ppm, None, false)?;
        lines.push(format!("TF\t{}", m.name));
        if let Some(a) = &m.altname {
            lines.push(format!("TF Name\t{a}"));
        }
        if let Some(f) = &m.family {
            lines.push(format!("Family\t{f}"));
        }
        if let Some(o) = &m.organism {
            lines.push(format!("Species\t{o}"));
        }
        let letters: Vec<String> = m.alphabet.letters().iter().map(|c| c.to_string()).collect();
        lines.push(format!("Pos\t{}", letters.join("\t")));
        for (j, c) in m.matrix.columns().enumerate() {
            let vals: Vec<String> = c.iter().map(|&v| format_shortest(v)).collect();
            lines.push(format!("{}\t{}", j + 1, vals.join("\t")));
        }
        lines.push(String::new());
    }
    Ok(join_lines(&lines))
}
