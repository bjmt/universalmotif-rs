//! JASPAR count matrices.

use super::{content_lines, fmt_err, join_lines, num};
use crate::alphabet::Alphabet;
use crate::convert::convert_type;
use crate::create::guess_matrix_alphabet;
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

/// Read JASPAR-formatted motifs.
///
/// Each motif starts with a `>name altname` header followed by one row per
/// letter, either labelled as in `A [ 1 2 3 ]` or as bare numbers. Labelled
/// rows may come in any order and the alphabet is taken from the labels; four
/// bare rows are read as DNA and twenty as amino acids. Motifs are read as PCMs
/// with `nsites` set to the largest column sum.
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
        let mut names = lines[h][1..].split_whitespace();
        let name = names.next().unwrap_or("").to_string();
        let altname = names.next().map(String::from);
        motifs.push(read_one(&lines[h + 1..end], name, altname)?);
    }
    Ok(motifs)
}

fn read_one(rows: &[&str], name: String, altname: Option<String>) -> Result<Motif> {
    let empty = || Error::Argument(format!("empty motif matrix for motif '{name}'"));
    let tokens: Vec<Vec<String>> = rows
        .iter()
        .map(|r| {
            r.replacen('[', " ", 1)
                .replacen(']', " ", 1)
                .split_whitespace()
                .map(String::from)
                .collect()
        })
        .collect();
    let labelled = tokens.iter().any(|t| {
        t.first()
            .is_some_and(|f| f.len() == 1 && f.chars().all(|c| c.is_ascii_uppercase()))
    });

    let mut letters = Vec::new();
    let mut values: Vec<Vec<f64>> = Vec::new();
    for t in &tokens {
        let t = if labelled {
            let Some((label, rest)) = t.split_first() else {
                continue;
            };
            let mut c = label.chars();
            match (c.next(), c.next()) {
                (Some(l), None) => letters.push(l),
                _ => return fmt_err(format!("unexpected row label '{label}' in motif '{name}'")),
            }
            rest
        } else {
            &t[..]
        };
        values.push(t.iter().map(|v| num(v, "a count")).collect::<Result<_>>()?);
    }
    let width = values.first().map(Vec::len).unwrap_or(0);
    if values.is_empty() || width == 0 {
        return Err(empty());
    }
    if values.iter().any(|r| r.len() != width) {
        return fmt_err(format!("rows of motif '{name}' have different lengths"));
    }

    let (alphabet, order): (Alphabet, Vec<usize>) = if labelled {
        let mut order: Vec<usize> = (0..letters.len()).collect();
        order.sort_by_key(|&i| letters[i]);
        let sorted: Vec<char> = order.iter().map(|&i| letters[i]).collect();
        (guess_matrix_alphabet(&sorted)?, order)
    } else {
        let a = match values.len() {
            4 => Alphabet::Dna,
            20 => Alphabet::Aa,
            n => {
                return fmt_err(format!(
                    "motif '{name}' has {n} unlabelled rows; label the rows to give its alphabet"
                ));
            }
        };
        (a, (0..values.len()).collect())
    };

    let data: Vec<f64> = (0..width)
        .flat_map(|j| order.iter().map(|&i| values[i][j]).collect::<Vec<_>>())
        .collect();
    let matrix = Matrix::from_vec(order.len(), width, data)?;
    let nsites = matrix
        .col_sums()
        .into_iter()
        .fold(f64::NEG_INFINITY, f64::max);
    Motif::builder(matrix, alphabet)
        .name(name)
        .altname(altname)
        .motif_type(Some(MotifType::Pcm))
        .nsites(Some(nsites))
        .build()
}

/// Write motifs in JASPAR format, as count matrices. Only DNA, RNA and amino
/// acid motifs can be written.
pub fn write(motifs: &[Motif]) -> Result<String> {
    let mut lines = Vec::new();
    for m in motifs {
        if matches!(m.alphabet, Alphabet::Custom(_)) {
            return Err(Error::Unsupported(
                "JASPAR format only supports DNA, RNA and AA motifs".into(),
            ));
        }
        let m = convert_type(m, MotifType::Pcm, None, false)?;
        lines.push(match &m.altname {
            Some(a) => format!(">{} {}", m.name, a),
            None => format!(">{}", m.name),
        });
        let nsites = m.nsites.unwrap_or(100.0);
        let width = if nsites > 99999.0 {
            6
        } else if nsites > 9999.0 {
            5
        } else if nsites > 999.0 {
            4
        } else if nsites > 99.0 {
            3
        } else {
            2
        };
        for (i, l) in m.alphabet.letters().iter().enumerate() {
            let vals: Vec<String> = (0..m.width())
                .map(|j| format!("{:>width$.0}", m.matrix.get(i, j)))
                .collect();
            lines.push(format!("{l} [ {} ]", vals.join(" ")));
        }
        lines.push(String::new());
    }
    Ok(join_lines(&lines))
}
