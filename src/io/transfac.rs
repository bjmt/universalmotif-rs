//! TRANSFAC matrices (`read_transfac()` / `write_transfac()`).

use super::{content_lines, fmt_err, format_r_number, join_lines, num, parse_num};
use crate::alphabet::Alphabet;
use crate::convert::{convert_type, pcm_to_ppm};
use crate::error::{Error, Result};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};

/// Read TRANSFAC motifs.
///
/// Records are separated by `//` lines. The matrix follows a `P0` (or `PO`)
/// line, and the `AC`, `ID`, `NA`, `HC` and `OS` fields give the name,
/// alternate name, family and organism. The name is the accession if there is
/// one, with the identifier or factor name as the alternate name. Records
/// without a matrix, such as the `VV` header of TRANSFAC releases, are
/// skipped. Motifs are returned as PCMs; positions with fewer counts than the
/// largest are scaled up to it.
pub fn read(text: &str) -> Result<Vec<Motif>> {
    let lines: Vec<&str> = content_lines(text)
        .into_iter()
        .filter(|l| !l.starts_with("XX"))
        .collect();
    let seps: Vec<usize> = (0..lines.len())
        .filter(|&i| lines[i].starts_with("//"))
        .collect();
    if seps.is_empty() {
        return fmt_err("no motifs found (no '//' separators)");
    }
    let mut motifs = Vec::new();
    let mut start = 0;
    for &s in &seps {
        if let Some(m) = read_record(&lines[start..s])? {
            motifs.push(m);
        }
        start = s + 1;
    }
    if motifs.is_empty() {
        return fmt_err("no motifs found");
    }
    Ok(motifs)
}

fn read_record(lines: &[&str]) -> Result<Option<Motif>> {
    let Some(p0) = lines
        .iter()
        .position(|l| l.starts_with("P0") || l.starts_with("PO"))
    else {
        return Ok(None);
    };

    let mut fields: Vec<(&str, String)> = Vec::new();
    for l in &lines[..p0] {
        let (key, rest) = l.split_once(char::is_whitespace).unwrap_or((l, ""));
        if matches!(key, "AC" | "ID" | "NA" | "HC" | "OS") && !fields.iter().any(|(k, _)| *k == key)
        {
            fields.push((key, rest.trim().to_string()));
        }
    }
    let get = |k: &str| fields.iter().find(|(f, _)| *f == k).map(|(_, v)| v.clone());
    let (ac, id, na) = (get("AC"), get("ID"), get("NA"));
    let (name, altname) = match (&ac, &id, &na) {
        (Some(a), Some(_), _) => (Some(a.clone()), id.clone()),
        (Some(a), None, Some(_)) => (Some(a.clone()), na.clone()),
        (None, Some(i), Some(_)) => (Some(i.clone()), na.clone()),
        _ => (ac.clone().or(id.clone()).or(na.clone()), None),
    };

    // Column order from the P0 header, if it names the four bases.
    let header: Vec<char> = lines[p0]
        .split_whitespace()
        .skip(1)
        .filter_map(|t| {
            let mut c = t.chars();
            match (c.next(), c.next()) {
                (Some(l), None) => Some(l),
                _ => None,
            }
        })
        .collect();
    let mut order = [0usize, 1, 2, 3];
    let mut sorted = header.clone();
    sorted.sort_unstable();
    if header.len() == 4 && sorted == ['A', 'C', 'G', 'T'] {
        for (i, l) in ['A', 'C', 'G', 'T'].iter().enumerate() {
            order[i] = header.iter().position(|h| h == l).unwrap();
        }
    }

    let mut counts: Vec<[f64; 4]> = Vec::new();
    for l in &lines[p0 + 1..] {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.first().and_then(|v| parse_num(v)).is_none() {
            continue;
        }
        if t.len() < 5 {
            return fmt_err(format!(
                "TRANSFAC matrix row '{l}' has fewer than four counts"
            ));
        }
        let raw = [
            num(t[1], "a count")?,
            num(t[2], "a count")?,
            num(t[3], "a count")?,
            num(t[4], "a count")?,
        ];
        counts.push([raw[order[0]], raw[order[1]], raw[order[2]], raw[order[3]]]);
    }
    if counts.is_empty() {
        return fmt_err(format!(
            "empty motif matrix for motif '{}'",
            name.unwrap_or_default()
        ));
    }

    let nsites = counts
        .iter()
        .map(|c| c.iter().sum::<f64>())
        .fold(f64::NEG_INFINITY, f64::max);
    let mut data = Vec::with_capacity(counts.len() * 4);
    for c in &counts {
        data.extend(pcm_to_ppm(c, 0.0));
    }
    let matrix = Matrix::from_vec(4, counts.len(), data)?;
    let ppm = Motif::builder(matrix, Alphabet::Dna)
        .name(name.unwrap_or_default())
        .altname(altname)
        .family(get("HC"))
        .organism(get("OS"))
        .motif_type(Some(MotifType::Ppm))
        .nsites(Some(nsites))
        .build()?;
    convert_type(&ppm, MotifType::Pcm, None, false).map(Some)
}

/// Options for [`write()`].
#[derive(Debug, Clone, PartialEq)]
pub struct WriteOptions {
    /// Field used for the motif name.
    pub name_tag: String,
    /// Field used for the alternate name.
    pub altname_tag: String,
}

impl Default for WriteOptions {
    fn default() -> Self {
        WriteOptions {
            name_tag: "ID".into(),
            altname_tag: "NA".into(),
        }
    }
}

/// Write DNA or RNA motifs in TRANSFAC format, as count matrices.
pub fn write(motifs: &[Motif], opts: &WriteOptions) -> Result<String> {
    let mut lines = Vec::new();
    for m in motifs {
        if !m.alphabet.is_nucleotide() {
            return Err(Error::Unsupported(
                "TRANSFAC format only supports DNA and RNA motifs".into(),
            ));
        }
        let m = convert_type(m, MotifType::Pcm, None, false)?;
        lines.push(format!("{} {}", opts.name_tag, m.name));
        if let Some(a) = &m.altname {
            lines.push(format!("{} {}", opts.altname_tag, a));
        }
        if let Some(f) = &m.family {
            lines.push(format!("HC {f}"));
        }
        if let Some(o) = &m.organism {
            lines.push(format!("OS {o}"));
        }
        lines.push("P0\tA\tC\tG\tT".into());
        let consensus: Vec<char> = m.consensus.as_deref().unwrap_or("").chars().collect();
        for (j, c) in m.matrix.columns().enumerate() {
            let vals: Vec<String> = c.iter().map(|&v| format_r_number(v)).collect();
            let letter = consensus.get(j).map(|c| c.to_string()).unwrap_or_default();
            lines.push(format!("{:02}\t{}\t{}", j + 1, vals.join("\t"), letter));
        }
        lines.push("XX".into());
        lines.push("//".into());
    }
    Ok(join_lines(&lines))
}
