//! FASTA sequence files.

use super::fmt_err;
use crate::error::Result;

/// Read FASTA records as `(name, sequence)` pairs. The name is the whole
/// header line after `>`, and whitespace inside sequences is removed. Lines
/// starting with `;` are comments.
pub fn read(text: &str) -> Result<Vec<(String, String)>> {
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            out.push((name.trim().to_string(), String::new()));
        } else if line.starts_with(';') || line.trim().is_empty() {
            continue;
        } else {
            match out.last_mut() {
                Some((_, seq)) => seq.extend(line.chars().filter(|c| !c.is_whitespace())),
                None => return fmt_err("FASTA text must start with a '>' header line"),
            }
        }
    }
    Ok(out)
}

/// Write FASTA records, wrapping sequences at `width` characters (no
/// wrapping if `width` is 0).
pub fn write<S: AsRef<str>, T: AsRef<str>>(records: &[(S, T)], width: usize) -> String {
    let mut out = String::new();
    for (name, seq) in records {
        out.push('>');
        out.push_str(name.as_ref());
        out.push('\n');
        let seq = seq.as_ref();
        if width == 0 || seq.is_empty() {
            out.push_str(seq);
            out.push('\n');
        } else {
            for chunk in seq.as_bytes().chunks(width) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
        }
    }
    out
}
