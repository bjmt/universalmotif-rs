//! FASTA sequence files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use super::fmt_err;
use crate::error::{Error, Result};

/// Read FASTA records as `(name, sequence)` pairs. The name is the whole
/// header line after `>`, and whitespace inside sequences is removed. Lines
/// starting with `;` are comments.
pub fn read(text: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for line in text.lines() {
        add_line(&mut out, line)?;
    }
    Ok(out)
}

/// Read FASTA records from a reader line by line (see [`read`]), so the
/// file's text is never held in memory as a whole. The first `skip` lines
/// are ignored.
pub fn read_reader<R: BufRead>(mut reader: R, skip: usize) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let mut line = String::new();
    let mut n = 0;
    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|e| Error::Argument(format!("could not read the FASTA file: {e}")))?;
        if read == 0 {
            break;
        }
        n += 1;
        if n > skip {
            add_line(&mut out, line.trim_end_matches(['\n', '\r']))?;
        }
    }
    Ok(out)
}

/// Read a FASTA file (see [`read_reader`]).
pub fn read_path(path: impl AsRef<Path>, skip: usize) -> Result<Vec<(String, String)>> {
    let path = path.as_ref();
    let file = File::open(path)
        .map_err(|e| Error::Argument(format!("could not open {}: {e}", path.display())))?;
    read_reader(BufReader::with_capacity(1 << 20, file), skip)
}

fn add_line(out: &mut Vec<(String, String)>, line: &str) -> Result<()> {
    if let Some(name) = line.strip_prefix('>') {
        out.push((name.trim().to_string(), String::new()));
    } else if !(line.starts_with(';') || line.trim().is_empty()) {
        match out.last_mut() {
            Some((_, seq)) => {
                if line.is_ascii() {
                    // Whitespace is rare inside sequence lines; copy runs
                    // between it.
                    seq.extend(line.split_ascii_whitespace());
                } else {
                    seq.extend(line.chars().filter(|c| !c.is_whitespace()));
                }
            }
            None => return fmt_err("FASTA text must start with a '>' header line"),
        }
    }
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streaming_matches_text() {
        let text = "junk line\r\n>one  first\r\nACG T\r\n\r\n; comment\r\nacgt\n>two\nNNAC\tGT\n>empty\n>αβ\nαβγ\n";
        let skipped = text.split_once('\n').unwrap().1;
        let from_text = read(skipped).unwrap();
        let streamed = read_reader(text.as_bytes(), 1).unwrap();
        assert_eq!(streamed, from_text);
        assert_eq!(streamed[0], ("one  first".into(), "ACGTacgt".into()));
        assert_eq!(streamed[1].1, "NNACGT");
        assert_eq!(streamed[3].1, "αβγ");
        assert!(read_reader("ACGT\n".as_bytes(), 0).is_err());
    }
}
