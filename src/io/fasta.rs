//! FASTA sequence files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use super::fmt_err;
use crate::alphabet::Alphabet;
use crate::error::{Error, Result};
use crate::packed::{EncodedSequences, RecordEncoder};

/// Read FASTA records as `(name, sequence)` pairs. The name is the whole
/// header line after `>`, and whitespace inside sequences is removed. Lines
/// starting with `;` are comments.
pub fn read(text: &str) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for line in text.lines() {
        parse_line(line, &mut |l| collect(&mut out, l))?;
    }
    Ok(out)
}

/// Read FASTA records from a reader line by line (see [`read`]), so the
/// file's text is never held in memory as a whole. The first `skip` lines
/// are ignored.
pub fn read_reader<R: BufRead>(reader: R, skip: usize) -> Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    for_each_line(reader, skip, |l| collect(&mut out, l))?;
    Ok(out)
}

/// Read a FASTA file (see [`read_reader`]), which may be gzipped.
pub fn read_path(path: impl AsRef<Path>, skip: usize) -> Result<Vec<(String, String)>> {
    read_reader(open(path)?, skip)
}

/// Read a FASTA file (which may be gzipped) straight into
/// [`EncodedSequences`] for `alphabet` (see [`read_reader`]), without ever
/// holding the sequences as text.
/// Returns the names and the sequences.
pub fn read_encoded(
    path: impl AsRef<Path>,
    alphabet: &Alphabet,
    skip: usize,
) -> Result<(Vec<String>, EncodedSequences)> {
    let mut out = (Vec::new(), EncodedSequences::empty(alphabet)?);
    read_encoded_batches(open(path)?, alphabet, skip, usize::MAX, |names, seqs| {
        out = (names, seqs);
        Ok(())
    })?;
    out.1.shrink_to_fit();
    Ok(out)
}

/// Read FASTA records from a reader into [`EncodedSequences`] for
/// `alphabet`, handing them to `f` in batches of whole records holding at
/// least `batch_size` letters (the last batch may hold fewer). Only the
/// batch being filled is held in memory.
pub fn read_encoded_batches<R: BufRead>(
    reader: R,
    alphabet: &Alphabet,
    skip: usize,
    batch_size: usize,
    f: impl FnMut(Vec<String>, EncodedSequences) -> Result<()>,
) -> Result<()> {
    read_encoded_batches_reusing(reader, alphabet, skip, batch_size, || None, f)
}

/// [`read_encoded_batches`], storing each record's letters in a buffer from
/// `spare` when it has one.
pub(crate) fn read_encoded_batches_reusing<R: BufRead>(
    reader: R,
    alphabet: &Alphabet,
    skip: usize,
    batch_size: usize,
    mut spare: impl FnMut() -> Option<Vec<u8>>,
    mut f: impl FnMut(Vec<String>, EncodedSequences) -> Result<()>,
) -> Result<()> {
    let mut names: Vec<String> = Vec::new();
    let mut batch = EncodedSequences::empty(alphabet)?;
    let mut current: Option<RecordEncoder> = None;
    let mut size = 0;
    for_each_line(reader, skip, |line| {
        match line {
            Line::Header(name) => {
                if let Some(r) = current.take() {
                    size += r.len();
                    batch.push(r);
                    if size >= batch_size {
                        let full =
                            std::mem::replace(&mut batch, EncodedSequences::empty(alphabet)?);
                        f(std::mem::take(&mut names), full)?;
                        size = 0;
                    }
                }
                names.push(name.to_string());
                current = Some(batch.record_reusing(spare().unwrap_or_default()));
            }
            // for_each_line reports sequence before any header as an error.
            Line::Sequence(piece) => current.as_mut().unwrap().push(piece),
        }
        Ok(())
    })?;
    if let Some(r) = current.take() {
        batch.push(r);
    }
    if !names.is_empty() {
        f(names, batch)?;
    }
    Ok(())
}

/// Open a FASTA file for reading, decompressing it if it is gzipped
/// (bgzip included), which is recognised from its first bytes rather than
/// its name.
pub(crate) fn open(path: impl AsRef<Path>) -> Result<Box<dyn BufRead + Send>> {
    let path = path.as_ref();
    let file = File::open(path)
        .map_err(|e| Error::Argument(format!("could not open {}: {e}", path.display())))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let head = reader
        .fill_buf()
        .map_err(|e| Error::Argument(format!("could not read {}: {e}", path.display())))?;
    if head.starts_with(&[0x1f, 0x8b]) {
        // bgzip files are many gzip members one after another.
        let gz = MultiGzDecoder::new(reader);
        Ok(Box::new(BufReader::with_capacity(1 << 20, gz)))
    } else {
        Ok(Box::new(reader))
    }
}

/// A line of a FASTA file: a record's name, or a piece of its sequence.
enum Line<'a> {
    Header(&'a str),
    Sequence(&'a [u8]),
}

fn collect(out: &mut Vec<(String, String)>, line: Line<'_>) -> Result<()> {
    match line {
        Line::Header(name) => out.push((name.to_string(), String::new())),
        // parse_line only passes valid UTF-8, split at whitespace.
        Line::Sequence(piece) => match out.last_mut() {
            Some((_, seq)) => seq.push_str(std::str::from_utf8(piece).unwrap()),
            None => return fmt_err("FASTA text must start with a '>' header line"),
        },
    }
    Ok(())
}

/// Feed the lines of FASTA text from `reader` (after the first `skip`) to
/// `f` as headers and pieces of sequence.
fn for_each_line<R: BufRead>(
    mut reader: R,
    skip: usize,
    mut f: impl FnMut(Line<'_>) -> Result<()>,
) -> Result<()> {
    let mut line = String::new();
    let mut n = 0;
    let mut started = false;
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
            parse_line(line.trim_end_matches(['\n', '\r']), &mut |l| {
                match l {
                    Line::Header(_) => started = true,
                    Line::Sequence(_) if !started => {
                        return fmt_err("FASTA text must start with a '>' header line");
                    }
                    Line::Sequence(_) => {}
                }
                f(l)
            })?;
        }
    }
    Ok(())
}

fn parse_line(line: &str, f: &mut impl FnMut(Line<'_>) -> Result<()>) -> Result<()> {
    if let Some(name) = line.strip_prefix('>') {
        f(Line::Header(name.trim()))
    } else if line.starts_with(';') || line.trim().is_empty() {
        Ok(())
    } else if line.is_ascii() {
        // Whitespace is rare inside sequence lines; pass on the runs
        // between it.
        line.split_ascii_whitespace()
            .try_for_each(|piece| f(Line::Sequence(piece.as_bytes())))
    } else {
        let kept: String = line.chars().filter(|c| !c.is_whitespace()).collect();
        f(Line::Sequence(kept.as_bytes()))
    }
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
pub(crate) mod tests {
    use super::*;

    /// `text` gzipped as two members split at `split`, as bgzip writes
    /// blocks.
    pub(crate) fn gzip_members(text: &[u8], split: usize) -> Vec<u8> {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::io::Write;
        let mut out = Vec::new();
        for part in [&text[..split], &text[split..]] {
            let mut e = GzEncoder::new(Vec::new(), Compression::fast());
            e.write_all(part).unwrap();
            out.extend(e.finish().unwrap());
        }
        out
    }

    #[test]
    fn gzipped_files_read_like_plain_ones() {
        let text = ">a x\nACGTN\nacgt\n>b\nGGCC\n>c\nTTTT\n";
        let dir = std::env::temp_dir().join(format!("um-fasta-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (plain, gz) = (dir.join("x.fa"), dir.join("x.fa.gz"));
        std::fs::write(&plain, text).unwrap();
        std::fs::write(&gz, gzip_members(text.as_bytes(), 15)).unwrap();
        assert_eq!(read_path(&gz, 0).unwrap(), read(text).unwrap());
        assert_eq!(read_path(&gz, 3).unwrap(), read_path(&plain, 3).unwrap());
        assert_eq!(
            read_encoded(&gz, &Alphabet::Dna, 0).unwrap(),
            read_encoded(&plain, &Alphabet::Dna, 0).unwrap()
        );
        // A gzip file cut short is an error, not a short read.
        let whole = gzip_members(text.as_bytes(), 15);
        std::fs::write(&gz, &whole[..whole.len() - 6]).unwrap();
        assert!(read_path(&gz, 0).is_err());
        std::fs::write(&gz, b"").unwrap();
        assert!(read_path(&gz, 0).unwrap().is_empty());
        std::fs::remove_dir_all(&dir).unwrap();
    }

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
        assert!(read("ACGT\n").is_err());
    }

    #[test]
    fn encoded_batches_match_text() {
        let text = ">a\nACGTN\nacgt\n>b\n\n>c x\nGG CC\n>d\nTTTT\n";
        let records = read(text).unwrap();
        let seqs: Vec<&str> = records.iter().map(|r| r.1.as_str()).collect();
        let whole = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
        for batch_size in [0, 4, 9, usize::MAX] {
            let mut names = Vec::new();
            let mut decoded = Vec::new();
            let mut batches = 0;
            read_encoded_batches(text.as_bytes(), &Alphabet::Dna, 0, batch_size, |n, b| {
                batches += 1;
                decoded.extend((0..b.len()).map(|i| b.decode(i)));
                names.extend(n);
                Ok(())
            })
            .unwrap();
            assert_eq!(names, ["a", "b", "c x", "d"]);
            let want: Vec<String> = (0..whole.len()).map(|i| whole.decode(i)).collect();
            assert_eq!(decoded, want);
            assert_eq!(
                batches,
                match batch_size {
                    0 => 4,
                    4 => 3,
                    9 => 2,
                    _ => 1,
                }
            );
        }
        assert_eq!(decoded_first(text), "ACGTNACGT");
    }

    fn decoded_first(text: &str) -> String {
        let mut first = None;
        read_encoded_batches(text.as_bytes(), &Alphabet::Dna, 0, 0, |_, b| {
            first.get_or_insert_with(|| b.decode(0));
            Ok(())
        })
        .unwrap();
        first.unwrap()
    }
}
