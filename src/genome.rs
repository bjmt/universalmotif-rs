//! Random access to genome sequences on disk.
//!
//! [`Genome`] reads any part of a sequence without loading the rest, from a
//! FASTA file, plain or compressed with bgzip, or a 2bit file. FASTA files are
//! indexed as `samtools faidx` does, using the `.fai` (and for bgzip the
//! `.gzi`) file next to the FASTA file if there is one, and otherwise building
//! the index in memory with one pass through the file.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

use std::io;

use crate::error::{Error, Result, arg_err};
use crate::ranges::GRange;

/// The names, lengths and records of a FASTA index.
type FaiIndex = (Vec<String>, Vec<u64>, Vec<FaiRecord>);

/// Where a FASTA record's sequence starts, and how its lines are laid out, as
/// in a `.fai` index.
#[derive(Debug, Clone, PartialEq)]
struct FaiRecord {
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// A 2bit record: where its packed bases start, and its runs of N and of
/// soft-masked letters, as (start, length).
#[derive(Debug, Clone, PartialEq)]
struct TwoBitRecord {
    dna_offset: u64,
    n_blocks: Vec<(u64, u64)>,
    mask_blocks: Vec<(u64, u64)>,
}

#[derive(Debug, Clone)]
enum Source {
    Fasta(Vec<FaiRecord>),
    /// Records, and the start of each BGZF block as (compressed offset,
    /// uncompressed offset).
    Bgzf(Vec<FaiRecord>, Vec<(u64, u64)>),
    TwoBit(Vec<TwoBitRecord>),
}

/// A genome on disk, read a piece at a time.
///
/// # Example
///
/// ```no_run
/// use universalmotif::genome::Genome;
///
/// let genome = Genome::open("genome.fa.gz")?;
/// // The first 100 bases of chromosome 1, counting from 0.
/// let start = genome.fetch("1", 0, 100)?;
/// # Ok::<(), universalmotif::Error>(())
/// ```
#[derive(Debug, Clone)]
pub struct Genome {
    path: PathBuf,
    names: Vec<String>,
    lengths: Vec<u64>,
    index: HashMap<String, usize>,
    source: Source,
}

impl Genome {
    /// Open a FASTA file, plain or compressed with bgzip, or a 2bit file,
    /// telling them apart by their first bytes. Sequence names are the first
    /// word of each FASTA header.
    pub fn open(path: impl AsRef<Path>) -> Result<Genome> {
        let path = path.as_ref();
        Genome::open_io(path).map_err(|e| {
            Error::Argument(format!("could not read the genome {}: {e}", path.display()))
        })
    }

    fn open_io(path: &Path) -> io::Result<Genome> {
        let path = path.to_path_buf();
        let mut magic = [0u8; 4];
        let n = File::open(&path)?.read(&mut magic)?;
        let (names, lengths, source) = if n == 4 && is_twobit(&magic) {
            let (names, lengths, records) = read_twobit_index(&path)?;
            (names, lengths, Source::TwoBit(records))
        } else if n >= 2 && magic[..2] == [0x1f, 0x8b] {
            let blocks = bgzf_blocks(&path)?;
            let (names, lengths, records) = match read_fai(&fai_path(&path))? {
                Some(fai) => fai,
                None => build_fai(BufReader::new(MultiGzDecoder::new(File::open(&path)?)))?,
            };
            (names, lengths, Source::Bgzf(records, blocks))
        } else {
            let (names, lengths, records) = match read_fai(&fai_path(&path))? {
                Some(fai) => fai,
                None => build_fai(BufReader::new(File::open(&path)?))?,
            };
            (names, lengths, Source::Fasta(records))
        };
        let mut index = HashMap::with_capacity(names.len());
        for (i, n) in names.iter().enumerate() {
            if index.insert(n.clone(), i).is_some() {
                return Err(bad(format!(
                    "the sequence name '{n}' appears more than once"
                )));
            }
        }
        Ok(Genome {
            path,
            names,
            lengths,
            index,
            source,
        })
    }

    /// Sequence names, in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Sequence lengths, in file order.
    pub fn lengths(&self) -> &[u64] {
        &self.lengths
    }

    /// The number of sequences.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Whether there are no sequences.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// The position of a sequence in the file, if it is there.
    pub fn position(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    /// Positions `start` to `end` of the sequence called `name`, counted from 0
    /// with `end` excluded. Letters are returned as they are stored, so
    /// soft-masked letters stay in lower case.
    pub fn fetch(&self, name: &str, start: u64, end: u64) -> Result<String> {
        let Some(i) = self.position(name) else {
            return arg_err(format!("no sequence called '{name}' in the genome"));
        };
        self.fetch_at(i, start, end)
    }

    /// Many pieces at once, read in parallel on `nthreads` threads (0 for all
    /// available). Each range's `chrom` is the position of a sequence in the
    /// genome, and strands are ignored. The pieces are read in order of
    /// position, so pieces close together in a bgzip file share the work of
    /// decompressing it.
    pub fn fetch_many(&self, ranges: &[GRange], nthreads: usize) -> Result<Vec<String>> {
        use rayon::prelude::*;
        let mut tuples = Vec::with_capacity(ranges.len());
        for r in ranges {
            if r.start < 0 || r.end < r.start {
                return arg_err(format!(
                    "positions {} to {} are not a valid range",
                    r.start, r.end
                ));
            }
            let (s, e) = (r.start as u64, r.end as u64);
            self.check(r.chrom, s, e)?;
            tuples.push((r.chrom, s, e));
        }
        let ranges = &tuples;
        let mut order: Vec<usize> = (0..ranges.len()).collect();
        order.sort_by_key(|&k| (ranges[k].0, ranges[k].1));
        let pieces = crate::pvalue::with_threads(nthreads, || {
            let chunk = order
                .len()
                .div_ceil(8 * rayon::current_num_threads())
                .max(1);
            order
                .par_chunks(chunk)
                .map(|ks| {
                    let mut cursor = Cursor::default();
                    ks.iter()
                        .map(|&k| {
                            let (i, s, e) = ranges[k];
                            Ok((k, self.read_with(&mut cursor, i, s, e)?))
                        })
                        .collect::<Result<Vec<(usize, String)>>>()
                })
                .collect::<Result<Vec<_>>>()
        })??;
        let mut out = vec![String::new(); ranges.len()];
        for (k, seq) in pieces.into_iter().flatten() {
            out[k] = seq;
        }
        Ok(out)
    }

    /// [`Genome::fetch`] for the sequence at position `i`.
    pub fn fetch_at(&self, i: usize, start: u64, end: u64) -> Result<String> {
        self.check(i, start, end)?;
        self.read_with(&mut Cursor::default(), i, start, end)
    }

    fn check(&self, i: usize, start: u64, end: u64) -> Result<()> {
        let Some(&len) = self.lengths.get(i) else {
            return arg_err(format!("sequence {i} is out of range"));
        };
        if start > end || end > len {
            return arg_err(format!(
                "positions {start} to {end} are outside {}, which has {len} letters",
                self.names[i]
            ));
        }
        Ok(())
    }

    /// Read a checked piece, keeping the file open and the last stretch
    /// decompressed in `cursor` for the next piece.
    fn read_with(&self, cursor: &mut Cursor, i: usize, start: u64, end: u64) -> Result<String> {
        if start == end {
            return Ok(String::new());
        }
        self.read(cursor, i, start, end)
            .map_err(|e| Error::Argument(format!("could not read {}: {e}", self.path.display())))
    }

    fn read(&self, cursor: &mut Cursor, i: usize, start: u64, end: u64) -> io::Result<String> {
        if cursor.file.is_none() {
            cursor.file = Some(File::open(&self.path)?);
        }
        let f = cursor.file.as_mut().unwrap();
        match &self.source {
            Source::Fasta(records) => {
                let (a, b) = byte_range(&records[i], start, end);
                f.seek(SeekFrom::Start(a))?;
                let mut buf = vec![0u8; (b - a) as usize];
                f.read_exact(&mut buf)?;
                Ok(strip_newlines(buf))
            }
            Source::Bgzf(records, blocks) => {
                let (a, b) = byte_range(&records[i], start, end);
                let cached = matches!(&cursor.span,
                    Some((u, data)) if *u <= a && b <= u + data.len() as u64);
                if !cached {
                    // Decompress whole blocks, from the one holding the first
                    // byte to the one holding the last.
                    let k = blocks.partition_point(|&(_, u)| u <= a) - 1;
                    let next = blocks.partition_point(|&(_, u)| u < b);
                    let (coff, uoff) = blocks[k];
                    f.seek(SeekFrom::Start(coff))?;
                    let mut compressed = Vec::new();
                    match blocks.get(next) {
                        Some(&(cend, _)) => {
                            compressed.resize((cend - coff) as usize, 0);
                            f.read_exact(&mut compressed)?;
                        }
                        None => {
                            f.read_to_end(&mut compressed)?;
                        }
                    }
                    let mut data = Vec::new();
                    MultiGzDecoder::new(&compressed[..]).read_to_end(&mut data)?;
                    if (data.len() as u64) < b - uoff {
                        return Err(bad("the file ends before the sequence does"));
                    }
                    cursor.span = Some((uoff, data));
                }
                let (u, data) = cursor.span.as_ref().unwrap();
                Ok(strip_newlines(
                    data[(a - u) as usize..(b - u) as usize].to_vec(),
                ))
            }
            Source::TwoBit(records) => read_twobit(f, &records[i], start, end),
        }
    }
}

/// What a reader keeps from one piece to the next: the open file, and for
/// bgzip files the last stretch decompressed, with its uncompressed offset.
#[derive(Default)]
struct Cursor {
    file: Option<File>,
    span: Option<(u64, Vec<u8>)>,
}

fn fai_path(path: &Path) -> PathBuf {
    let mut p = path.as_os_str().to_owned();
    p.push(".fai");
    PathBuf::from(p)
}

/// The bytes holding positions `start` to `end` of a record, as a half-open
/// range of offsets in the (uncompressed) file.
fn byte_range(r: &FaiRecord, start: u64, end: u64) -> (u64, u64) {
    let at = |pos: u64| r.offset + (pos / r.line_bases) * r.line_width + pos % r.line_bases;
    (at(start), at(end - 1) + 1)
}

fn strip_newlines(mut buf: Vec<u8>) -> String {
    buf.retain(|&b| b != b'\n' && b != b'\r');
    String::from_utf8(buf).unwrap_or_else(|e| String::from_utf8_lossy(e.as_bytes()).into_owned())
}

/// Read a `.fai` index, if there is one.
fn read_fai(path: &Path) -> io::Result<Option<FaiIndex>> {
    let Ok(f) = File::open(path) else {
        return Ok(None);
    };
    let (mut names, mut lengths, mut records) = (Vec::new(), Vec::new(), Vec::new());
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let num = |k: usize| -> io::Result<u64> {
            fields
                .get(k)
                .and_then(|v| v.trim().parse().ok())
                .ok_or_else(|| bad(format!("bad line in {}: {line}", path.display())))
        };
        let (length, line_bases, line_width) = (num(1)?, num(3)?, num(4)?);
        if length > 0 && (line_bases == 0 || line_width < line_bases) {
            return Err(bad(format!("bad line in {}: {line}", path.display())));
        }
        names.push(fields[0].to_string());
        lengths.push(length);
        records.push(FaiRecord {
            offset: num(2)?,
            line_bases,
            line_width,
        });
    }
    Ok(Some((names, lengths, records)))
}

/// Index a FASTA file in one pass, as `samtools faidx` would. Every line of a
/// sequence but its last must have the same length.
fn build_fai<R: BufRead>(mut reader: R) -> io::Result<FaiIndex> {
    let (mut names, mut lengths, mut records) = (Vec::new(), Vec::new(), Vec::new());
    let mut offset = 0u64;
    let mut line = Vec::new();
    // For the record being read: its line layout so far and whether a short
    // line, which must be its last, has been seen.
    let mut short_seen = false;
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)? as u64;
        if n == 0 {
            break;
        }
        let content = line
            .iter()
            .rposition(|&b| b != b'\n' && b != b'\r')
            .map_or(0, |p| p + 1) as u64;
        if line.first() == Some(&b'>') {
            let header = String::from_utf8_lossy(&line[1..content as usize]);
            names.push(header.split_whitespace().next().unwrap_or("").to_string());
            lengths.push(0);
            records.push(FaiRecord {
                offset: offset + n,
                line_bases: 0,
                line_width: 0,
            });
            short_seen = false;
        } else if content > 0 {
            let Some(r) = records.last_mut() else {
                return Err(bad("sequence before the first FASTA header"));
            };
            let name = names.last().unwrap();
            if r.line_bases == 0 {
                r.line_bases = content;
                r.line_width = n;
            } else if short_seen
                || content > r.line_bases
                // Only the last line of the file can lack its line ending.
                || (n - content != r.line_width - r.line_bases && line.last() == Some(&b'\n'))
            {
                return Err(bad(format!(
                    "the lines of {name} differ in length, so it cannot be read in pieces; \
                     rewrap the file, for example with `samtools faidx`"
                )));
            } else if content < r.line_bases {
                short_seen = true;
            }
            *lengths.last_mut().unwrap() += content;
        }
        offset += n;
    }
    Ok((names, lengths, records))
}

/// The start of every BGZF block of a bgzip file, as (compressed offset,
/// uncompressed offset), from its `.gzi` index or from the block headers.
fn bgzf_blocks(path: &Path) -> io::Result<Vec<(u64, u64)>> {
    let mut gzi = path.as_os_str().to_owned();
    gzi.push(".gzi");
    if let Ok(mut f) = File::open(PathBuf::from(gzi)) {
        let mut buf = Vec::new();
        f.read_to_end(&mut buf)?;
        let u64_at = |k: usize| u64::from_le_bytes(buf[k..k + 8].try_into().unwrap());
        if buf.len() >= 8 {
            let n = u64_at(0) as usize;
            if buf.len() >= 8 + 16 * n {
                let mut out = vec![(0, 0)];
                out.extend((0..n).map(|i| (u64_at(8 + 16 * i), u64_at(16 + 16 * i))));
                return Ok(out);
            }
        }
        return Err(bad("the .gzi index is truncated"));
    }
    let mut f = File::open(path)?;
    let size = f.metadata()?.len();
    let (mut coff, mut uoff) = (0u64, 0u64);
    let mut out = Vec::new();
    while coff < size {
        let mut head = [0u8; 12];
        f.seek(SeekFrom::Start(coff))?;
        f.read_exact(&mut head)?;
        if head[..4] != [0x1f, 0x8b, 8, 4] {
            return Err(bad(NOT_BGZIP));
        }
        let xlen = u16::from_le_bytes([head[10], head[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        f.read_exact(&mut extra)?;
        let mut bsize = None;
        let mut k = 0;
        while k + 4 <= xlen {
            let slen = u16::from_le_bytes([extra[k + 2], extra[k + 3]]) as usize;
            if extra[k] == b'B' && extra[k + 1] == b'C' && slen == 2 && k + 6 <= xlen {
                bsize = Some(u16::from_le_bytes([extra[k + 4], extra[k + 5]]) as u64 + 1);
            }
            k += 4 + slen;
        }
        // A block holds at least its header, extra field and 8-byte footer.
        let Some(bsize) = bsize.filter(|&b| b >= 20 + xlen as u64) else {
            return Err(bad(NOT_BGZIP));
        };
        let mut isize = [0u8; 4];
        f.seek(SeekFrom::Start(coff + bsize - 4))?;
        f.read_exact(&mut isize)?;
        out.push((coff, uoff));
        coff += bsize;
        uoff += u32::from_le_bytes(isize) as u64;
    }
    if out.is_empty() {
        out.push((0, 0));
    }
    Ok(out)
}

const NOT_BGZIP: &str = "the file is gzipped but not with bgzip, so it cannot be read in \
                         pieces; recompress it with `bgzip`";

/// An error for a file that cannot be read as it should be.
fn bad(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

fn is_twobit(magic: &[u8; 4]) -> bool {
    *magic == [0x43, 0x27, 0x41, 0x1a] || *magic == [0x1a, 0x41, 0x27, 0x43]
}

/// Reads little- or big-endian integers from a 2bit file.
struct TwoBitReader {
    f: BufReader<File>,
    big: bool,
}

impl TwoBitReader {
    fn u32(&mut self) -> io::Result<u32> {
        let mut b = [0u8; 4];
        self.f.read_exact(&mut b)?;
        Ok(if self.big {
            u32::from_be_bytes(b)
        } else {
            u32::from_le_bytes(b)
        })
    }

    fn u64(&mut self) -> io::Result<u64> {
        let mut b = [0u8; 8];
        self.f.read_exact(&mut b)?;
        Ok(if self.big {
            u64::from_be_bytes(b)
        } else {
            u64::from_le_bytes(b)
        })
    }

    fn u32s(&mut self, n: usize) -> io::Result<Vec<u64>> {
        (0..n).map(|_| self.u32().map(u64::from)).collect()
    }
}

fn read_twobit_index(path: &Path) -> io::Result<(Vec<String>, Vec<u64>, Vec<TwoBitRecord>)> {
    let mut f = BufReader::new(File::open(path)?);
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic)?;
    let mut r = TwoBitReader {
        f,
        big: magic == [0x1a, 0x41, 0x27, 0x43],
    };
    let version = r.u32()?;
    if version > 1 {
        return Err(bad(format!("unknown 2bit version {version}")));
    }
    let count = r.u32()? as usize;
    r.u32()?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        let mut n = [0u8; 1];
        r.f.read_exact(&mut n)?;
        let mut name = vec![0u8; n[0] as usize];
        r.f.read_exact(&mut name)?;
        let offset = if version == 1 {
            r.u64()?
        } else {
            u64::from(r.u32()?)
        };
        entries.push((String::from_utf8_lossy(&name).into_owned(), offset));
    }
    let (mut names, mut lengths, mut records) = (Vec::new(), Vec::new(), Vec::new());
    for (name, offset) in entries {
        r.f.seek(SeekFrom::Start(offset))?;
        let dna_size = u64::from(r.u32()?);
        let n_count = r.u32()? as usize;
        let n_starts = r.u32s(n_count)?;
        let n_sizes = r.u32s(n_count)?;
        let m_count = r.u32()? as usize;
        let m_starts = r.u32s(m_count)?;
        let m_sizes = r.u32s(m_count)?;
        r.u32()?;
        let dna_offset = offset + 4 * (4 + 2 * n_count as u64 + 2 * m_count as u64);
        names.push(name);
        lengths.push(dna_size);
        let blocks = |starts: Vec<u64>, sizes: Vec<u64>| {
            let mut b: Vec<(u64, u64)> = starts.into_iter().zip(sizes).collect();
            b.sort_unstable();
            b
        };
        records.push(TwoBitRecord {
            dna_offset,
            n_blocks: blocks(n_starts, n_sizes),
            mask_blocks: blocks(m_starts, m_sizes),
        });
    }
    Ok((names, lengths, records))
}

fn read_twobit(f: &mut File, r: &TwoBitRecord, start: u64, end: u64) -> io::Result<String> {
    let (a, b) = (start / 4, (end - 1) / 4 + 1);
    f.seek(SeekFrom::Start(r.dna_offset + a))?;
    let mut packed = vec![0u8; (b - a) as usize];
    f.read_exact(&mut packed)?;
    const BASES: [u8; 4] = *b"TCAG";
    let mut out: Vec<u8> = (start..end)
        .map(|p| {
            let byte = packed[(p / 4 - a) as usize];
            BASES[((byte >> (6 - 2 * (p % 4))) & 3) as usize]
        })
        .collect();
    // Blocks are in order and do not overlap, so only those from the first
    // ending after `start` up to the first starting at `end` are needed.
    let mut apply = |blocks: &[(u64, u64)], f: &dyn Fn(u8) -> u8| {
        let first = blocks.partition_point(|&(s, n)| s + n <= start);
        for &(s, n) in blocks[first..].iter().take_while(|&&(s, _)| s < end) {
            let (lo, hi) = (s.max(start), (s + n).min(end));
            for p in lo..hi.max(lo) {
                let k = (p - start) as usize;
                out[k] = f(out[k]);
            }
        }
    };
    apply(&r.n_blocks, &|_| b'N');
    apply(&r.mask_blocks, &|c| c.to_ascii_lowercase());
    Ok(String::from_utf8(out).expect("bases are ASCII"))
}
