//! Sequences encoded for scanning.
//!
//! DNA and RNA are packed at two bits per position, four positions per byte,
//! the first in the two most significant bits, with letters coded 0 to 3 in the
//! order A, C, G, T or U. Other alphabets take a byte per position, holding the
//! letter's index in the alphabet. Anything outside the alphabet, such as `N`,
//! is coded as the first letter and recorded in runs of one repeated byte, so
//! decoding gives it back exactly; only the case of the alphabet's own letters
//! is lost when they are matched in either case.
//!
//! Sequences can be encoded whole or fed in pieces, such as a line of a FASTA
//! file at a time, with [`EncodedSequences::record`], so they never need to be
//! held as text.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::error::{Error, Result};
use crate::pvalue::with_threads;

/// A run of positions `start..end` holding the byte `letter`, which is not
/// in the alphabet.
pub(crate) type UnknownRun = (usize, usize, u8);

/// Record `b` at position `i` as unknown, extending the last run if it
/// holds the same byte and ends at `i`.
#[inline]
fn add_unknown(unknown: &mut Vec<UnknownRun>, i: usize, b: u8) {
    match unknown.last_mut() {
        Some(run) if run.1 == i && run.2 == b => run.1 = i + 1,
        _ => unknown.push((i, i + 1, b)),
    }
}

/// Overwrite the unknown positions of `out`, which holds positions
/// `start..end`, with their original bytes.
fn patch_unknown(unknown: &[UnknownRun], start: usize, end: usize, out: &mut [u8]) {
    let k = unknown.partition_point(|r| r.1 <= start);
    for &(a, b, letter) in &unknown[k..] {
        if a >= end {
            break;
        }
        out[a.max(start) - start..b.min(end) - start].fill(letter);
    }
}

/// The stretches between runs of unknown letters.
fn segments(unknown: &[UnknownRun], len: usize) -> impl Iterator<Item = (usize, usize)> + '_ {
    let mut start = 0;
    unknown
        .iter()
        .map(|&(a, b, _)| (a, b))
        .chain(std::iter::once((len, len)))
        .filter_map(move |(a, b)| {
            let seg = (start, a);
            start = b;
            (seg.1 > seg.0).then_some(seg)
        })
}

/// Letter codes for every byte: the letter's index, or `u8::MAX` for bytes
/// outside the alphabet.
fn code_table(letters: &[u8], case_insensitive: bool) -> [u8; 256] {
    let mut code = [u8::MAX; 256];
    for (i, &l) in letters.iter().enumerate() {
        code[l as usize] = i as u8;
        if case_insensitive {
            code[l.to_ascii_uppercase() as usize] = i as u8;
            code[l.to_ascii_lowercase() as usize] = i as u8;
        }
    }
    code
}

/// A packed DNA or RNA sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PackedSeq {
    /// Packed letters, followed by two zero bytes so that [`PackedSeq::quad`]
    /// can read past the end.
    bytes: Vec<u8>,
    len: usize,
    unknown: Vec<UnknownRun>,
}

impl PackedSeq {
    #[cfg(test)]
    /// Pack a sequence of the four letters `letters`, matched in either case.
    pub fn new(seq: &[u8], letters: [u8; 4]) -> PackedSeq {
        let mut b = PackedSeqBuilder::with_capacity(letters, seq.len());
        b.push(seq);
        b.finish()
    }

    /// The number of letters.
    pub fn len(&self) -> usize {
        self.len
    }

    #[cfg(test)]
    /// Runs of positions that are not one of the four letters.
    pub fn unknown_runs(&self) -> &[UnknownRun] {
        &self.unknown
    }

    /// The letter code at `i`, from 0 to 3; unknown letters read as 0.
    pub fn get(&self, i: usize) -> u8 {
        (self.bytes[i >> 2] >> (6 - 2 * (i & 3))) & 3
    }

    /// The four letters starting at `i` as one byte, first letter in the
    /// high bits. Letters past the end read as 0.
    #[inline(always)]
    pub fn quad(&self, i: usize) -> u8 {
        let k = i >> 2;
        let w = ((self.bytes[k] as u16) << 8) | self.bytes[k + 1] as u16;
        (w >> (8 - 2 * (i & 3))) as u8
    }

    /// Half-open stretches of the sequence without unknown letters.
    pub fn known_segments(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        segments(&self.unknown, self.len)
    }

    /// Append positions `start..end` to `out`, writing the alphabet's
    /// letters as `letters` and unknown letters as they were.
    pub fn decode_into(&self, start: usize, end: usize, letters: &[u8; 4], out: &mut Vec<u8>) {
        let base = out.len();
        out.extend((start..end).map(|i| letters[self.get(i) as usize]));
        patch_unknown(&self.unknown, start, end, &mut out[base..]);
    }

    /// Memory held, in bytes.
    pub fn heap_bytes(&self) -> usize {
        self.bytes.capacity() + self.unknown.capacity() * std::mem::size_of::<UnknownRun>()
    }
}

/// Builds a [`PackedSeq`] from pieces of a sequence.
#[derive(Debug, Clone)]
pub(crate) struct PackedSeqBuilder {
    code: [u8; 256],
    bytes: Vec<u8>,
    len: usize,
    unknown: Vec<UnknownRun>,
}

impl PackedSeqBuilder {
    #[cfg(test)]
    /// A builder for sequences of the four letters `letters`, matched in either
    /// case.
    pub fn new(letters: [u8; 4]) -> PackedSeqBuilder {
        PackedSeqBuilder::with_capacity(letters, 0)
    }

    #[cfg(test)]
    /// As [`PackedSeqBuilder::new`], with room for `len` letters.
    pub fn with_capacity(letters: [u8; 4], len: usize) -> PackedSeqBuilder {
        PackedSeqBuilder {
            code: code_table(&letters, true),
            bytes: Vec::with_capacity(len.div_ceil(4) + 2),
            len: 0,
            unknown: Vec::new(),
        }
    }

    /// The number of letters added so far.
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    fn push_one(&mut self, b: u8) {
        let i = self.len;
        if i & 3 == 0 {
            self.bytes.push(0);
        }
        let c = self.code[b as usize];
        if c == u8::MAX {
            add_unknown(&mut self.unknown, i, b);
        } else {
            self.bytes[i >> 2] |= c << (6 - 2 * (i & 3));
        }
        self.len += 1;
    }

    /// Add the letters of `piece` to the end of the sequence.
    pub fn push(&mut self, piece: &[u8]) {
        let mut rest = piece;
        while self.len & 3 != 0 && !rest.is_empty() {
            self.push_one(rest[0]);
            rest = &rest[1..];
        }
        let (quads, remainder) = rest.as_chunks::<4>();
        for q in quads {
            let c = [
                self.code[q[0] as usize],
                self.code[q[1] as usize],
                self.code[q[2] as usize],
                self.code[q[3] as usize],
            ];
            if c.iter().all(|&x| x < 4) {
                self.bytes.push(c[0] << 6 | c[1] << 4 | c[2] << 2 | c[3]);
                self.len += 4;
            } else {
                q.iter().for_each(|&b| self.push_one(b));
            }
        }
        remainder.iter().for_each(|&b| self.push_one(b));
    }

    /// Finish the sequence, ready to be scanned.
    pub fn finish(mut self) -> PackedSeq {
        self.bytes.extend([0, 0]);
        PackedSeq {
            bytes: self.bytes,
            len: self.len,
            unknown: self.unknown,
        }
    }
}

/// A sequence over any alphabet of at most 255 letters, one byte per
/// letter: its index in the alphabet. Letters outside it are coded as 0 and
/// recorded as runs, as for [`PackedSeq`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CodedSeq {
    /// Letter codes, followed by one zero so pairs can be read to the end.
    codes: Vec<u8>,
    unknown: Vec<UnknownRun>,
}

impl CodedSeq {
    #[cfg(test)]
    /// Code a sequence over `letters`, of which there can be at most 255. With
    /// `case_insensitive`, letters match in either case.
    pub fn new(seq: &[u8], letters: &[u8], case_insensitive: bool) -> CodedSeq {
        let mut b = CodedSeqBuilder::with_capacity(letters, case_insensitive, seq.len());
        b.push(seq);
        b.finish()
    }

    /// The number of letters.
    pub fn len(&self) -> usize {
        self.codes.len() - 1
    }

    /// The letter code at `i`. Unknown letters read as 0, as does the position
    /// just past the end.
    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        self.codes[i]
    }

    #[cfg(test)]
    /// Runs of positions that are not letters of the alphabet.
    pub fn unknown_runs(&self) -> &[UnknownRun] {
        &self.unknown
    }

    /// Half-open stretches of the sequence without unknown letters.
    pub fn known_segments(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        segments(&self.unknown, self.len())
    }

    /// Append positions `start..end` to `out`, writing the alphabet's
    /// letters as `letters` and unknown letters as they were.
    pub fn decode_into(&self, start: usize, end: usize, letters: &[u8], out: &mut Vec<u8>) {
        let base = out.len();
        out.extend(self.codes[start..end].iter().map(|&c| letters[c as usize]));
        patch_unknown(&self.unknown, start, end, &mut out[base..]);
    }

    /// Memory held, in bytes.
    pub fn heap_bytes(&self) -> usize {
        self.codes.capacity() + self.unknown.capacity() * std::mem::size_of::<UnknownRun>()
    }
}

/// Builds a [`CodedSeq`] from pieces of a sequence.
#[derive(Debug, Clone)]
pub(crate) struct CodedSeqBuilder {
    code: [u8; 256],
    codes: Vec<u8>,
    unknown: Vec<UnknownRun>,
}

impl CodedSeqBuilder {
    #[cfg(test)]
    /// A builder for sequences over `letters`, of which there can be at most
    /// 255. With `case_insensitive`, letters match in either case.
    pub fn new(letters: &[u8], case_insensitive: bool) -> CodedSeqBuilder {
        CodedSeqBuilder::with_capacity(letters, case_insensitive, 0)
    }

    #[cfg(test)]
    /// As [`CodedSeqBuilder::new`], with room for `len` letters.
    pub fn with_capacity(letters: &[u8], case_insensitive: bool, len: usize) -> CodedSeqBuilder {
        CodedSeqBuilder {
            code: code_table(letters, case_insensitive),
            codes: Vec::with_capacity(len + 1),
            unknown: Vec::new(),
        }
    }

    /// The number of letters added so far.
    pub fn len(&self) -> usize {
        self.codes.len()
    }

    /// Add the letters of `piece` to the end of the sequence.
    pub fn push(&mut self, piece: &[u8]) {
        for &b in piece {
            let c = self.code[b as usize];
            if c == u8::MAX {
                add_unknown(&mut self.unknown, self.codes.len(), b);
                self.codes.push(0);
            } else {
                self.codes.push(c);
            }
        }
    }

    /// Finish the sequence, ready to be scanned.
    pub fn finish(mut self) -> CodedSeq {
        self.codes.push(0);
        CodedSeq {
            codes: self.codes,
            unknown: self.unknown,
        }
    }
}

/// Encoded sequences, one way or the other.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Seqs {
    Packed(Vec<PackedSeq>),
    Coded(Vec<CodedSeq>),
}

/// A set of sequences encoded for one alphabet: packed at two bits per
/// position for DNA and RNA, a byte per position otherwise. A human genome
/// packs into about 0.8 GB, a quarter of its size as text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodedSequences {
    alphabet: Alphabet,
    letters: Vec<u8>,
    seqs: Seqs,
}

/// One sequence being encoded a piece at a time, started by
/// [`EncodedSequences::record`].
#[derive(Debug, Clone)]
pub struct RecordEncoder(RecordKind);

#[derive(Debug, Clone)]
enum RecordKind {
    Packed(PackedSeqBuilder),
    Coded(CodedSeqBuilder),
}

impl RecordEncoder {
    /// Make room for `additional` more letters.
    pub fn reserve(&mut self, additional: usize) {
        match &mut self.0 {
            RecordKind::Packed(b) => b
                .bytes
                .reserve((b.len + additional).div_ceil(4) + 2 - b.bytes.len()),
            RecordKind::Coded(b) => b.codes.reserve(additional + 1),
        }
    }

    /// Add the letters of `piece` to the end of the sequence.
    pub fn push(&mut self, piece: &[u8]) {
        match &mut self.0 {
            RecordKind::Packed(b) => b.push(piece),
            RecordKind::Coded(b) => b.push(piece),
        }
    }

    /// The number of letters added so far.
    pub fn len(&self) -> usize {
        match &self.0 {
            RecordKind::Packed(b) => b.len(),
            RecordKind::Coded(b) => b.len(),
        }
    }

    /// Whether no letters have been added.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl EncodedSequences {
    /// No sequences yet, for `alphabet`. Letters of DNA, RNA and amino acid
    /// sequences will match in either case, and letters of custom alphabets
    /// exactly. Alphabets of more than 255 letters, or with letters outside
    /// ASCII, cannot be encoded.
    pub fn empty(alphabet: &Alphabet) -> Result<EncodedSequences> {
        if alphabet.len() > 255 {
            return Err(Error::Unsupported(
                "alphabets of more than 255 letters cannot be scanned".into(),
            ));
        }
        if !alphabet.letters().iter().all(char::is_ascii) {
            return Err(Error::Unsupported(
                "alphabets with letters outside ASCII cannot be scanned".into(),
            ));
        }
        let letters: Vec<u8> = match alphabet {
            Alphabet::Dna => b"ACGT".to_vec(),
            Alphabet::Rna => b"ACGU".to_vec(),
            _ => alphabet.letters().iter().map(|&c| c as u8).collect(),
        };
        let seqs = if alphabet.is_nucleotide() {
            Seqs::Packed(Vec::new())
        } else {
            Seqs::Coded(Vec::new())
        };
        Ok(EncodedSequences {
            alphabet: alphabet.clone(),
            letters,
            seqs,
        })
    }

    /// Encode `sequences` for `alphabet` with `nthreads` threads, or all
    /// available for 0. [`EncodedSequences::empty`] describes how letters are
    /// matched.
    pub fn new<S: AsRef<str> + Sync>(
        sequences: &[S],
        alphabet: &Alphabet,
        nthreads: usize,
    ) -> Result<EncodedSequences> {
        let mut out = EncodedSequences::empty(alphabet)?;
        let encode = |s: &S| {
            let mut r = out.record();
            r.reserve(s.as_ref().len());
            r.push(s.as_ref().as_bytes());
            r
        };
        let records: Vec<RecordEncoder> =
            with_threads(nthreads, || sequences.par_iter().map(encode).collect())?;
        for r in records {
            out.push(r);
        }
        Ok(out)
    }

    /// A new, empty sequence to feed pieces to and then add with
    /// [`EncodedSequences::push`].
    pub fn record(&self) -> RecordEncoder {
        self.record_reusing(Vec::new())
    }

    /// Like [`EncodedSequences::record`], but storing the letters in `buf`
    /// after clearing it, so its memory is used again.
    pub(crate) fn record_reusing(&self, mut buf: Vec<u8>) -> RecordEncoder {
        buf.clear();
        RecordEncoder(match &self.seqs {
            Seqs::Packed(_) => RecordKind::Packed(PackedSeqBuilder {
                code: code_table(&self.letters[..4], true),
                bytes: buf,
                len: 0,
                unknown: Vec::new(),
            }),
            Seqs::Coded(_) => RecordKind::Coded(CodedSeqBuilder {
                code: code_table(&self.letters, !matches!(self.alphabet, Alphabet::Custom(_))),
                codes: buf,
                unknown: Vec::new(),
            }),
        })
    }

    /// Release memory reserved beyond what the sequences need.
    pub fn shrink_to_fit(&mut self) {
        match &mut self.seqs {
            Seqs::Packed(v) => v.iter_mut().for_each(|s| s.bytes.shrink_to_fit()),
            Seqs::Coded(v) => v.iter_mut().for_each(|s| s.codes.shrink_to_fit()),
        }
    }

    /// The buffers holding the letters, to be used again with
    /// [`EncodedSequences::record_reusing`].
    pub(crate) fn into_buffers(self) -> Vec<Vec<u8>> {
        match self.seqs {
            Seqs::Packed(v) => v.into_iter().map(|s| s.bytes).collect(),
            Seqs::Coded(v) => v.into_iter().map(|s| s.codes).collect(),
        }
    }

    /// Add a sequence made with [`EncodedSequences::record`] on this set or on
    /// one for the same alphabet.
    pub fn push(&mut self, record: RecordEncoder) {
        match (&mut self.seqs, record.0) {
            (Seqs::Packed(v), RecordKind::Packed(b)) => v.push(b.finish()),
            (Seqs::Coded(v), RecordKind::Coded(b)) => v.push(b.finish()),
            _ => panic!("a record encoded for another alphabet"),
        }
    }

    /// The alphabet the sequences were encoded for.
    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    /// The number of sequences.
    pub fn len(&self) -> usize {
        match &self.seqs {
            Seqs::Packed(v) => v.len(),
            Seqs::Coded(v) => v.len(),
        }
    }

    /// Whether there are no sequences.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The length of sequence `i`.
    pub fn seq_len(&self, i: usize) -> usize {
        match &self.seqs {
            Seqs::Packed(v) => v[i].len(),
            Seqs::Coded(v) => v[i].len(),
        }
    }

    /// The total length of the sequences.
    pub fn total_len(&self) -> usize {
        (0..self.len()).map(|i| self.seq_len(i)).sum()
    }

    /// Append the 0-based, half-open positions `start..end` of sequence `i` to
    /// `out`. Letters of the alphabet are written in upper case unless it is
    /// custom, and other letters as they were.
    pub fn decode_into(&self, i: usize, start: usize, end: usize, out: &mut Vec<u8>) {
        match &self.seqs {
            Seqs::Packed(v) => {
                v[i].decode_into(start, end, self.letters[..4].try_into().unwrap(), out)
            }
            Seqs::Coded(v) => v[i].decode_into(start, end, &self.letters, out),
        }
    }

    /// Sequence `i` as text, decoded as by [`EncodedSequences::decode_into`].
    pub fn decode(&self, i: usize) -> String {
        let mut out = Vec::with_capacity(self.seq_len(i));
        self.decode_into(i, 0, self.seq_len(i), &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    /// Memory held, in bytes.
    pub fn heap_bytes(&self) -> usize {
        match &self.seqs {
            Seqs::Packed(v) => v.iter().map(PackedSeq::heap_bytes).sum(),
            Seqs::Coded(v) => v.iter().map(CodedSeq::heap_bytes).sum(),
        }
    }

    pub(crate) fn seqs(&self) -> &Seqs {
        &self.seqs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_and_reads_back() {
        let seq = b"ACGTNNacgtTGCAn";
        let p = PackedSeq::new(seq, *b"ACGT");
        assert_eq!(p.len(), 15);
        assert_eq!(p.unknown_runs(), &[(4, 6, b'N'), (14, 15, b'n')]);
        let codes: Vec<u8> = (0..15).map(|i| p.get(i)).collect();
        assert_eq!(codes, [0, 1, 2, 3, 0, 0, 0, 1, 2, 3, 3, 2, 1, 0, 0]);
        for i in 0..15 {
            let want = (0..4).fold(0u8, |q, j| {
                (q << 2) | if i + j < 15 { codes[i + j] } else { 0 }
            });
            assert_eq!(p.quad(i), want, "quad at {i}");
        }
        let segs: Vec<_> = p.known_segments().collect();
        assert_eq!(segs, [(0, 4), (6, 14)]);
        let rna = PackedSeq::new(b"ACGU", *b"ACGU");
        assert_eq!(rna.quad(0), 0b00011011);
        assert!(
            PackedSeq::new(b"", *b"ACGT")
                .known_segments()
                .next()
                .is_none()
        );
        let mut out = Vec::new();
        p.decode_into(2, 15, b"ACGT", &mut out);
        assert_eq!(out, b"GTNNACGTTGCAn");
    }

    #[test]
    fn pieces_pack_like_the_whole() {
        // Every split point and piece size, so the four-letter fast path
        // meets every alignment and unknown letters fall across pieces.
        let seq: &[u8] = b"ACGTRYacgtNNNNnnACGTTGCAxACGTACGTAC-GT";
        let whole = PackedSeq::new(seq, *b"ACGT");
        for size in 1..12 {
            let mut b = PackedSeqBuilder::new(*b"ACGT");
            for piece in seq.chunks(size) {
                b.push(piece);
            }
            assert_eq!(b.finish(), whole, "pieces of {size}");
        }
        let mut out = Vec::new();
        whole.decode_into(0, seq.len(), b"ACGT", &mut out);
        // Upper case for the alphabet's letters, the rest as they were.
        assert_eq!(out, b"ACGTRYACGTNNNNnnACGTTGCAxACGTACGTAC-GT");
        let coded = CodedSeq::new(seq, b"ACGT", true);
        let mut cb = CodedSeqBuilder::new(b"ACGT", true);
        for piece in seq.chunks(5) {
            cb.push(piece);
        }
        assert_eq!(cb.finish(), coded);
    }

    #[test]
    fn codes_any_alphabet() {
        let aa = CodedSeq::new(b"MKvXXlA*", b"ACDEFGHIKLMNPQRSTVWY", true);
        assert_eq!(aa.len(), 8);
        assert_eq!(aa.unknown_runs(), &[(3, 5, b'X'), (7, 8, b'*')]);
        assert_eq!(
            (aa.get(0), aa.get(1), aa.get(2), aa.get(5), aa.get(6)),
            (10, 8, 17, 9, 0)
        );
        assert_eq!(aa.known_segments().collect::<Vec<_>>(), [(0, 3), (5, 7)]);
        let mut out = Vec::new();
        aa.decode_into(1, 8, b"ACDEFGHIKLMNPQRSTVWY", &mut out);
        assert_eq!(out, b"KVXXLA*");
        // Custom letters match exactly.
        let custom = CodedSeq::new(b"xyXy", b"xy", false);
        assert_eq!(custom.unknown_runs(), &[(2, 3, b'X')]);
    }

    #[test]
    fn encoded_sets() {
        let seqs = ["acgtNNACGT", "", "TTTT"];
        let e = EncodedSequences::new(&seqs, &Alphabet::Dna, 2).unwrap();
        assert_eq!((e.len(), e.seq_len(0), e.total_len()), (3, 10, 14));
        assert_eq!(e.decode(0), "ACGTNNACGT");
        assert_eq!(e.decode(1), "");
        let mut by_piece = EncodedSequences::empty(&Alphabet::Dna).unwrap();
        for s in seqs {
            let mut r = by_piece.record();
            for piece in s.as_bytes().chunks(3) {
                r.push(piece);
            }
            by_piece.push(r);
        }
        assert_eq!(by_piece, e);
        let rna = EncodedSequences::new(&["ACGUT"], &Alphabet::Rna, 1).unwrap();
        assert_eq!(rna.decode(0), "ACGUT");
        let aa = EncodedSequences::new(&["mkvXb"], &Alphabet::Aa, 1).unwrap();
        assert_eq!(aa.decode(0), "MKVXb");
        let custom = Alphabet::custom("aBc".chars()).unwrap();
        assert_eq!(
            EncodedSequences::new(&["abcABC"], &custom, 1)
                .unwrap()
                .decode(0),
            "abcABC"
        );
        let wide = Alphabet::custom("aé".chars()).unwrap();
        assert!(EncodedSequences::empty(&wide).is_err());
    }
}
