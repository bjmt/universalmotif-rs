//! Nucleotide sequences packed at two bits per position.
//!
//! The layout follows tangermeme's `PackedSequence`: four positions per
//! byte, the first in the two most significant bits, with letters coded in
//! alphabet order (A, C, G, T or U). Anything else (such as `N`) is stored
//! as a run of positions and packed as `A`.

/// A packed DNA or RNA sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedSeq {
    /// Packed letters, followed by two zero bytes so that [`PackedSeq::quad`]
    /// can read past the end.
    bytes: Vec<u8>,
    len: usize,
    /// Half-open runs of positions holding other letters.
    unknown: Vec<(usize, usize)>,
}

impl PackedSeq {
    /// Pack a sequence of the four letters `letters` (matched in either
    /// case).
    pub fn new(seq: &[u8], letters: [u8; 4]) -> PackedSeq {
        let mut code = [4u8; 256];
        for (i, &l) in letters.iter().enumerate() {
            code[l.to_ascii_uppercase() as usize] = i as u8;
            code[l.to_ascii_lowercase() as usize] = i as u8;
        }
        let mut bytes = vec![0u8; seq.len().div_ceil(4) + 2];
        let mut unknown: Vec<(usize, usize)> = Vec::new();
        for (i, &b) in seq.iter().enumerate() {
            let c = code[b as usize];
            if c == 4 {
                match unknown.last_mut() {
                    Some(run) if run.1 == i => run.1 = i + 1,
                    _ => unknown.push((i, i + 1)),
                }
            } else {
                bytes[i >> 2] |= c << (6 - 2 * (i & 3));
            }
        }
        PackedSeq {
            bytes,
            len: seq.len(),
            unknown,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Runs of positions that are not one of the four letters.
    pub fn unknown_runs(&self) -> &[(usize, usize)] {
        &self.unknown
    }

    /// The letter code (0 to 3) at `i`; unknown letters read as 0.
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
        let mut start = 0;
        self.unknown
            .iter()
            .copied()
            .chain(std::iter::once((self.len, self.len)))
            .filter_map(move |(a, b)| {
                let seg = (start, a);
                start = b;
                (seg.1 > seg.0).then_some(seg)
            })
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
        assert_eq!(p.unknown_runs(), &[(4, 6), (14, 15)]);
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
    }
}
