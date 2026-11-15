//! Masking sequences by pattern with [`mask_seqs`], by position with
//! [`mask_ranges`], and in low-complexity regions with [`mask_complexity`].

use std::collections::VecDeque;

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::error::{Result, arg_err};
use crate::pvalue::with_threads;
use crate::scan::reverse_complement_seq;

fn check_letter(letter: char) -> Result<u8> {
    if letter.is_ascii() {
        Ok(letter as u8)
    } else {
        arg_err("the masking letter must be a single ASCII character")
    }
}

/// Mask every occurrence of `pattern` in `seq`, overlapping ones included,
/// matching letters in either case.
fn mask_pattern(seq: &mut [u8], pattern: &[u8], letter: u8) {
    let m = pattern.len();
    if m == 0 || seq.len() < m {
        return;
    }
    let mut masked = vec![false; seq.len()];
    for i in 0..=seq.len() - m {
        if seq[i..i + m].eq_ignore_ascii_case(pattern) {
            masked[i..i + m].fill(true);
        }
    }
    for (c, m) in seq.iter_mut().zip(masked) {
        if m {
            *c = letter;
        }
    }
}

/// Replace every occurrence of `pattern` in the sequences with `letter`.
/// Overlapping occurrences are all masked, and letters match in either case.
/// For DNA and RNA, `rc` masks the occurrences of the reverse complement as
/// well.
pub fn mask_seqs<S: AsRef<str>>(
    seqs: &[S],
    pattern: &str,
    rc: bool,
    letter: char,
    alphabet: &Alphabet,
) -> Result<Vec<String>> {
    if pattern.is_empty() {
        return arg_err("the pattern must not be empty");
    }
    let letter = check_letter(letter)?;
    if rc && !alphabet.is_nucleotide() {
        return arg_err("rc is only valid for DNA and RNA sequences");
    }
    let rc_pattern = rc.then(|| reverse_complement_seq(pattern, *alphabet == Alphabet::Rna));
    Ok(seqs
        .iter()
        .map(|s| {
            let mut bytes = s.as_ref().as_bytes().to_vec();
            mask_pattern(&mut bytes, pattern.as_bytes(), letter);
            if let Some(p) = &rc_pattern {
                mask_pattern(&mut bytes, p.as_bytes(), letter);
            }
            String::from_utf8_lossy(&bytes).into_owned()
        })
        .collect())
}

/// Replace the positions in `ranges` of `seq` with `letter`, where each range
/// covers `start..end`, counting from 0. Ranges may overlap; one reaching
/// outside the sequence is an error.
pub fn mask_ranges(seq: &str, ranges: &[(i64, i64)], letter: char) -> Result<String> {
    let letter = check_letter(letter)?;
    let mut bytes = seq.as_bytes().to_vec();
    let n = bytes.len() as i64;
    for &(start, end) in ranges {
        if end <= start {
            continue;
        }
        if start < 0 || end > n {
            return arg_err(format!(
                "the range {start}-{end} is outside the sequence (length {n})"
            ));
        }
        bytes[start as usize..end as usize].fill(letter);
    }
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// Options for [`mask_complexity`] and [`dust_intervals`].
#[derive(Debug, Clone, PartialEq)]
pub struct DustOptions {
    /// The window DUST scores, in bases; 64 by default, as in dustmasker.
    pub window: usize,
    /// The score above which a region is masked; 20 by default, as in
    /// dustmasker. Lower values mask more.
    pub threshold: u32,
    /// Also mask lowercase letters, as in soft-masked genomes.
    pub lowercase: bool,
    /// The masking letter; `-` by default.
    pub letter: char,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for DustOptions {
    fn default() -> Self {
        DustOptions {
            window: 64,
            threshold: 20,
            lowercase: false,
            letter: '-',
            nthreads: 1,
        }
    }
}

// Symmetric DUST, ported from sdust.c of minimap2
// (https://github.com/lh3/minimap2), which is under the MIT License:
// Copyright (c) 2018- Dana-Farber Cancer Institute, 2017-2018 Broad
// Institute, Inc.
const DUST_WORD: usize = 3;

/// An interval whose score is above the threshold and no lower than that of
/// any interval inside it, with its score `r` over `l` word pairs.
struct Perfect {
    start: usize,
    finish: usize,
    r: i64,
    l: i64,
}

/// The state of symmetric DUST over one run of A, C, G and T.
struct Dust {
    window: usize,
    threshold: i64,
    words: VecDeque<usize>,
    // Words of the window, and of its suffix of L words in which no word
    // appears too often, with their counts and scores.
    big_l: usize,
    rw: i64,
    rv: i64,
    cw: [i64; 64],
    cv: [i64; 64],
    perfect: VecDeque<Perfect>,
}

impl Dust {
    fn new(window: usize, threshold: u32) -> Self {
        Dust {
            window,
            threshold: i64::from(threshold),
            words: VecDeque::new(),
            big_l: 0,
            rw: 0,
            rv: 0,
            cw: [0; 64],
            cv: [0; 64],
            perfect: VecDeque::new(),
        }
    }

    fn shift_window(&mut self, t: usize) {
        if self.words.len() > self.window - DUST_WORD {
            let s = self.words.pop_front().unwrap();
            self.cw[s] -= 1;
            self.rw -= self.cw[s];
            if self.big_l > self.words.len() {
                self.big_l -= 1;
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
            }
        }
        self.words.push_back(t);
        self.big_l += 1;
        self.rw += self.cw[t];
        self.cw[t] += 1;
        self.rv += self.cv[t];
        self.cv[t] += 1;
        if self.cv[t] * 10 > self.threshold * 2 {
            loop {
                let s = self.words[self.words.len() - self.big_l];
                self.cv[s] -= 1;
                self.rv -= self.cv[s];
                self.big_l -= 1;
                if s == t {
                    break;
                }
            }
        }
    }

    /// Add the perfect intervals ending at the newest word, keeping them
    /// sorted by descending start.
    fn find_perfect(&mut self, start: usize) {
        let mut c = self.cv;
        let mut r = self.rv;
        let (mut max_r, mut max_l) = (0i64, 0i64);
        let n = self.words.len();
        // sdust rescans the intervals from the first for every word, but the
        // best score so far is a running maximum, so the scan can carry on
        // from where it stopped.
        let mut j = 0;
        for i in (0..n - self.big_l).rev() {
            let t = self.words[i];
            r += c[t];
            c[t] += 1;
            let (new_r, new_l) = (r, (n - i - 1) as i64);
            if new_r * 10 > self.threshold * new_l {
                while j < self.perfect.len() && self.perfect[j].start >= i + start {
                    let p = &self.perfect[j];
                    if max_r == 0 || p.r * max_l > max_r * p.l {
                        (max_r, max_l) = (p.r, p.l);
                    }
                    j += 1;
                }
                if max_r == 0 || new_r * max_l >= max_r * new_l {
                    (max_r, max_l) = (new_r, new_l);
                    self.perfect.insert(
                        j,
                        Perfect {
                            start: i + start,
                            finish: n + DUST_WORD - 1 + start,
                            r: new_r,
                            l: new_l,
                        },
                    );
                    j += 1;
                }
            }
        }
    }

    /// Move the perfect interval starting first into `out` once the window
    /// has passed its start, and drop those it has left behind.
    fn save(&mut self, start: usize, out: &mut Vec<(usize, usize)>) {
        let Some(p) = self.perfect.back() else { return };
        if p.start >= start {
            return;
        }
        match out.last_mut() {
            Some(last) if p.start <= last.1 => last.1 = last.1.max(p.finish),
            _ => out.push((p.start, p.finish)),
        }
        let keep = self
            .perfect
            .iter()
            .rposition(|q| q.start >= start)
            .map_or(0, |i| i + 1);
        self.perfect.truncate(keep);
    }

    /// Run over the bases of `run`, coded 0 to 3, which start at `offset`.
    fn run(mut self, run: &[u8], offset: usize, out: &mut Vec<(usize, usize)>) {
        let mut t = 0usize;
        for (k, &b) in run.iter().enumerate() {
            let l = k + 1;
            t = ((t << 2) | b as usize) & 63;
            if l >= DUST_WORD {
                let start = l.saturating_sub(self.window) + offset;
                self.save(start, out);
                self.shift_window(t);
                if self.rw * 10 > (self.big_l as i64) * self.threshold {
                    self.find_perfect(start);
                }
            }
        }
        // As in sdust, the flush starts one past where the window would.
        let mut start = (run.len() + 1).saturating_sub(self.window) + offset + 1;
        while !self.perfect.is_empty() {
            self.save(start, out);
            start += 1;
        }
    }
}

fn dust_code(c: u8) -> Option<u8> {
    match c {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

fn check_dust(opts: &DustOptions) -> Result<()> {
    if opts.window <= DUST_WORD {
        return arg_err(format!(
            "the DUST window must be longer than {DUST_WORD} bases"
        ));
    }
    if opts.threshold == 0 {
        return arg_err("the DUST threshold must be above 0");
    }
    Ok(())
}

fn dust_unchecked(seq: &[u8], opts: &DustOptions) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut codes = Vec::with_capacity(seq.len());
    let mut i = 0;
    while i < seq.len() {
        // Letters other than A, C, G, T and U split the sequence into runs
        // that are scored on their own.
        if dust_code(seq[i]).is_none() {
            i += 1;
            continue;
        }
        let begin = i;
        codes.clear();
        while i < seq.len() {
            match dust_code(seq[i]) {
                Some(b) => codes.push(b),
                None => break,
            }
            i += 1;
        }
        Dust::new(opts.window, opts.threshold).run(&codes, begin, &mut out);
    }
    if opts.lowercase {
        let mut i = 0;
        while i < seq.len() {
            if seq[i].is_ascii_lowercase() {
                let begin = i;
                while i < seq.len() && seq[i].is_ascii_lowercase() {
                    i += 1;
                }
                out.push((begin, i));
            } else {
                i += 1;
            }
        }
        out.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::with_capacity(out.len());
        for (s, e) in out {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        out = merged;
    }
    out
}

/// The low-complexity regions of a DNA or RNA sequence, found with the
/// symmetric DUST algorithm of Morgulis et al. (2006) as implemented in
/// sdust by Heng Li, which gives the same regions as dustmasker. Each region
/// covers `start..end`, counting from 0, and the regions are sorted and do
/// not overlap. Letters other than A, C, G, T and U, in either case, split
/// the sequence into parts that are scored separately. With
/// `opts.lowercase`, lowercase letters are included as well.
///
/// # References
///
/// - Morgulis A, Gertz EM, Schäffer AA, Agarwala R (2006). A fast and symmetric
///   DUST implementation to mask low-complexity DNA sequences. *Journal of
///   Computational Biology* 13(5):1028–1040.
///   <https://doi.org/10.1089/cmb.2006.13.1028>
pub fn dust_intervals(seq: &str, opts: &DustOptions) -> Result<Vec<(usize, usize)>> {
    check_dust(opts)?;
    Ok(dust_unchecked(seq.as_bytes(), opts))
}

/// Replace the low-complexity regions of DNA or RNA sequences, found by
/// [`dust_intervals`], with `opts.letter`.
pub fn mask_complexity<S: AsRef<str> + Sync>(
    seqs: &[S],
    opts: &DustOptions,
) -> Result<Vec<String>> {
    check_dust(opts)?;
    let letter = check_letter(opts.letter)?;
    with_threads(opts.nthreads, || {
        seqs.par_iter()
            .map(|s| {
                let mut bytes = s.as_ref().as_bytes().to_vec();
                for (start, end) in dust_unchecked(&bytes, opts) {
                    bytes[start..end].fill(letter);
                }
                String::from_utf8_lossy(&bytes).into_owned()
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn masks_overlapping_and_reverse_complement_matches() {
        let seqs = ["ACGTACGTTACG", "aaaaa"];
        let m = mask_seqs(&seqs, "ACG", false, '-', &Alphabet::Dna).unwrap();
        assert_eq!(m, ["---T---TT---", "aaaaa"]);
        let m = mask_seqs(&seqs, "aa", false, '-', &Alphabet::Dna).unwrap();
        assert_eq!(m[1], "-----");
        let m = mask_seqs(&seqs, "CGT", true, 'N', &Alphabet::Dna).unwrap();
        assert_eq!(m[0], "ANNNANNNTNNN");
        assert!(mask_seqs(&seqs, "KV", true, 'X', &Alphabet::Aa).is_err());
        assert_eq!(
            mask_ranges("ACGTACGT", &[(0, 2), (1, 3), (7, 8)], 'N').unwrap(),
            "NNNTACGN"
        );
        assert!(mask_ranges("ACGT", &[(2, 5)], 'N').is_err());
    }

    const DUST_A: &str = "ACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAACGTTGCAGTCGATCGATGCATGCATCAGTCAGCTAGCTAGCATCGACTAGCTAGCATGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTGACCAGT";
    const DUST_B: &str = "GGTTGTACTGTACAGGAACTTATAGAATAATTCTGGATGTATATATATATATATATATATATATATATATTCCAAACGAGCCTCCCGATAGTATCGAAATGGGGGGGGGGGGGATTTCCGTAGTTATACGTCTGCCAA";
    const DUST_C: &str = "AGTAATGCAGTCCATACTAGTGTACATCCTTTCCGCCCACTTGCCACTTCCAGCATCAGCAGCATCAGCAGCAGCATCAGCAGTTGGTCGCCTGCTGGAGCGCTGTCGGTCTTGTTGGATTCGGGATGGCGTG";

    #[test]
    fn dust_matches_dustmasker() {
        // Reference intervals from dustmasker 1.0.0 (BLAST 2.16.0), whose
        // inclusive ends are one less than these.
        let opts = DustOptions::default();
        assert_eq!(dust_intervals(DUST_A, &opts).unwrap(), [(8, 33), (50, 127)]);
        assert_eq!(
            dust_intervals(DUST_B, &opts).unwrap(),
            [(39, 70), (100, 113)]
        );
        assert_eq!(dust_intervals(DUST_C, &opts).unwrap(), [(50, 83)]);
        let opts = DustOptions {
            window: 32,
            threshold: 10,
            ..DustOptions::default()
        };
        assert_eq!(dust_intervals(DUST_A, &opts).unwrap(), [(8, 33), (63, 127)]);
        // Case and U do not change the regions.
        let lower = DUST_B.to_lowercase().replace('t', "u");
        assert_eq!(
            dust_intervals(&lower, &DustOptions::default()).unwrap(),
            [(39, 70), (100, 113)]
        );
    }

    #[test]
    fn dust_splits_at_other_letters_and_masks_lowercase() {
        let opts = DustOptions::default();
        // The halves of an AT repeat either side of Ns are scored apart, and
        // the Ns are left as they are.
        let seq = format!("{}NNN{}", &DUST_B[..55], &DUST_B[55..]);
        assert_eq!(
            dust_intervals(&seq, &opts).unwrap(),
            [(39, 55), (58, 73), (103, 116)]
        );
        assert!(
            dust_intervals("ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGT", &opts)
                .unwrap()
                .is_empty()
        );
        // Lowercase runs join the DUST regions, merging where they touch.
        let soft = format!("acgt{}gattaca", &DUST_B[4..120]);
        let opts = DustOptions {
            lowercase: true,
            ..DustOptions::default()
        };
        assert_eq!(
            dust_intervals(&soft, &opts).unwrap(),
            [(0, 4), (39, 70), (100, 113), (120, 127)]
        );
        let masked = mask_complexity(&[DUST_B, "acgt"], &opts).unwrap();
        assert_eq!(&masked[0][36..73], "ATG-------------------------------TCC");
        assert_eq!(masked[1], "----");
        let opts = DustOptions {
            letter: 'N',
            ..DustOptions::default()
        };
        assert_eq!(
            &mask_complexity(&[DUST_C], &opts).unwrap()[0][48..85],
            "TCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTT"
        );
        assert!(
            dust_intervals(
                "ACGT",
                &DustOptions {
                    window: 3,
                    ..DustOptions::default()
                }
            )
            .is_err()
        );
        assert!(
            dust_intervals(
                "ACGT",
                &DustOptions {
                    threshold: 0,
                    ..DustOptions::default()
                }
            )
            .is_err()
        );
        assert!(
            mask_complexity(
                &["ACGT"],
                &DustOptions {
                    letter: 'é',
                    ..DustOptions::default()
                }
            )
            .is_err()
        );
    }
}
