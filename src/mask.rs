//! Masking sequences by pattern with [`mask_seqs`] and by position with
//! [`mask_ranges`].

use crate::alphabet::Alphabet;
use crate::error::{Result, arg_err};
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

/// Replace the positions in `ranges` of `seq` with `letter`, where ranges are
/// 1-based and inclusive. Ranges may overlap; one reaching outside the sequence
/// is an error.
pub fn mask_ranges(seq: &str, ranges: &[(i64, i64)], letter: char) -> Result<String> {
    let letter = check_letter(letter)?;
    let mut bytes = seq.as_bytes().to_vec();
    let n = bytes.len() as i64;
    for &(start, end) in ranges {
        if end < start {
            continue;
        }
        if start < 1 || end > n {
            return arg_err(format!(
                "the range {start}-{end} is outside the sequence (length {n})"
            ));
        }
        bytes[start as usize - 1..end as usize].fill(letter);
    }
    Ok(String::from_utf8_lossy(&bytes).into_owned())
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
            mask_ranges("ACGTACGT", &[(1, 2), (2, 3), (8, 8)], 'N').unwrap(),
            "NNNTACGN"
        );
        assert!(mask_ranges("ACGT", &[(3, 5)], 'N').is_err());
    }
}
