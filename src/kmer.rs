//! k-lets, k-let counts and background frequencies of sequences.

use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::pvalue::with_threads;
use crate::scan::reverse_complement_seq;

/// All k-lets of `letters`, in lexicographic order of the letters given, so the
/// first position varies slowest.
pub fn klets(letters: &[char], k: usize) -> Vec<String> {
    let mut out = vec![String::new()];
    for _ in 0..k {
        out = out
            .iter()
            .flat_map(|p| letters.iter().map(move |l| format!("{p}{l}")))
            .collect();
    }
    out
}

/// Sorted, unique letters of some sequences.
pub fn sequence_letters<S: AsRef<str>>(seqs: &[S]) -> Vec<char> {
    let mut ascii = [false; 128];
    let mut other = std::collections::BTreeSet::new();
    for s in seqs {
        let s = s.as_ref();
        if s.is_ascii() {
            s.bytes().for_each(|b| ascii[b as usize] = true);
        } else {
            for c in s.chars() {
                if c.is_ascii() {
                    ascii[c as usize] = true;
                } else {
                    other.insert(c);
                }
            }
        }
    }
    (0u8..128)
        .filter(|&b| ascii[b as usize])
        .map(char::from)
        .chain(other)
        .collect()
}

/// Counts in a sequence of each k-let of `alphabet`, in [`klets`] order. k-lets
/// containing letters outside the alphabet are not counted.
pub fn count_klets(seq: &str, k: usize, alphabet: &[char]) -> Vec<u64> {
    let n = alphabet.len();
    let mut counts = vec![0u64; n.pow(k as u32)];
    if k == 0 {
        return counts;
    }
    // A rolling code of the last k letters; `run` counts the letters since
    // the last one outside the alphabet.
    let keep = n.pow(k as u32 - 1);
    let (mut code, mut run) = (0usize, 0usize);
    let mut feed = |idx: Option<usize>| match idx {
        None => run = 0,
        Some(i) => {
            code = (code % keep) * n + i;
            run += 1;
            if run >= k {
                counts[code] += 1;
            }
        }
    };
    if seq.is_ascii() && alphabet.iter().all(char::is_ascii) {
        let mut table = [usize::MAX; 128];
        for (i, &a) in alphabet.iter().enumerate() {
            table[a as usize] = i;
        }
        for b in seq.bytes() {
            let i = table[b as usize];
            feed((i != usize::MAX).then_some(i));
        }
    } else {
        for c in seq.chars() {
            feed(alphabet.iter().position(|&a| a == c));
        }
    }
    counts
}

/// The part of `s` from 1-based letter `start` to `stop`, inclusive. For
/// non-ASCII strings, `bounds` holds the byte offset of every letter and of the
/// end.
fn letters<'a>(s: &'a str, bounds: Option<&[usize]>, start: usize, stop: usize) -> &'a str {
    match bounds {
        None => &s[start - 1..stop],
        Some(b) => &s[b[start - 1]..b[stop]],
    }
}

/// 1-based start and stop coordinates of sliding windows over a sequence of
/// length `len`. The last window may be shorter than `window`, unless
/// `incomplete` is false. Returns no windows
/// if `window` exceeds `len` or `overlap` is not below `window`.
pub fn windows(len: usize, window: usize, overlap: usize, incomplete: bool) -> Vec<(usize, usize)> {
    if window > len || overlap >= window || window == 0 {
        return Vec::new();
    }
    let mut out = vec![(1, window)];
    loop {
        let (start, stop) = *out.last().unwrap();
        if start == len {
            break;
        }
        let next_start = start + window - overlap;
        let next_stop = next_start + window - 1;
        if next_stop > len {
            let next_start = next_start.min(len);
            if len != stop {
                out.push((next_start, len));
            }
            break;
        }
        out.push((next_start, next_stop));
    }
    if !incomplete {
        // Only a last window shorter than `window` is dropped (see
        // DIVERGENCES.md).
        let (start, stop) = *out.last().unwrap();
        if stop - start + 1 < window {
            out.pop();
        }
    }
    out
}

/// Windowing for [`get_bkg`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowOptions {
    /// Window width; values below 1 are a fraction of each sequence's length.
    pub size: f64,
    /// Window overlap; values between 0 and 1 are a fraction of each
    /// sequence's length.
    pub overlap: f64,
}

/// Options for [`get_bkg`].
#[derive(Debug, Clone, PartialEq)]
pub struct BkgOptions {
    /// k-let sizes to count.
    pub k: Vec<usize>,
    /// Added to every count.
    pub pseudocount: u64,
    /// Letters to count; by default the letters found in the sequences.
    pub alphabet: Option<Vec<char>>,
    /// Also count the reverse complements, for DNA and RNA only.
    pub rc: bool,
    /// Combine all sequences into one set of counts.
    pub merge: bool,
    /// Count in sliding windows.
    pub window: Option<WindowOptions>,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for BkgOptions {
    fn default() -> Self {
        BkgOptions {
            k: vec![1, 2, 3],
            pseudocount: 0,
            alphabet: None,
            rc: false,
            merge: true,
            window: None,
            nthreads: 1,
        }
    }
}

/// One row of [`get_bkg`] output. `sequence`, `start` and `stop` are set only
/// for per-sequence or windowed counts.
#[derive(Debug, Clone, PartialEq)]
pub struct BkgRow {
    pub sequence: Option<usize>,
    pub start: Option<usize>,
    pub stop: Option<usize>,
    pub klet: String,
    pub count: u64,
    pub probability: f64,
}

fn probabilities(counts: &[u64]) -> Vec<f64> {
    let total: u64 = counts.iter().sum();
    counts
        .iter()
        .map(|&c| {
            if total == 0 {
                0.0
            } else {
                c as f64 / total as f64
            }
        })
        .collect()
}

/// k-let counts and frequencies of sequences.
///
/// Rows are ordered by k, k-let and sequence for per-sequence counts; by k-let
/// across all k and then window start for merged windows; and by sequence name,
/// k-let and window start for per-sequence windows. Names and k-lets are
/// compared bytewise, and `names` gives the sequence names to order by.
pub fn get_bkg<S: AsRef<str> + Sync>(
    seqs: &[S],
    names: &[String],
    opts: &BkgOptions,
) -> Result<Vec<BkgRow>> {
    if seqs.is_empty() {
        return arg_err("no sequences given");
    }
    if opts.k.is_empty() || opts.k.contains(&0) {
        return arg_err("k must be at least 1");
    }
    let mut ks = opts.k.clone();
    ks.sort_unstable();
    ks.dedup();
    // Borrowed unless reverse complements are added below.
    let mut seqs: Vec<std::borrow::Cow<str>> = seqs.iter().map(|s| s.as_ref().into()).collect();
    let alphabet = match &opts.alphabet {
        Some(a) => {
            let mut a = a.clone();
            a.sort_unstable();
            a.dedup();
            a
        }
        None => sequence_letters(&seqs),
    };
    if opts.rc {
        let rna = alphabet.contains(&'U');
        let rcs: Vec<String> = seqs
            .iter()
            .map(|s| reverse_complement_seq(s, rna))
            .collect();
        seqs.extend(rcs.into_iter().map(std::borrow::Cow::Owned));
    }
    let names: Vec<String> = if names.len() == seqs.len() {
        names.to_vec()
    } else {
        let base: Vec<String> = if names.is_empty() {
            (1..=seqs.len()).map(|i| i.to_string()).collect()
        } else {
            names.to_vec()
        };
        (0..seqs.len())
            .map(|i| base[i % base.len()].clone())
            .collect()
    };
    let kmax = *ks.last().unwrap();
    let min_len = seqs.iter().map(|s| s.chars().count()).min().unwrap();
    if min_len < kmax {
        return arg_err("k must not be greater than the length of the shortest sequence");
    }

    let Some(win) = opts.window else {
        let mut rows = Vec::new();
        for &k in &ks {
            let kl = klets(&alphabet, k);
            let per_seq: Vec<Vec<u64>> = with_threads(opts.nthreads, || {
                seqs.par_iter()
                    .map(|s| count_klets(s, k, &alphabet))
                    .collect()
            })?;
            if opts.merge {
                let counts: Vec<u64> = (0..kl.len())
                    .map(|j| per_seq.iter().map(|c| c[j]).sum::<u64>() + opts.pseudocount)
                    .collect();
                let probs = probabilities(&counts);
                for (j, klet) in kl.into_iter().enumerate() {
                    rows.push(BkgRow {
                        sequence: None,
                        start: None,
                        stop: None,
                        klet,
                        count: counts[j],
                        probability: probs[j],
                    });
                }
            } else {
                let counts: Vec<Vec<u64>> = per_seq
                    .iter()
                    .map(|c| c.iter().map(|x| x + opts.pseudocount).collect())
                    .collect();
                let probs: Vec<Vec<f64>> = counts.iter().map(|c| probabilities(c)).collect();
                for (j, klet) in kl.iter().enumerate() {
                    for (si, c) in counts.iter().enumerate() {
                        rows.push(BkgRow {
                            sequence: Some(si),
                            start: None,
                            stop: None,
                            klet: klet.clone(),
                            count: c[j],
                            probability: probs[si][j],
                        });
                    }
                }
            }
        }
        return Ok(rows);
    };

    // Windowed counts.
    let mut wins: Vec<Vec<(usize, usize)>> = Vec::with_capacity(seqs.len());
    for s in &seqs {
        let len = s.chars().count();
        let size = if win.size < 1.0 {
            (win.size * len as f64) as usize
        } else {
            win.size as usize
        };
        let overlap = if win.overlap > 0.0 && win.overlap < 1.0 {
            (win.overlap * len as f64) as usize
        } else {
            win.overlap as usize
        };
        if size == 0 {
            return arg_err("the window size is too small for these sequences (resolved to 0)");
        }
        if size <= overlap {
            return arg_err("the window overlap must be smaller than the window size");
        }
        if size < kmax {
            return arg_err("the window size must be at least k");
        }
        let w = windows(len, size, overlap, true);
        if w.is_empty() {
            return arg_err("incorrect window and/or overlap size");
        }
        wins.push(w);
    }
    if opts.merge && seqs.iter().any(|s| s.chars().count() != min_len) {
        return arg_err("merged windowed counts need all sequences to be of equal length");
    }
    let bounds: Vec<Option<Vec<usize>>> = seqs
        .iter()
        .map(|s| {
            (!s.is_ascii()).then(|| {
                s.char_indices()
                    .map(|(i, _)| i)
                    .chain(std::iter::once(s.len()))
                    .collect()
            })
        })
        .collect();
    if wins
        .iter()
        .flatten()
        .any(|&(start, stop)| stop - start + 1 < kmax)
    {
        return arg_err("k must not be greater than the smallest window");
    }

    let mut rows = Vec::new();
    for &k in &ks {
        let kl = klets(&alphabet, k);
        // counts[seq][window][klet]
        let counts: Vec<Vec<Vec<u64>>> = with_threads(opts.nthreads, || {
            seqs.par_iter()
                .zip(&bounds)
                .zip(&wins)
                .map(|((s, b), ws)| {
                    ws.iter()
                        .map(|&(start, stop)| {
                            count_klets(letters(s, b.as_deref(), start, stop), k, &alphabet)
                        })
                        .collect()
                })
                .collect()
        })?;
        if opts.merge {
            let ws = &wins[0];
            let summed: Vec<Vec<u64>> = (0..ws.len())
                .map(|w| {
                    (0..kl.len())
                        .map(|j| counts.iter().map(|c| c[w][j]).sum::<u64>() + opts.pseudocount)
                        .collect()
                })
                .collect();
            let probs: Vec<Vec<f64>> = summed.iter().map(|c| probabilities(c)).collect();
            for (j, klet) in kl.iter().enumerate() {
                for (w, &(start, stop)) in ws.iter().enumerate() {
                    rows.push(BkgRow {
                        sequence: None,
                        start: Some(start),
                        stop: Some(stop),
                        klet: klet.clone(),
                        count: summed[w][j],
                        probability: probs[w][j],
                    });
                }
            }
        } else {
            let mut order: Vec<usize> = (0..seqs.len()).collect();
            order.sort_by(|&a, &b| names[a].as_bytes().cmp(names[b].as_bytes()));
            for &si in &order {
                let c: Vec<Vec<u64>> = counts[si]
                    .iter()
                    .map(|w| w.iter().map(|x| x + opts.pseudocount).collect())
                    .collect();
                let probs: Vec<Vec<f64>> = c.iter().map(|w| probabilities(w)).collect();
                for (j, klet) in kl.iter().enumerate() {
                    for (w, &(start, stop)) in wins[si].iter().enumerate() {
                        rows.push(BkgRow {
                            sequence: Some(si),
                            start: Some(start),
                            stop: Some(stop),
                            klet: klet.clone(),
                            count: c[w][j],
                            probability: probs[w][j],
                        });
                    }
                }
            }
        }
    }
    if opts.merge {
        // Merged windows are ordered by k-let only, across all k.
        rows.sort_by(|a, b| a.klet.as_bytes().cmp(b.klet.as_bytes()));
    } else {
        rows.sort_by(|a, b| {
            names[a.sequence.unwrap()]
                .as_bytes()
                .cmp(names[b.sequence.unwrap()].as_bytes())
                .then(a.klet.as_bytes().cmp(b.klet.as_bytes()))
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn klets_counted_for_any_letters() {
        let ascii = count_klets("AACGNAAC", 2, &['A', 'C', 'G']);
        // AA, AC, CG, then N breaks the run, then AA, AC.
        assert_eq!(ascii, [2, 2, 0, 0, 0, 1, 0, 0, 0]);
        let greek = count_klets("αβαNββ", 2, &['α', 'β']);
        assert_eq!(greek, [0, 1, 1, 1]);
        let bounds: Vec<usize> = "αβγ".char_indices().map(|(i, _)| i).chain([6]).collect();
        assert_eq!(letters("αβγ", Some(&bounds), 2, 3), "βγ");
    }
}
