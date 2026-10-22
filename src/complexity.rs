//! Sequence complexity, of whole sequences with [`complexity`] and in sliding
//! windows with [`sequence_complexity`].

use std::collections::BTreeSet;

use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::kmer::windows;
use crate::pvalue::with_threads;

/// Complexity measure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    /// Wootton and Federhen's complexity, the probability of a sequence's
    /// letter composition.
    #[default]
    WoottonFederhen,
    /// The entropy approximation of Wootton and Federhen's complexity.
    WoottonFederhenFast,
    /// Trifonov's linguistic complexity: the product over word sizes of the
    /// fraction of possible words seen.
    Trifonov,
    /// The sum of words seen over the sum of possible words.
    TrifonovFast,
    /// The DUST score of repeated trinucleotides, for DNA only.
    Dust,
}

impl Method {
    pub fn parse(s: &str) -> Result<Method> {
        Ok(match s {
            "WoottonFederhen" => Method::WoottonFederhen,
            "WoottonFederhenFast" => Method::WoottonFederhenFast,
            "Trifonov" => Method::Trifonov,
            "TrifonovFast" => Method::TrifonovFast,
            "DUST" => Method::Dust,
            _ => {
                return arg_err(format!(
                    "unknown complexity method '{s}' (use WoottonFederhen, WoottonFederhenFast, \
                     Trifonov, TrifonovFast or DUST)"
                ));
            }
        })
    }
}

fn letter_counts(seq: &[u8], alph: &[u8]) -> Vec<f64> {
    alph.iter()
        .map(|a| seq.iter().filter(|c| *c == a).count() as f64)
        .collect()
}

fn distinct_words(seq: &[u8], n: usize) -> f64 {
    seq.windows(n).collect::<BTreeSet<_>>().len() as f64
}

/// The complexity of one sequence. `alph` gives the possible letters, by
/// default those in the sequence, and `max_word` limits Trifonov word sizes.
pub fn complexity(seq: &str, method: Method, alph: Option<&str>, max_word: usize) -> f64 {
    let s = seq.as_bytes();
    let alph: Vec<u8> = match alph {
        Some(a) if !a.is_empty() => a.bytes().collect(),
        _ => s
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect(),
    };
    let len = s.len() as f64;
    let n = alph.len() as f64;
    match method {
        Method::WoottonFederhen => {
            // log(L! / prod(n_i!)) with log-gamma, so long sequences do not
            // overflow.
            let counts = letter_counts(s, &alph);
            let log_pi = libm::lgamma(len + 1.0)
                - counts.iter().map(|&c| libm::lgamma(c + 1.0)).sum::<f64>();
            (log_pi / n.ln()) / len
        }
        Method::WoottonFederhenFast => letter_counts(s, &alph)
            .iter()
            .filter(|&&c| c > 0.0)
            .fold(0.0, |acc, &c| acc - (c / len) * ((c / len).ln() / n.ln())),
        Method::Trifonov | Method::TrifonovFast => {
            let maxw = max_word.min(s.len());
            let (mut seen, mut possible) = (Vec::new(), Vec::new());
            for i in 0..maxw {
                seen.push(distinct_words(s, i + 1));
                possible.push(f64::min(n.powi(i as i32 + 1), len - i as f64));
            }
            if method == Method::Trifonov {
                seen.iter()
                    .zip(&possible)
                    .fold(1.0, |acc, (a, b)| acc * (a / b))
            } else {
                seen.iter().sum::<f64>() / possible.iter().sum::<f64>()
            }
        }
        Method::Dust => {
            let l = len - 2.0;
            let mut counts: std::collections::BTreeMap<&[u8], f64> = Default::default();
            for w in s.windows(3) {
                *counts.entry(w).or_default() += 1.0;
            }
            counts.values().map(|c| c * (c - 1.0) / 2.0).sum::<f64>() / (l - 1.0)
        }
    }
}

/// Options for [`sequence_complexity`].
#[derive(Debug, Clone, PartialEq)]
pub struct ComplexityOptions {
    pub method: Method,
    /// Window width; values below 1 are a fraction of each sequence's length.
    pub window: f64,
    /// Window overlap; by default half the window. Values between 0 and 1
    /// are a fraction of each sequence's length.
    pub overlap: Option<f64>,
    /// Possible letters; by default those in the sequences.
    pub alphabet: Option<String>,
    pub max_word: usize,
    pub nthreads: usize,
}

impl Default for ComplexityOptions {
    fn default() -> Self {
        ComplexityOptions {
            method: Method::WoottonFederhen,
            window: 20.0,
            overlap: None,
            alphabet: None,
            max_word: 7,
            nthreads: 1,
        }
    }
}

/// Complexity of one window.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowComplexity {
    pub sequence: usize,
    pub start: usize,
    pub stop: usize,
    pub complexity: f64,
}

/// Complexity in sliding windows along sequences. Windows longer than a
/// sequence shrink to it.
pub fn sequence_complexity<S: AsRef<str> + Sync>(
    seqs: &[S],
    opts: &ComplexityOptions,
) -> Result<Vec<WindowComplexity>> {
    let mut letters: BTreeSet<u8> = seqs.iter().flat_map(|s| s.as_ref().bytes()).collect();
    if let Some(a) = &opts.alphabet {
        letters.extend(a.bytes());
    }
    let alph: String = letters.into_iter().map(|b| b as char).collect();
    if alph.len() <= 1 {
        return arg_err(format!(
            "can't calculate complexity when the alphabet is a single letter [{alph}]"
        ));
    }
    let default_overlap = (opts.window / 2.0).round_ties_even();
    let mut jobs = Vec::new();
    for (si, s) in seqs.iter().enumerate() {
        let len = s.as_ref().len();
        let (mut size, mut overlap) = (opts.window, opts.overlap.unwrap_or(default_overlap));
        if size > len as f64 {
            overlap = 0.0;
            size = len as f64;
        }
        if size == 0.0 {
            return arg_err("the window size must be greater than 0");
        }
        let size = if size < 1.0 {
            (size * len as f64) as usize
        } else {
            size as usize
        };
        let overlap = if overlap > 0.0 && overlap < 1.0 {
            (overlap * len as f64) as usize
        } else {
            overlap as usize
        };
        if size <= overlap {
            return arg_err("the window overlap must be smaller than the window size");
        }
        let max_word = opts.max_word.min(size);
        for (start, stop) in windows(len, size, overlap, true) {
            jobs.push((si, start, stop, max_word));
        }
    }
    with_threads(opts.nthreads, || {
        jobs.par_iter()
            .map(|&(si, start, stop, max_word)| WindowComplexity {
                sequence: si,
                start,
                stop,
                complexity: complexity(
                    &seqs[si].as_ref()[start - 1..stop],
                    opts.method,
                    Some(&alph),
                    max_word,
                ),
            })
            .collect()
    })
}
