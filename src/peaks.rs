//! Positional enrichment of motif hits, tested by [`motif_peaks`].
//!
//! This is a CentriMo-style test (Bailey and Machanick 2012). If hit centres
//! were spread uniformly along sequences of equal length, the number falling in
//! a window would be binomial. Each motif's best hit per sequence is kept,
//! windows of several widths are tested, centred on the middle of the sequences
//! or anywhere along them, and the best window's binomial P-value is
//! Bonferroni-corrected for the number of windows tried.

use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::pvalue::with_threads;
use crate::stats::{p_adjust_bh, pbinom_upper};

/// Which windows [`motif_peaks`] tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeaksMode {
    /// Windows centred on the middle of the sequences.
    Central,
    /// Windows anywhere along the sequences.
    Local,
}

impl PeaksMode {
    /// The name [`str::parse`] reads, such as `"central"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            PeaksMode::Central => "central",
            PeaksMode::Local => "local",
        }
    }
}

impl std::str::FromStr for PeaksMode {
    type Err = crate::error::Error;

    /// Parse `"central"` or `"local"`.
    fn from_str(s: &str) -> Result<PeaksMode> {
        match s {
            "central" => Ok(PeaksMode::Central),
            "local" => Ok(PeaksMode::Local),
            _ => arg_err(format!("unknown mode '{s}'; expected 'central' or 'local'")),
        }
    }
}

impl std::fmt::Display for PeaksMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Options for [`motif_peaks`].
#[derive(Debug, Clone, PartialEq)]
pub struct PeaksOptions {
    /// Where windows are placed; centred on the sequences by default.
    pub mode: PeaksMode,
    /// Narrowest window; 10 by default.
    pub min_window: usize,
    /// Widest window; by default half the sequence length.
    pub max_window: Option<usize>,
    /// Step between window widths; 10 by default.
    pub window_step: usize,
    /// Step between window centres in local mode.
    pub position_step: usize,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for PeaksOptions {
    fn default() -> Self {
        PeaksOptions {
            mode: PeaksMode::Central,
            min_window: 10,
            max_window: None,
            window_step: 10,
            position_step: 5,
            nthreads: 1,
        }
    }
}

/// A motif hit for [`motif_peaks`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeakHit {
    /// Motif identifier, of any value; results come in increasing order of it.
    pub motif: usize,
    /// Sequence identifier, of any value.
    pub sequence: usize,
    /// Centre of the hit, `(start + end - 1) / 2` for a hit covering
    /// `start..end`, which on a sequence read on the minus strand is
    /// reflected to `length - 1 - centre`.
    pub centre: f64,
    /// The hit's score; only the best hit of a motif in each sequence is used.
    pub score: f64,
}

/// The best window for one motif.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct MotifPeak {
    /// The motif, as given in the hits.
    pub motif: usize,
    /// Sequences with a hit; only the best hit per sequence is used.
    pub nhits: usize,
    /// Width and centre of the best window; `None` if no window had a
    /// P-value below 1.
    pub best_window: Option<usize>,
    /// Centre of the best window.
    pub best_centre: Option<f64>,
    /// Hits inside the best window.
    pub hits_in: usize,
    /// Hits outside the best window.
    pub hits_out: usize,
    /// Hits expected in the best window if they were uniform.
    pub expected_in: Option<f64>,
    /// The observed hits in the best window divided by the expected number.
    pub enrichment: Option<f64>,
    /// The base 2 logarithm of the enrichment.
    pub log2_enrichment: Option<f64>,
    /// Bonferroni-corrected binomial P-value of the best window.
    pub pvalue: f64,
    /// Benjamini-Hochberg q-value across motifs.
    pub qvalue: f64,
    /// Centres of the hits used, rounded half to even, in input order.
    pub centres: Vec<i64>,
}

fn test_motif(
    motif: usize,
    centres: &[f64],
    len: usize,
    widths: &[usize],
    opts: &PeaksOptions,
) -> MotifPeak {
    let n = centres.len();
    let lenf = len as f64;
    // Positions count from 0, so the middle of the sequence, taken as half
    // its length when counting from 1, sits one position lower.
    let mid = lenf / 2.0 - 1.0;
    let (mut best_p, mut best_w, mut best_c, mut best_k) = (1.0, None, None, 0);
    let mut n_tests = 0usize;
    let count = |c: f64, half: f64| centres.iter().filter(|&&x| (x - c).abs() <= half).count();
    for &w in widths {
        let half = w as f64 / 2.0;
        let p = w as f64 / lenf;
        let cs: Vec<f64> = match opts.mode {
            PeaksMode::Central => vec![mid],
            PeaksMode::Local => {
                let (lo, hi) = (half.ceil() as i64, (lenf - half).floor() as i64);
                if hi < lo {
                    continue;
                }
                (lo..=hi)
                    .step_by(opts.position_step)
                    .map(|c| c as f64 - 1.0)
                    .collect()
            }
        };
        for c in cs {
            let k = count(c, half);
            let pv = pbinom_upper(k as u64, n as u64, p);
            n_tests += 1;
            if pv < best_p {
                (best_p, best_w, best_c, best_k) = (pv, Some(w), Some(c), k);
            }
        }
    }
    let pvalue = (best_p * n_tests.max(1) as f64).min(1.0);
    let expected_in = best_w.map(|w| n as f64 * w as f64 / lenf);
    let enrichment = expected_in.filter(|&e| e > 0.0).map(|e| best_k as f64 / e);
    MotifPeak {
        motif,
        nhits: n,
        best_window: best_w,
        best_centre: best_c,
        hits_in: best_k,
        hits_out: n - best_k,
        expected_in,
        enrichment,
        log2_enrichment: enrichment.filter(|&e| e > 0.0).map(f64::log2),
        pvalue,
        qvalue: 0.0,
        // Rounded as the positions would be when counting from 1, with ties
        // going to the even one.
        centres: centres
            .iter()
            .map(|c| (c + 1.0).round_ties_even() as i64 - 1)
            .collect(),
    }
}

/// Test whether each motif's hits cluster along sequences of length
/// `seq_length`. Returns one result per motif with hits, in increasing order of motif, unfiltered.
pub fn motif_peaks(
    hits: &[PeakHit],
    seq_length: usize,
    opts: &PeaksOptions,
) -> Result<Vec<MotifPeak>> {
    if seq_length < 2 {
        return arg_err("the sequence length must be at least 2");
    }
    if opts.min_window < 1 || opts.window_step < 1 || opts.position_step < 1 {
        return arg_err("min_window, window_step and position_step must be at least 1");
    }
    let max_window = opts.max_window.unwrap_or(seq_length / 2).min(seq_length);
    if max_window < opts.min_window {
        return arg_err(format!(
            "the widest window ({max_window}) is narrower than min_window ({})",
            opts.min_window
        ));
    }
    let widths: Vec<usize> = (opts.min_window..=max_window)
        .step_by(opts.window_step)
        .collect();

    // The best-scoring hit (the first, for ties) per motif and sequence.
    let mut best: Vec<(usize, usize, usize)> = Vec::new(); // (motif, sequence, hit)
    let mut index: std::collections::HashMap<(usize, usize), usize> = Default::default();
    for (i, h) in hits.iter().enumerate() {
        match index.get(&(h.motif, h.sequence)) {
            Some(&b) => {
                if h.score > hits[best[b].2].score {
                    best[b].2 = i;
                }
            }
            None => {
                index.insert((h.motif, h.sequence), best.len());
                best.push((h.motif, h.sequence, i));
            }
        }
    }
    let mut kept: Vec<usize> = best.iter().map(|b| b.2).collect();
    kept.sort_unstable();
    let mut motif_ids: Vec<usize> = kept.iter().map(|&i| hits[i].motif).collect();
    motif_ids.sort_unstable();
    motif_ids.dedup();
    let groups: Vec<(usize, Vec<f64>)> = motif_ids
        .iter()
        .map(|&m| {
            let c = kept
                .iter()
                .filter(|&&i| hits[i].motif == m)
                .map(|&i| hits[i].centre)
                .collect();
            (m, c)
        })
        .collect();
    let mut out: Vec<MotifPeak> = with_threads(opts.nthreads, || {
        groups
            .par_iter()
            .map(|(m, c)| test_motif(*m, c, seq_length, &widths, opts))
            .collect()
    })?;
    let q = p_adjust_bh(&out.iter().map(|p| p.pvalue).collect::<Vec<_>>());
    for (p, q) in out.iter_mut().zip(q) {
        p.qvalue = q;
    }
    Ok(out)
}
