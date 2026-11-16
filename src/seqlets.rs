//! Calling seqlets, the short stretches of high contribution in a sequence
//! model's attribution scores, with [`call_seqlets`].
//!
//! A stretch is a seqlet when its total contribution is unlikely under a
//! null built from all the scores given, and so is the total of every
//! stretch inside it of at least the minimum length, so that a seqlet holds
//! no run of unremarkable positions. The null for each length comes from
//! the scores' own distribution, convolved with itself, so a set share of
//! stretches passes even in pure noise; to correct for this, the same calls
//! are made on the scores shuffled within each sequence, and seqlets are
//! kept down to the total at which these decoys make up at most a given
//! share of those passing.

use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::progress::{self, Progress};
use crate::pvalue::with_threads;

/// Which signs of contribution to call seqlets for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SeqletSign {
    /// Positive and negative seqlets, each from its own null.
    #[default]
    Both,
    /// Stretches that raise the prediction.
    Positive,
    /// Stretches that lower it.
    Negative,
}

/// Options for [`call_seqlets`].
#[derive(Debug, Clone, PartialEq)]
pub struct SeqletOptions {
    /// The P-value every stretch inside a seqlet must reach; 0.01 by default.
    pub threshold: f64,
    /// The shortest seqlet, and the shortest stretch whose P-value is
    /// checked inside one; 4 by default.
    pub min_len: usize,
    /// The longest seqlet; 25 by default.
    pub max_len: usize,
    /// Positions added on either side of each seqlet, within its sequence;
    /// 3 by default.
    pub flank: usize,
    /// With a value, keep seqlets down to the total at which decoys are at
    /// most this share of those passing; 0.2 by default. Decoys are called
    /// in the decoy scores given to [`call_seqlets`], or else in each
    /// sequence's scores shuffled. `None` keeps every seqlet.
    pub fdr: Option<f64>,
    /// The signs to call.
    pub sign: SeqletSign,
    /// The bins the scores are counted in to build the null; 1000 by
    /// default.
    pub n_bins: usize,
    /// Seed for the shuffles.
    pub seed: u64,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for SeqletOptions {
    fn default() -> Self {
        SeqletOptions {
            threshold: 0.01,
            min_len: 4,
            max_len: 25,
            flank: 3,
            fdr: Some(0.2),
            sign: SeqletSign::Both,
            n_bins: 1000,
            seed: 1,
            nthreads: 1,
        }
    }
}

/// A seqlet.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Seqlet {
    /// Index of the sequence.
    pub sequence: usize,
    /// First position, flanks included, counting from 0.
    pub start: usize,
    /// Position just past the seqlet, flanks included.
    pub end: usize,
    /// First position of the seqlet without its flanks.
    pub core_start: usize,
    /// Position just past the seqlet without its flanks.
    pub core_end: usize,
    /// Whether the seqlet raises the prediction.
    pub positive: bool,
    /// The total contribution of the seqlet without its flanks, negative for
    /// negative seqlets.
    pub score: f64,
    /// The P-value of that total under the null for its length.
    pub pvalue: f64,
}

/// Call seqlets in the observed contribution scores of each sequence, one
/// score per position. Missing scores (NaN) count as 0.
///
/// Seqlets are called from the longest down: a stretch from `min_len` to
/// `max_len` long is taken when it and every stretch inside it of at least
/// `min_len` positions have a P-value below `threshold`, in order of P-value,
/// unless it overlaps a seqlet already taken. Seqlets come back in order of
/// sequence and position.
///
/// `decoys`, the scores of sequences with nothing to find, such as a model's
/// contributions for shuffled sequences, set the false discovery rate when
/// given: seqlets are called in them the same way, and their number, scaled
/// to as many positions as the scores have, stands for the false seqlets
/// among those called. Otherwise each sequence's own scores, shuffled,
/// stand in for them, which keeps their values but not the model's view of
/// sequence without sites.
///
/// The scores can be any slices of `f64` or `f32`, such as `Vec<f64>` or
/// rows borrowed from a larger array; each sequence's are copied only while it is called, so the
/// memory used beyond the scores stays small.
pub fn call_seqlets<T: Copy + Into<f64> + Sync, S: AsRef<[T]> + Sync>(
    scores: &[S],
    decoys: Option<&[S]>,
    opts: &SeqletOptions,
) -> Result<Vec<Seqlet>> {
    if !(opts.threshold > 0.0 && opts.threshold <= 1.0) {
        return arg_err("threshold must be above 0 and at most 1");
    }
    if opts.min_len == 0 || opts.min_len > opts.max_len {
        return arg_err("min_len must be at least 1 and at most max_len");
    }
    if opts.n_bins < 2 {
        return arg_err("n_bins must be at least 2");
    }
    if let Some(f) = opts.fdr
        && !(f > 0.0 && f <= 1.0)
    {
        return arg_err("fdr must be above 0 and at most 1");
    }
    let signs: &[f64] = match opts.sign {
        SeqletSign::Both => &[1.0, -1.0],
        SeqletSign::Positive => &[1.0],
        SeqletSign::Negative => &[-1.0],
    };
    let positions = |v: &[S]| v.iter().map(|s| s.as_ref().len()).sum::<usize>();
    let total = positions(scores);
    let work = progress::current();
    progress::begin(&work, "calling seqlets", (total * signs.len()) as u64);
    // Each sequence's scores are signed, cleaned and, for decoys, shuffled
    // only as it is called, so that no copy of all the scores is kept.
    let out = with_threads(opts.nthreads, || -> Result<Vec<Seqlet>> {
        let mut out = Vec::new();
        for &sign in signs {
            let Some(null) = Null::new(scores, sign, opts.max_len, opts.n_bins) else {
                continue;
            };
            let mut real = call_all(
                scores.len(),
                |i| signed(scores[i].as_ref(), sign),
                &null,
                opts,
                &work,
            )?;
            if let Some(fdr) = opts.fdr {
                // Decoy scores given, such as the model's on shuffled
                // sequences, or else each sequence's scores shuffled.
                let (called, n_decoy) = match decoys {
                    Some(d) => (
                        call_all(d.len(), |i| signed(d[i].as_ref(), sign), &null, opts, &None)?,
                        positions(d),
                    ),
                    None => {
                        let mut rng =
                            StdRng::seed_from_u64(opts.seed ^ if sign > 0.0 { 0 } else { 0x5eed });
                        let seeds: Vec<u64> = (0..scores.len()).map(|_| rng.next_u64()).collect();
                        let shuffled = |i: usize| {
                            let mut v = signed(scores[i].as_ref(), sign);
                            v.shuffle(&mut StdRng::seed_from_u64(seeds[i]));
                            v
                        };
                        (call_all(scores.len(), shuffled, &null, opts, &None)?, total)
                    }
                };
                // Decoys are counted as if there were as many positions of
                // them as of the scores.
                let scale = total.max(1) as f64 / n_decoy.max(1) as f64;
                let cut = fdr_cut(
                    &real.iter().map(|s| s.score).collect::<Vec<_>>(),
                    &called.iter().map(|s| s.score).collect::<Vec<_>>(),
                    fdr,
                    scale,
                );
                real.retain(|s| s.score >= cut);
            }
            for s in &mut real {
                s.positive = sign > 0.0;
                s.score *= sign;
            }
            out.extend(real);
        }
        Ok(out)
    })??;
    let mut out = out;
    out.sort_by_key(|s| (s.sequence, s.core_start, !s.positive));
    Ok(out)
}

/// A sequence's scores times `sign`, with missing ones as zero.
fn signed<T: Copy + Into<f64>>(s: &[T], sign: f64) -> Vec<f64> {
    s.iter()
        .map(|&v| {
            let v: f64 = v.into();
            if v.is_finite() { v * sign } else { 0.0 }
        })
        .collect()
}

/// The lowest total at which decoys at or above it, times `scale`, are at
/// most `fdr` of the real seqlets at or above it, or infinity when there is
/// none.
fn fdr_cut(real: &[f64], decoy: &[f64], fdr: f64, scale: f64) -> f64 {
    let mut v = real.to_vec();
    v.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let mut d = decoy.to_vec();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut cut = f64::INFINITY;
    for (i, &t) in v.iter().enumerate() {
        let n_decoy = d.len() - d.partition_point(|&x| x < t);
        if n_decoy as f64 * scale / (i + 1) as f64 <= fdr {
            cut = t;
        }
    }
    cut
}

/// The null distribution of the total of a stretch of each length, from the
/// binned scores convolved with themselves.
struct Null {
    min: f64,
    width: f64,
    /// For each length, the chance that a stretch's binned total is at least
    /// each value.
    surv: Vec<Vec<f64>>,
}

impl Null {
    fn new<T: Copy + Into<f64>, S: AsRef<[T]>>(
        x: &[S],
        sign: f64,
        max_len: usize,
        n_bins: usize,
    ) -> Option<Null> {
        let all = || {
            x.iter().flat_map(|s| s.as_ref().iter()).map(move |&v| {
                let v: f64 = v.into();
                if v.is_finite() { v * sign } else { 0.0 }
            })
        };
        let (lo, hi) = all().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| {
            (a.min(v), b.max(v))
        });
        // Empty, constant or all-missing scores give no null.
        if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        let width = (hi - lo) / (n_bins - 1) as f64;
        let mut f = vec![0.0; n_bins];
        let mut n = 0.0;
        for v in all() {
            f[(((v - lo) / width) as usize).min(n_bins - 1)] += 1.0;
            n += 1.0;
        }
        for v in &mut f {
            *v /= n;
        }
        let mut surv = vec![Vec::new(), survival(&f)];
        let mut pmf = f.clone();
        for _ in 2..=max_len {
            pmf = convolve(&pmf, &f);
            surv.push(survival(&pmf));
        }
        Some(Null {
            min: lo,
            width,
            surv,
        })
    }

    /// The P-value of a stretch of `len` positions totalling `sum`.
    fn pvalue(&self, len: usize, sum: f64) -> f64 {
        let s = &self.surv[len];
        let b = ((sum - self.min * len as f64) / self.width).floor();
        if b < 0.0 {
            1.0
        } else {
            s.get(b as usize).copied().unwrap_or(0.0)
        }
    }
}

fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    out.par_iter_mut().enumerate().for_each(|(k, o)| {
        let lo = k.saturating_sub(b.len() - 1);
        let hi = k.min(a.len() - 1);
        *o = (lo..=hi).map(|i| a[i] * b[k - i]).sum();
    });
    out
}

fn survival(pmf: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; pmf.len()];
    let mut acc = 0.0;
    for i in (0..pmf.len()).rev() {
        acc += pmf[i];
        out[i] = acc.min(1.0);
    }
    out
}

fn call_all<F: Fn(usize) -> Vec<f64> + Sync>(
    n: usize,
    get: F,
    null: &Null,
    opts: &SeqletOptions,
    work: &Option<Progress>,
) -> Result<Vec<Seqlet>> {
    let calls: Vec<Result<Vec<Seqlet>>> = (0..n)
        .into_par_iter()
        .map(|i| {
            progress::check(work)?;
            let s = get(i);
            let found = call_one(i, &s, null, opts);
            progress::add(work, s.len() as u64);
            Ok(found)
        })
        .collect();
    let mut out = Vec::new();
    for c in calls {
        out.extend(c?);
    }
    Ok(out)
}

/// The seqlets of one sequence.
fn call_one(seq: usize, x: &[f64], null: &Null, opts: &SeqletOptions) -> Vec<Seqlet> {
    let n = x.len();
    let (lo, hi) = (opts.min_len, opts.max_len.min(n));
    if n < lo {
        return Vec::new();
    }
    let mut csum = vec![0.0; n + 1];
    for i in 0..n {
        csum[i + 1] = csum[i] + x[i];
    }
    // own[len][k]: the P-value of the stretch at k of that length. worst[len][k]:
    // the largest P-value of any stretch of at least lo positions inside it,
    // from prefixes (worst over own[len-1][k]) and suffixes (worst[len-1][k+1]).
    let width = hi - lo + 1;
    let mut own = vec![Vec::new(); width];
    let mut worst = vec![Vec::new(); width];
    for (j, len) in (lo..=hi).enumerate() {
        let m = n - len + 1;
        own[j] = (0..m)
            .map(|k| null.pvalue(len, csum[k + len] - csum[k]))
            .collect();
        worst[j] = if j == 0 {
            own[0].clone()
        } else {
            (0..m)
                .map(|k| own[j][k].max(worst[j - 1][k]).max(worst[j - 1][k + 1]))
                .collect()
        };
    }
    let mut used = vec![false; n];
    let mut out = Vec::new();
    for j in (0..width).rev() {
        let len = lo + j;
        let mut cands: Vec<usize> = (0..worst[j].len())
            .filter(|&k| worst[j][k] < opts.threshold)
            .collect();
        cands.sort_by(|&a, &b| own[j][a].partial_cmp(&own[j][b]).unwrap().then(a.cmp(&b)));
        for k in cands {
            if used[k..k + len].iter().any(|&u| u) {
                continue;
            }
            used[k..k + len].iter_mut().for_each(|u| *u = true);
            out.push(Seqlet {
                sequence: seq,
                start: k.saturating_sub(opts.flank),
                end: (k + len + opts.flank).min(n),
                core_start: k,
                core_end: k + len,
                positive: true,
                score: csum[k + len] - csum[k],
                pvalue: own[j][k],
            });
        }
    }
    out
}
