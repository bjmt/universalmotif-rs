//! De novo motif discovery (`motif_finder()`).
//!
//! A port of universalmotif's `src/motif_finder.cpp`, itself a port of
//! `yamme` from yamtk. For each motif width, k-mer seeds enriched in the
//! positive sequences (by Fisher's exact test on per-sequence presence) are
//! grown into PPMs from their Hamming-1 neighbours, refined by rescanning,
//! and accepted if hits are enriched over the negative sequences. Accepted
//! motifs mask their sites before the next search. Motifs found at different
//! widths that cover the same sites are deduplicated, Benjamini-Hochberg
//! q-values are computed, and low-information flanks are trimmed.

use std::collections::{HashMap, HashSet};

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::create::{CreateOptions, from_matrix};
use crate::error::{Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};
use crate::pvalue::with_threads;
use crate::shuffle::{ShuffleMethod, ShuffleOptions, shuffle_sequences};

const HAMMING_MISMATCH: u64 = 1;
const REFINE_PASSES: usize = 2;
const TOP_K_SEEDS: usize = 4;
const MAX_MOTIF_WIDTH: usize = 30;
const MIN_REFINE_HITS: u64 = 20;
const MIN_IC_BITS: f64 = 0.5;
const AMBIGUITY_SCORE: i32 = -10_000_000;
const MAX_CDF_SIZE: u64 = 2_097_152;
const COMP4: [usize; 4] = [3, 2, 1, 0];

const fn char_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}
const CHAR2IDX: [u8; 256] = char_table();

#[derive(Clone, Default)]
struct FMotif {
    pwm: Vec<i32>,
    pwm_rc: Vec<i32>,
    probs: Vec<[f64; 4]>,
    cdf: Vec<f64>,
    threshold: i32,
    mn: i32,
    mx: i32,
    max_score: i32,
    min_score: i32,
    cdf_max: i32,
    cdf_offset: i32,
    cdf_size: u64,
    size: usize,
    nsites: u64,
}

impl FMotif {
    fn init(&mut self, w: usize) {
        *self = FMotif {
            pwm: vec![0; w * 5],
            pwm_rc: vec![0; w * 5],
            probs: vec![[0.0; 4]; w],
            size: w,
            ..Default::default()
        };
        for i in 0..w {
            self.pwm[i * 5 + 4] = AMBIGUITY_SCORE;
            self.pwm_rc[i * 5 + 4] = AMBIGUITY_SCORE;
        }
    }

    fn score(&self, idx: usize, pos: usize) -> i32 {
        self.pwm[idx + pos * 5]
    }
}

fn calc_score_ns(prob: f64, bkg: f64, ns: u64, pc: u32) -> i32 {
    let mut x = prob * ns as f64;
    x += pc as f64 / 4.0;
    x /= (ns + pc as u64) as f64;
    if x <= 0.0 {
        x = 1e-12;
    }
    ((x / bkg).log2() * 1000.0) as i32
}

fn fill_cdf(m: &mut FMotif, bkg: &[f64; 4]) -> bool {
    if m.cdf_size == 0 || m.cdf_size > MAX_CDF_SIZE {
        return false;
    }
    let size = m.cdf_size as usize;
    let cmax = m.cdf_max as usize;
    m.cdf = vec![1.0; size];
    let mut tmp = vec![0.0; size];
    for i in 0..m.size {
        let max_step = i * cmax;
        tmp.copy_from_slice(&m.cdf);
        let end = (max_step + cmax + 1).min(size);
        m.cdf[..end].fill(0.0);
        for (j, b) in bkg.iter().enumerate() {
            let s = (m.score(j, i) - m.mn) as usize;
            for (c, t) in m.cdf[s..=s + max_step].iter_mut().zip(&tmp) {
                *c += t * b;
            }
        }
    }
    let sum: f64 = m.cdf.iter().fold(0.0, |a, &b| a + b);
    if (sum - 1.0).abs() > 0.0001 {
        m.cdf.iter_mut().for_each(|v| *v /= sum);
    }
    for i in (0..size - 1).rev() {
        m.cdf[i] += m.cdf[i + 1];
    }
    true
}

fn set_threshold_pval(m: &mut FMotif, pval: f64) {
    let threshold_i = m
        .cdf
        .iter()
        .position(|&c| c < pval)
        .unwrap_or(m.cdf_size as usize) as i32;
    m.threshold = threshold_i + m.mn * m.size as i32;
    m.max_score = 0;
    m.min_score = 0;
    for i in 0..m.size {
        let col: Vec<i32> = (0..4).map(|j| m.score(j, i)).collect();
        m.max_score += *col.iter().max().unwrap();
        m.min_score += *col.iter().min().unwrap();
    }
    let min_pval = m.cdf[(m.max_score - m.cdf_offset) as usize];
    if min_pval / pval > 1.0001 {
        m.threshold = i32::MAX;
    }
}

fn convert_ppm_to_motif(
    m: &mut FMotif,
    ppm: &[[u64; 4]],
    nsites: u64,
    bkg: &[f64; 4],
    pc: u32,
    hit_pval: f64,
) {
    let w = ppm.len();
    m.init(w);
    m.nsites = nsites;
    for (j, col) in ppm.iter().enumerate() {
        let mut total = col.iter().sum::<u64>() as f64;
        if total <= 0.0 {
            total = 1.0;
        }
        for i in 0..4 {
            let raw = col[i] as f64 / total;
            m.probs[j][i] = (raw * nsites as f64 + pc as f64 / 4.0) / (nsites + pc as u64) as f64;
            m.pwm[i + j * 5] = calc_score_ns(raw, bkg[i], nsites, pc);
        }
    }
    m.mn = (0..w)
        .flat_map(|p| (0..4).map(move |l| (p, l)))
        .fold(0, |a, (p, l)| a.min(m.score(l, p)));
    m.mx = (0..w)
        .flat_map(|p| (0..4).map(move |l| (p, l)))
        .fold(0, |a, (p, l)| a.max(m.score(l, p)));
    m.cdf_offset = m.mn * w as i32;
    for pos in 0..w {
        for (idx, &comp) in COMP4.iter().enumerate() {
            m.pwm_rc[idx + (w - 1 - pos) * 5] = m.score(comp, pos);
        }
    }
    m.cdf_max = m.mx - m.mn;
    m.cdf_size = w as u64 * m.cdf_max as u64 + 1;
    if m.cdf_max == 0 || m.cdf_size > MAX_CDF_SIZE || !fill_cdf(m, bkg) {
        m.threshold = i32::MAX;
        return;
    }
    set_threshold_pval(m, hit_pval);
    if m.threshold == i32::MAX {
        set_threshold_pval(m, 0.01);
    }
}

fn kmer_encode(seq: &[u8], w: usize, off: usize) -> Option<u64> {
    let mut k = 0u64;
    for &c in &seq[off..off + w] {
        let i = CHAR2IDX[c as usize];
        if i >= 4 {
            return None;
        }
        k = k * 4 + i as u64;
    }
    Some(k)
}

fn rc_kmer(mut kmer: u64, w: usize) -> u64 {
    let mut rc = 0;
    for _ in 0..w {
        rc = rc * 4 + COMP4[(kmer & 3) as usize] as u64;
        kmer >>= 2;
    }
    rc
}

fn hamming(mut a: u64, mut b: u64, w: usize) -> u64 {
    let mut mm = 0;
    for _ in 0..w {
        if a & 3 != b & 3 {
            mm += 1;
        }
        a >>= 2;
        b >>= 2;
    }
    mm
}

fn log_choose(n: f64, k: f64) -> f64 {
    libm::lgamma(n + 1.0) - libm::lgamma(k + 1.0) - libm::lgamma(n - k + 1.0)
}

/// One-sided (greater) Fisher's exact test P-value for a 2x2 table.
fn fisher_greater(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let (r1, r2) = ((a + b) as f64, (c + d) as f64);
    let c1 = (a + c) as f64;
    let n = (a + b + c + d) as f64;
    let denom = log_choose(n, c1);
    let x_max = r1.min(c1) as u64;
    let mut log_sum = log_choose(r1, a as f64) + log_choose(r2, c1 - a as f64) - denom;
    for x in a + 1..=x_max {
        let term = log_choose(r1, x as f64) + log_choose(r2, c1 - x as f64) - denom;
        if term > log_sum {
            log_sum = term + (log_sum - term).exp().ln_1p();
        } else {
            log_sum += (term - log_sum).exp().ln_1p();
        }
    }
    log_sum.exp().min(1.0)
}

fn bh(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| p[a].total_cmp(&p[b]));
    let mut q = vec![1.0; n];
    let mut prev = 1.0f64;
    for i in (0..n).rev() {
        let qv = (p[order[i]] * n as f64 / (i + 1) as f64).min(prev);
        prev = qv;
        q[order[i]] = qv.min(1.0);
    }
    q
}

fn score_at(m: &FMotif, seq: &[u8], off: usize) -> (i32, i32) {
    let (mut s, mut src) = (0, 0);
    for i in 0..m.size {
        let idx = CHAR2IDX[seq[off + i] as usize] as usize;
        s += m.pwm[idx + i * 5];
        src += m.pwm_rc[idx + i * 5];
    }
    (s, src)
}

fn enumerate_seeds(pos: &[Vec<u8>], neg: &[Vec<u8>], w: usize, rc: bool) -> Vec<(u64, f64)> {
    let count = |seqs: &[Vec<u8>]| -> HashMap<u64, u64> {
        let mut h: HashMap<u64, u64> = HashMap::new();
        for s in seqs {
            if s.len() < w {
                continue;
            }
            let mut seen = HashSet::new();
            for off in 0..=s.len() - w {
                if let Some(k) = kmer_encode(s, w, off) {
                    seen.insert(if rc { k.min(rc_kmer(k, w)) } else { k });
                }
            }
            for k in seen {
                *h.entry(k).or_default() += 1;
            }
        }
        h
    };
    let (ph, nh) = (count(pos), count(neg));
    let (np, nn) = (pos.len() as u64, neg.len() as u64);
    let mut cands: Vec<(u64, f64)> = ph
        .iter()
        .map(|(&k, &pc)| {
            let nc = nh.get(&k).copied().unwrap_or(0);
            (
                k,
                fisher_greater(pc, nc, np.saturating_sub(pc), nn.saturating_sub(nc)),
            )
        })
        .collect();
    // Ties are broken by k-mer code rather than hash order, so results are
    // reproducible.
    cands.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    cands.truncate(TOP_K_SEEDS);
    cands
}

fn add_site(ppm: &mut [[u64; 4]], seq: &[u8], off: usize, w: usize, reverse: bool) {
    for i in 0..w {
        if reverse {
            let idx = CHAR2IDX[seq[off + w - 1 - i] as usize] as usize;
            if idx < 4 {
                ppm[i][COMP4[idx]] += 1;
            }
        } else {
            let idx = CHAR2IDX[seq[off + i] as usize] as usize;
            if idx < 4 {
                ppm[i][idx] += 1;
            }
        }
    }
}

fn build_ppm_from_seed(pos: &[Vec<u8>], seed: u64, w: usize, rc: bool) -> (Vec<[u64; 4]>, u64) {
    let rc_seed = rc.then(|| rc_kmer(seed, w));
    let mut ppm = vec![[0u64; 4]; w];
    let mut nsites = 0;
    for s in pos {
        if s.len() < w {
            continue;
        }
        for off in 0..=s.len() - w {
            let Some(k) = kmer_encode(s, w, off) else {
                continue;
            };
            if hamming(k, seed, w) <= HAMMING_MISMATCH {
                add_site(&mut ppm, s, off, w, false);
                nsites += 1;
            } else if rc_seed.is_some_and(|r| hamming(k, r, w) <= HAMMING_MISMATCH) {
                add_site(&mut ppm, s, off, w, true);
                nsites += 1;
            }
        }
    }
    (ppm, nsites)
}

fn refine(
    pos: &[Vec<u8>],
    m: &mut FMotif,
    w: usize,
    rc: bool,
    bkg: &[f64; 4],
    pc: u32,
    hit_pval: f64,
) -> bool {
    if m.threshold == i32::MAX {
        return false;
    }
    let thr = m.threshold - 1;
    let mut ppm = vec![[0u64; 4]; w];
    let mut nsites = 0;
    for s in pos {
        if s.len() < w {
            continue;
        }
        for off in 0..=s.len() - w {
            let (sc, src) = score_at(m, s, off);
            if sc > thr {
                add_site(&mut ppm, s, off, w, false);
                nsites += 1;
            }
            if rc && src > thr {
                add_site(&mut ppm, s, off, w, true);
                nsites += 1;
            }
        }
    }
    if nsites < MIN_REFINE_HITS {
        return false;
    }
    convert_ppm_to_motif(m, &ppm, nsites, bkg, pc, hit_pval);
    true
}

struct Evaluation {
    pvalue: f64,
    sites_pos: u64,
    sites_neg: u64,
    seqs_pos: u64,
    seqs_neg: u64,
    covered: Vec<Vec<bool>>,
}

fn evaluate(pos: &[Vec<u8>], neg: &[Vec<u8>], m: &FMotif, rc: bool) -> Evaluation {
    let thr = m.threshold - 1;
    let w = m.size;
    let mut covered: Vec<Vec<bool>> = pos.iter().map(|s| vec![false; s.len()]).collect();
    let (mut sites_pos, mut seqs_pos, mut sites_neg, mut seqs_neg) = (0, 0, 0, 0);
    for (si, s) in pos.iter().enumerate() {
        if s.len() < w {
            continue;
        }
        let mut hit = false;
        for off in 0..=s.len() - w {
            let (sc, src) = score_at(m, s, off);
            for found in [sc > thr, rc && src > thr] {
                if found {
                    hit = true;
                    sites_pos += 1;
                    covered[si][off..off + w].fill(true);
                }
            }
        }
        seqs_pos += hit as u64;
    }
    for s in neg {
        if s.len() < w {
            continue;
        }
        let mut hit = false;
        for off in 0..=s.len() - w {
            let (sc, src) = score_at(m, s, off);
            for found in [sc > thr, rc && src > thr] {
                if found {
                    hit = true;
                    sites_neg += 1;
                }
            }
        }
        seqs_neg += hit as u64;
    }
    let (np, nn) = (pos.len() as u64, neg.len() as u64);
    Evaluation {
        pvalue: fisher_greater(
            seqs_pos,
            seqs_neg,
            np.saturating_sub(seqs_pos),
            nn.saturating_sub(seqs_neg),
        ),
        sites_pos,
        sites_neg,
        seqs_pos,
        seqs_neg,
        covered,
    }
}

struct Found {
    motif: FMotif,
    eval: Evaluation,
    width: usize,
    order: usize,
}

#[allow(clippy::too_many_arguments)]
fn discover_for_width(
    pos: &[Vec<u8>],
    neg: &[Vec<u8>],
    w: usize,
    nmotifs: usize,
    rc: bool,
    bkg: &[f64; 4],
    pc: u32,
    hit_pval: f64,
    stop_pval: f64,
) -> Vec<Found> {
    let mut work = pos.to_vec();
    let mut out = Vec::new();
    for _ in 0..nmotifs {
        let seeds = enumerate_seeds(&work, neg, w, rc);
        if seeds.is_empty() {
            break;
        }
        let mut accepted = false;
        for &(seed, _) in &seeds {
            let (ppm, nsites) = build_ppm_from_seed(&work, seed, w, rc);
            if nsites < MIN_REFINE_HITS {
                continue;
            }
            let mut m = FMotif::default();
            convert_ppm_to_motif(&mut m, &ppm, nsites, bkg, pc, hit_pval);
            if m.threshold == i32::MAX {
                continue;
            }
            for _ in 0..REFINE_PASSES {
                if !refine(&work, &mut m, w, rc, bkg, pc, hit_pval) {
                    break;
                }
            }
            if m.threshold == i32::MAX {
                continue;
            }
            let eval = evaluate(&work, neg, &m, rc);
            if eval.pvalue > stop_pval {
                continue;
            }
            for (s, c) in work.iter_mut().zip(&eval.covered) {
                for (b, &cov) in s.iter_mut().zip(c) {
                    if cov {
                        *b = b'N';
                    }
                }
            }
            let order = out.len();
            out.push(Found {
                motif: m,
                eval,
                width: w,
                order,
            });
            accepted = true;
            break;
        }
        if !accepted {
            break;
        }
    }
    out
}

fn column_ic(p: &[f64; 4]) -> f64 {
    let ic = p
        .iter()
        .filter(|&&x| x > 0.0)
        .fold(2.0, |acc, &x| acc + x * x.log2());
    ic.max(0.0)
}

/// Options for [`motif_finder`].
#[derive(Debug, Clone, PartialEq)]
pub struct FinderOptions {
    pub min_width: usize,
    pub max_width: usize,
    /// Most motifs searched for at each width.
    pub nmotifs: usize,
    /// P-value for a site to count as a hit.
    pub hit_pvalue: f64,
    /// A motif is kept if its hits are enriched with at most this P-value.
    pub stop_pvalue: f64,
    /// Motifs with q-values above this are dropped.
    pub qvalue: f64,
    /// Fraction of shared site positions above which the less significant
    /// of two motifs is dropped.
    pub dedup_overlap: f64,
    /// Search both strands.
    pub rc: bool,
    /// Shuffling for negative sequences when none are given.
    pub shuffle_k: usize,
    pub seed: u64,
    pub pseudocount: u32,
    pub nthreads: usize,
}

impl Default for FinderOptions {
    fn default() -> Self {
        FinderOptions {
            min_width: 6,
            max_width: 15,
            nmotifs: 10,
            hit_pvalue: 1e-4,
            stop_pvalue: 1e-3,
            qvalue: 1e-3,
            dedup_overlap: 0.5,
            rc: true,
            shuffle_k: 2,
            seed: 1,
            pseudocount: 1,
            nthreads: 1,
        }
    }
}

/// A motif found by [`motif_finder`], with its discovery statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct FoundMotif {
    pub motif: Motif,
    /// Rank by P-value (1-based).
    pub rank: usize,
    /// Positive and negative sequences with at least one hit.
    pub seqs_pos: u64,
    pub seqs_neg: u64,
    /// Hits in positive and negative sequences.
    pub sites_pos: u64,
    pub sites_neg: u64,
    pub n_pos: usize,
    pub n_neg: usize,
    pub pvalue: f64,
    pub qvalue: f64,
}

/// The background used by [`motif_finder`]: A, C, G and T (or U)
/// frequencies over all sequences, floored at 0.001.
pub fn finder_background<S: AsRef<str>>(pos: &[S], neg: &[S]) -> [f64; 4] {
    let mut counts = [0f64; 4];
    for s in pos.iter().chain(neg) {
        for b in s.as_ref().bytes() {
            if CHAR2IDX[b as usize] < 4 {
                counts[CHAR2IDX[b as usize] as usize] += 1.0;
            }
        }
    }
    let total: f64 = counts.iter().sum();
    let mut bkg = counts.map(|c| if total > 0.0 { c / total } else { 0.25 });
    for b in bkg.iter_mut() {
        if *b < 1e-3 {
            *b = 1e-3;
        }
    }
    let s: f64 = bkg.iter().sum();
    bkg.map(|b| b / s)
}

/// Find motifs enriched in DNA or RNA sequences relative to negative
/// sequences (by default the positives shuffled, keeping `shuffle_k`-lets).
/// Motifs come back as PPMs named `motif_1`, `motif_2` and so on, in order
/// of P-value, with the q-value stored as their E-value.
pub fn motif_finder<S: AsRef<str> + Sync>(
    sequences: &[S],
    negatives: Option<&[S]>,
    opts: &FinderOptions,
) -> Result<Vec<FoundMotif>> {
    if sequences.is_empty() {
        return arg_err("no sequences given");
    }
    if opts.min_width < 3 || opts.max_width > MAX_MOTIF_WIDTH || opts.min_width > opts.max_width {
        return arg_err("widths must satisfy 3 <= min_width <= max_width <= 30");
    }
    if opts.nmotifs < 1 {
        return arg_err("nmotifs must be at least 1");
    }
    if !(opts.hit_pvalue > 0.0 && opts.hit_pvalue < 1.0)
        || !(opts.stop_pvalue > 0.0 && opts.stop_pvalue <= 1.0)
        || !(opts.qvalue > 0.0 && opts.qvalue <= 1.0)
        || !(0.0..=1.0).contains(&opts.dedup_overlap)
    {
        return arg_err(
            "hit_pvalue must be in (0, 1); stop_pvalue and qvalue in (0, 1]; dedup_overlap in [0, 1]",
        );
    }
    let rna = sequences.iter().any(|s| s.as_ref().contains(['U', 'u']));
    let alphabet = if rna { Alphabet::Rna } else { Alphabet::Dna };
    let shuffled: Vec<String>;
    let neg_strs: Vec<&str> = match negatives {
        Some(n) => n.iter().map(|s| s.as_ref()).collect(),
        None => {
            shuffled = shuffle_sequences(
                sequences,
                &ShuffleOptions {
                    k: opts.shuffle_k,
                    method: ShuffleMethod::Euler,
                    seed: opts.seed,
                    window: None,
                    nthreads: opts.nthreads,
                },
            )?;
            shuffled.iter().map(|s| s.as_str()).collect()
        }
    };
    if neg_strs.is_empty() {
        return arg_err("no negative sequences");
    }
    let pos_strs: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();
    let bkg = finder_background(&pos_strs, &neg_strs);
    let pos: Vec<Vec<u8>> = pos_strs.iter().map(|s| s.as_bytes().to_vec()).collect();
    let neg: Vec<Vec<u8>> = neg_strs.iter().map(|s| s.as_bytes().to_vec()).collect();

    let widths: Vec<usize> = (opts.min_width..=opts.max_width).collect();
    let per_width: Vec<Vec<Found>> = with_threads(opts.nthreads, || {
        widths
            .par_iter()
            .map(|&w| {
                discover_for_width(
                    &pos,
                    &neg,
                    w,
                    opts.nmotifs,
                    opts.rc,
                    &bkg,
                    opts.pseudocount,
                    opts.hit_pvalue,
                    opts.stop_pvalue,
                )
            })
            .collect()
    })?;
    let mut found: Vec<Found> = per_width.into_iter().flatten().collect();
    found.sort_by_key(|f| (f.width, f.order));
    if found.is_empty() {
        return Ok(Vec::new());
    }
    // Re-evaluate on the unmasked sequences so coverage is comparable.
    for f in found.iter_mut() {
        f.eval = evaluate(&pos, &neg, &f.motif, opts.rc);
    }
    found.sort_by(|a, b| a.eval.pvalue.total_cmp(&b.eval.pvalue));
    let mut dropped = vec![false; found.len()];
    for i in 0..found.len() {
        if dropped[i] {
            continue;
        }
        for j in i + 1..found.len() {
            if dropped[j] {
                continue;
            }
            let (mut inter, mut si, mut sj) = (0u64, 0u64, 0u64);
            for (ci, cj) in found[i].eval.covered.iter().zip(&found[j].eval.covered) {
                for (&a, &b) in ci.iter().zip(cj) {
                    inter += (a && b) as u64;
                    si += a as u64;
                    sj += b as u64;
                }
            }
            let smin = si.min(sj);
            if smin > 0 && inter as f64 / smin as f64 > opts.dedup_overlap {
                dropped[j] = true;
            }
        }
    }
    let survivors: Vec<Found> = found
        .into_iter()
        .zip(dropped)
        .filter(|(_, d)| !d)
        .map(|(f, _)| f)
        .collect();
    let q = bh(&survivors.iter().map(|f| f.eval.pvalue).collect::<Vec<_>>());

    let mut out = Vec::new();
    for (i, (f, qv)) in survivors.iter().zip(q).enumerate() {
        // Trim low-information flanks, keeping at least three columns.
        let probs = &f.motif.probs;
        let mut left = 0;
        while left < probs.len() && column_ic(&probs[left]) < MIN_IC_BITS {
            left += 1;
        }
        let mut right = probs.len();
        while right > left && column_ic(&probs[right - 1]) < MIN_IC_BITS {
            right -= 1;
        }
        let kept = if right - left < 3 {
            probs.as_slice()
        } else {
            &probs[left..right]
        };
        let cols: Vec<Vec<f64>> = kept.iter().map(|c| c.to_vec()).collect();
        if qv > opts.qvalue {
            continue;
        }
        let copts = CreateOptions {
            motif_type: MotifType::Ppm,
            name: format!("motif_{}", i + 1),
            nsites: Some(f.motif.nsites as f64),
            pseudocount: opts.pseudocount as f64,
            bkg: Some(bkg.to_vec()),
            pval: Some(f.eval.pvalue),
            eval: Some(qv),
            ..Default::default()
        };
        let motif = from_matrix(
            Matrix::from_columns(&cols)?,
            None,
            Some(alphabet.clone()),
            &copts,
        )?;
        out.push(FoundMotif {
            motif,
            rank: out.len() + 1,
            seqs_pos: f.eval.seqs_pos,
            seqs_neg: f.eval.seqs_neg,
            sites_pos: f.eval.sites_pos,
            sites_neg: f.eval.sites_neg,
            n_pos: pos.len(),
            n_neg: neg.len(),
            pvalue: f.eval.pvalue,
            qvalue: qv,
        });
    }
    Ok(out)
}
