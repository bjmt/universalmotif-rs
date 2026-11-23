//! De novo motif discovery with [`motif_finder`], a minimal STREME-like
//! search (Bailey 2021).
//!
//! For each motif width, k-mer seeds enriched in the positive sequences, by
//! Fisher's exact test on per-sequence presence, are grown into PPMs from their
//! Hamming-1 neighbours, refined by rescanning, and accepted if hits are
//! enriched over the negative sequences. Accepted motifs mask their sites
//! before the next search. Motifs found at different widths that cover the same
//! sites are deduplicated, Benjamini-Hochberg q-values are computed, and
//! low-information flanks are trimmed.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::create::{CreateOptions, from_matrix};
use crate::error::{Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType, Strand};
use crate::progress::{self, Progress};
use crate::pvalue::with_threads;
use crate::shuffle::{ShuffleMethod, ShuffleOptions, shuffled_backgrounds};

const HAMMING_MISMATCH: u64 = 1;
const REFINE_PASSES: usize = 2;
const TOP_K_SEEDS: usize = 4;
const MAX_MOTIF_WIDTH: usize = 30;
const MIN_REFINE_HITS: u64 = 20;
const MIN_IC_BITS: f64 = 0.5;
const AMBIGUITY_SCORE: i32 = -10_000_000;
const MAX_CDF_SIZE: u64 = 2_097_152;
const MAX_SPACED_PAIRS: usize = 65_536;
const COMP4: [usize; 4] = [3, 2, 1, 0];
/// Sequences handed to a thread at a time when a search step is spread over
/// the sequences.
const PAR_CHUNK: usize = 512;

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

thread_local! {
    /// Buffers for score distributions, reused because they can take 16 MB
    /// each and differ in size from motif to motif, which leaves the system
    /// allocator holding on to many freed blocks if each is allocated anew.
    static CDF_BUFFERS: std::cell::RefCell<(Vec<f64>, Vec<f64>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

fn fill_cdf(m: &FMotif, bkg: &[f64; 4], cdf: &mut Vec<f64>, tmp: &mut Vec<f64>) -> bool {
    if m.cdf_size == 0 || m.cdf_size > MAX_CDF_SIZE {
        return false;
    }
    let size = m.cdf_size as usize;
    let cmax = m.cdf_max as usize;
    cdf.clear();
    cdf.resize(size, 1.0);
    tmp.clear();
    tmp.resize(size, 0.0);
    for i in 0..m.size {
        let max_step = i * cmax;
        tmp.copy_from_slice(cdf);
        let end = (max_step + cmax + 1).min(size);
        cdf[..end].fill(0.0);
        for (j, b) in bkg.iter().enumerate() {
            let s = (m.score(j, i) - m.mn) as usize;
            for (c, t) in cdf[s..=s + max_step].iter_mut().zip(tmp.iter()) {
                *c += t * b;
            }
        }
    }
    let sum: f64 = cdf.iter().fold(0.0, |a, &b| a + b);
    if (sum - 1.0).abs() > 0.0001 {
        cdf.iter_mut().for_each(|v| *v /= sum);
    }
    for i in (0..size - 1).rev() {
        cdf[i] += cdf[i + 1];
    }
    true
}

fn set_threshold_pval(m: &mut FMotif, cdf: &[f64], pval: f64) {
    let threshold_i = cdf
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
    let min_pval = cdf[(m.max_score - m.cdf_offset) as usize];
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
    CDF_BUFFERS.with(|b| {
        let (cdf, tmp) = &mut *b.borrow_mut();
        if m.cdf_max == 0 || m.cdf_size > MAX_CDF_SIZE || !fill_cdf(m, bkg, cdf, tmp) {
            m.threshold = i32::MAX;
            return;
        }
        set_threshold_pval(m, cdf, hit_pval);
        if m.threshold == i32::MAX {
            set_threshold_pval(m, cdf, 0.01);
        }
    });
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

/// The number of letters two k-mer codes differ at.
fn hamming(a: u64, b: u64) -> u64 {
    let x = a ^ b;
    ((x | (x >> 1)) & 0x5555_5555_5555_5555).count_ones() as u64
}

/// `ln(n!)` for every `n` up to the number of sequences, so Fisher tests need
/// no `lgamma` calls of their own; the values are the same.
struct LogFact(Vec<f64>);

impl LogFact {
    fn new(n: usize) -> LogFact {
        LogFact((0..=n).map(|i| libm::lgamma(i as f64 + 1.0)).collect())
    }

    fn choose(&self, n: u64, k: u64) -> f64 {
        self.0[n as usize] - self.0[k as usize] - self.0[(n - k) as usize]
    }
}

/// One-sided Fisher's exact test P-value for a 2x2 table, with the "greater"
/// alternative.
fn fisher_greater(lf: &LogFact, a: u64, b: u64, c: u64, d: u64) -> f64 {
    ln_fisher_greater(lf, a, b, c, d).exp()
}

/// The natural log of [`fisher_greater`], which stays finite for P-values too
/// small for a double.
fn ln_fisher_greater(lf: &LogFact, a: u64, b: u64, c: u64, d: u64) -> f64 {
    let (r1, r2) = (a + b, c + d);
    let c1 = a + c;
    let n = a + b + c + d;
    let denom = lf.choose(n, c1);
    let x_max = r1.min(c1);
    let mut log_sum = lf.choose(r1, a) + lf.choose(r2, c1 - a) - denom;
    for x in a + 1..=x_max {
        let term = lf.choose(r1, x) + lf.choose(r2, c1 - x) - denom;
        if term > log_sum {
            log_sum = term + (log_sum - term).exp().ln_1p();
        } else {
            let add = (term - log_sum).exp().ln_1p();
            // A term too small to change the sum, even eight times over, lies
            // well past the mode, where the terms only shrink, so none of the
            // rest would change it either.
            if log_sum + 8.0 * add == log_sum {
                break;
            }
            log_sum += add;
        }
    }
    log_sum.min(0.0)
}

/// A motif's scores for two positions at a time, which halves the lookups in
/// scoring a site. Scoring stops after the first half of the positions if
/// neither strand can pass with the most the rest could add. Scores are whole
/// numbers, so their sums are the same in any order.
struct Scorer<'a> {
    m: &'a FMotif,
    /// The scores on both strands of each pair of letters at each pair of
    /// positions, 25 to a pair of positions, at 5 times the first letter
    /// plus the second.
    pairs: Vec<[i32; 2]>,
    /// The pairs of positions scored before checking.
    check_at: usize,
    /// The most the positions after them can add on each strand.
    rest: [i32; 2],
}

impl Scorer<'_> {
    fn new(m: &FMotif) -> Scorer<'_> {
        let n = m.size / 2;
        let mut pairs = Vec::with_capacity(n * 25);
        for p in 0..n {
            let (i, j) = (2 * p * 5, (2 * p + 1) * 5);
            for a in 0..5 {
                for b in 0..5 {
                    pairs.push([
                        m.pwm[i + a] + m.pwm[j + b],
                        m.pwm_rc[i + a] + m.pwm_rc[j + b],
                    ]);
                }
            }
        }
        let check_at = n.div_ceil(2);
        let best = |pwm: &[i32]| -> i32 {
            (2 * check_at..m.size)
                .map(|i| *pwm[i * 5..i * 5 + 4].iter().max().unwrap())
                .sum()
        };
        Scorer {
            m,
            pairs,
            check_at,
            rest: [best(&m.pwm), best(&m.pwm_rc)],
        }
    }

    /// Whether the site at `off` scores above `thr` on the forward strand,
    /// and with `rc` on the reverse strand.
    #[inline(always)]
    fn hits_at(&self, seq: &[u8], off: usize, thr: i32, rc: bool) -> (bool, bool) {
        let (mut s, mut src) = (0, 0);
        let site = &seq[off..off + self.m.size];
        let pair =
            |p: usize| self.pairs[p * 25 + site[2 * p] as usize * 5 + site[2 * p + 1] as usize];
        for p in 0..self.check_at {
            let [f, r] = pair(p);
            s += f;
            src += r;
        }
        if s + self.rest[0] <= thr && (!rc || src + self.rest[1] <= thr) {
            return (false, false);
        }
        for p in self.check_at..self.m.size / 2 {
            let [f, r] = pair(p);
            s += f;
            src += r;
        }
        if self.m.size % 2 == 1 {
            let i = self.m.size - 1;
            let idx = site[i] as usize;
            s += self.m.pwm[idx + i * 5];
            src += self.m.pwm_rc[idx + i * 5];
        }
        (s > thr, rc && src > thr)
    }
}

/// The distinct k-mers of a sequence, or the smaller of each and its reverse
/// complement, written sorted to the start of `out`, which must have room for
/// every word; returns how many there are.
fn distinct_words(s: &[u8], w: usize, rc: bool, out: &mut [u64]) -> usize {
    let mask = if w == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * w)) - 1
    };
    let mut k = 0;
    // Rolling codes of the last w letters on each strand.
    let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
    for &idx in s {
        if idx >= 4 {
            run = 0;
            continue;
        }
        fwd = ((fwd << 2) | idx as u64) & mask;
        rev = (rev >> 2) | ((COMP4[idx as usize] as u64) << (2 * (w - 1)));
        run += 1;
        if run >= w {
            out[k] = if rc { fwd.min(rev) } else { fwd };
            k += 1;
        }
    }
    // Each word once.
    let words = &mut out[..k];
    words.sort_unstable();
    let mut n = 0;
    for i in 0..words.len() {
        if n == 0 || words[n - 1] != words[i] {
            words[n] = words[i];
            n += 1;
        }
    }
    n
}

/// Every k-mer, or the smaller of it and its reverse complement, found in the
/// sequences, sorted, with the number of sequences containing each.
fn kmer_presence(seqs: &[Vec<u8>], w: usize, rc: bool) -> (Vec<u64>, Vec<u32>) {
    // The distinct words of each sequence, gathered a chunk of sequences at a
    // time on the thread pool. Each chunk writes into its own part of one
    // buffer, sized for every word, so the words are never held twice.
    let room: Vec<usize> = seqs
        .chunks(PAR_CHUNK)
        .map(|chunk| chunk.iter().map(|s| (s.len() + 1).saturating_sub(w)).sum())
        .collect();
    let mut all: Vec<u64> = vec![0; room.iter().sum()];
    let mut parts: Vec<&mut [u64]> = Vec::with_capacity(room.len());
    let mut rest = all.as_mut_slice();
    for &r in &room {
        let (part, tail) = rest.split_at_mut(r);
        parts.push(part);
        rest = tail;
    }
    let used: Vec<usize> = seqs
        .par_chunks(PAR_CHUNK)
        .zip(parts.into_par_iter())
        .map(|(chunk, out)| {
            let mut k = 0;
            for s in chunk {
                k += distinct_words(s, w, rc, &mut out[k..]);
            }
            k
        })
        .collect();
    // Close the gaps the parts left.
    let (mut dst, mut src) = (0, 0);
    for (&u, &r) in used.iter().zip(&room) {
        all.copy_within(src..src + u, dst);
        dst += u;
        src += r;
    }
    all.truncate(dst);
    all.par_sort_unstable();
    // Collapse runs of the same k-mer in place, counting them, with the counts
    // sized for the distinct k-mers alone.
    let distinct = all.chunk_by(|a, b| a == b).count();
    let mut counts: Vec<u32> = Vec::with_capacity(distinct);
    let mut n = 0;
    for i in 0..all.len() {
        if n > 0 && all[n - 1] == all[i] {
            counts[n - 1] += 1;
        } else {
            all[n] = all[i];
            counts.push(1);
            n += 1;
        }
    }
    all.truncate(n);
    all.shrink_to_fit();
    (all, counts)
}

/// The best few seeds by P-value from the words of the positives, as from
/// [`kmer_presence`], where words in no sequence are left out, and the number
/// of words tested. Many k-mers share the same counts, so each table is
/// tested once, its P-value kept in `memo` for the searches that follow.
fn enumerate_seeds(
    pos_words: &(Vec<u64>, Vec<u32>),
    np: u64,
    neg_words: &(Vec<u64>, Vec<u32>),
    nn: u64,
    lf: &LogFact,
    memo: &mut HashMap<(u64, u64), f64>,
) -> (Vec<(u64, f64)>, u64) {
    let (pk, pn) = (&pos_words.0, &pos_words.1);
    let (nk, nn_) = (&neg_words.0, &neg_words.1);
    let mut best: Vec<(u64, f64)> = Vec::with_capacity(TOP_K_SEEDS + 1);
    let mut tested = 0;
    let mut j = 0;
    for (&k, &pc) in pk.iter().zip(pn) {
        if pc == 0 {
            continue;
        }
        tested += 1;
        let pc = pc as u64;
        while j < nk.len() && nk[j] < k {
            j += 1;
        }
        let nc = if j < nk.len() && nk[j] == k {
            nn_[j] as u64
        } else {
            0
        };
        let p = *memo.entry((pc, nc)).or_insert_with(|| {
            fisher_greater(lf, pc, nc, np.saturating_sub(pc), nn.saturating_sub(nc))
        });
        // Keep the best few by P-value, ties going to the lower k-mer code
        // (k-mers arrive in increasing order, so a tie never displaces).
        if best.len() < TOP_K_SEEDS || p < best[best.len() - 1].1 {
            let at = best.partition_point(|&(_, q)| q <= p);
            best.insert(at, (k, p));
            best.truncate(TOP_K_SEEDS);
        }
    }
    (best, tested)
}

fn add_site(ppm: &mut [[u64; 4]], seq: &[u8], off: usize, w: usize, reverse: bool) {
    for i in 0..w {
        if reverse {
            let idx = seq[off + w - 1 - i] as usize;
            if idx < 4 {
                ppm[i][COMP4[idx]] += 1;
            }
        } else {
            let idx = seq[off + i] as usize;
            if idx < 4 {
                ppm[i][idx] += 1;
            }
        }
    }
}

fn build_ppm_from_seed(pos: &[Vec<u8>], seed: u64, w: usize, rc: bool) -> (Vec<[u64; 4]>, u64) {
    let rc_seed = rc.then(|| rc_kmer(seed, w));
    let mask = (1u64 << (2 * w)) - 1;
    sum_sites(pos, w, |s, ppm, nsites| {
        // The rolling code of the last w letters.
        let (mut k, mut run) = (0u64, 0usize);
        for (i, &idx) in s.iter().enumerate() {
            if idx >= 4 {
                run = 0;
                continue;
            }
            k = ((k << 2) | idx as u64) & mask;
            run += 1;
            if run < w {
                continue;
            }
            let off = i + 1 - w;
            if hamming(k, seed) <= HAMMING_MISMATCH {
                add_site(ppm, s, off, w, false);
                *nsites += 1;
            } else if rc_seed.is_some_and(|r| hamming(k, r) <= HAMMING_MISMATCH) {
                add_site(ppm, s, off, w, true);
                *nsites += 1;
            }
        }
    })
}

/// The letter counts of the sites `visit` adds from each sequence, and how
/// many sites it added, summed over chunks of sequences on the thread pool.
/// The counts are whole numbers, so the order of the chunks doesn't matter.
fn sum_sites<F>(seqs: &[Vec<u8>], w: usize, visit: F) -> (Vec<[u64; 4]>, u64)
where
    F: Fn(&[u8], &mut [[u64; 4]], &mut u64) + Sync,
{
    seqs.par_chunks(PAR_CHUNK)
        .map(|chunk| {
            let mut ppm = vec![[0u64; 4]; w];
            let mut nsites = 0;
            for s in chunk {
                visit(s, &mut ppm, &mut nsites);
            }
            (ppm, nsites)
        })
        .reduce(
            || (vec![[0u64; 4]; w], 0),
            |(mut a, na), (b, nb)| {
                for (x, y) in a.iter_mut().zip(&b) {
                    for k in 0..4 {
                        x[k] += y[k];
                    }
                }
                (a, na + nb)
            },
        )
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
    let scored: &FMotif = m;
    let scorer = Scorer::new(scored);
    let (ppm, nsites) = sum_sites(pos, w, |s, ppm, nsites| {
        if s.len() < w {
            return;
        }
        for off in 0..=s.len() - w {
            let (fwd, rev) = scorer.hits_at(s, off, thr, rc);
            if fwd {
                add_site(ppm, s, off, w, false);
                *nsites += 1;
            }
            if rev {
                add_site(ppm, s, off, w, true);
                *nsites += 1;
            }
        }
    });
    if nsites < MIN_REFINE_HITS {
        return false;
    }
    convert_ppm_to_motif(m, &ppm, nsites, bkg, pc, hit_pval);
    true
}

struct Evaluation {
    /// Natural log of the P-value.
    ln_pvalue: f64,
    sites_pos: u64,
    sites_neg: u64,
    seqs_pos: u64,
    seqs_neg: u64,
    covered: Coverage,
}

/// The positions of a set of sequences covered by a motif's sites, a bit
/// each.
struct Coverage {
    /// Where each sequence's bits begin.
    starts: Vec<usize>,
    bits: Vec<u64>,
}

impl Coverage {
    fn new(seqs: &[Vec<u8>]) -> Coverage {
        let mut starts = Vec::with_capacity(seqs.len());
        let mut total = 0;
        for s in seqs {
            starts.push(total);
            total += s.len();
        }
        Coverage {
            starts,
            bits: vec![0; total.div_ceil(64)],
        }
    }

    fn fill(&mut self, seq: usize, from: usize, len: usize) {
        for p in self.starts[seq] + from..self.starts[seq] + from + len {
            self.bits[p / 64] |= 1 << (p % 64);
        }
    }

    fn contains(&self, seq: usize, pos: usize) -> bool {
        let p = self.starts[seq] + pos;
        self.bits[p / 64] & (1 << (p % 64)) != 0
    }

    fn count(&self) -> u64 {
        self.bits.iter().map(|w| w.count_ones() as u64).sum()
    }

    fn shared(&self, other: &Coverage) -> u64 {
        self.bits
            .iter()
            .zip(&other.bits)
            .map(|(a, b)| (a & b).count_ones() as u64)
            .sum()
    }
}

fn evaluate(pos: &[Vec<u8>], neg: &[Vec<u8>], m: &FMotif, rc: bool, lf: &LogFact) -> Evaluation {
    let thr = m.threshold - 1;
    let w = m.size;
    let scorer = Scorer::new(m);
    // The sites of one sequence, at most two at each start (one per strand).
    let sites_of = |s: &[u8], out: &mut Vec<usize>| {
        if s.len() < w {
            return;
        }
        for off in 0..=s.len() - w {
            let (fwd, rev) = scorer.hits_at(s, off, thr, rc);
            for found in [fwd, rev] {
                if found {
                    out.push(off);
                }
            }
        }
    };
    // Each chunk's sites in the sequences, by sequence; the background only
    // needs counting.
    let pos_sites: Vec<Vec<Vec<usize>>> = pos
        .par_chunks(PAR_CHUNK)
        .map(|chunk| {
            chunk
                .iter()
                .map(|s| {
                    let mut v = Vec::new();
                    sites_of(s, &mut v);
                    v
                })
                .collect()
        })
        .collect();
    let (sites_neg, seqs_neg) = neg
        .par_chunks(PAR_CHUNK)
        .map(|chunk| {
            let mut v = Vec::new();
            let (mut sites, mut seqs) = (0u64, 0u64);
            for s in chunk {
                v.clear();
                sites_of(s, &mut v);
                sites += v.len() as u64;
                seqs += !v.is_empty() as u64;
            }
            (sites, seqs)
        })
        .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    let mut covered = Coverage::new(pos);
    let (mut sites_pos, mut seqs_pos) = (0u64, 0u64);
    for (si, offs) in pos_sites.iter().flatten().enumerate() {
        sites_pos += offs.len() as u64;
        seqs_pos += !offs.is_empty() as u64;
        for &off in offs {
            covered.fill(si, off, w);
        }
    }
    let (np, nn) = (pos.len() as u64, neg.len() as u64);
    Evaluation {
        ln_pvalue: ln_fisher_greater(
            lf,
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
    /// Words tested in the search that found the motif.
    candidates: u64,
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
    lf: &LogFact,
    progress: &Option<Progress>,
) -> (Vec<Found>, u64) {
    let mut work = pos.to_vec();
    // The background is never masked, so its words are counted once.
    let neg_words = kmer_presence(neg, w, rc);
    // The positives' words are counted once as well, and a word's count
    // lowered for each sequence that masking takes it from.
    let mut pos_words = kmer_presence(&work, w, rc);
    // The numbers of sequences stay the same, so the P-value of each count
    // of a word in them does too.
    let mut memo = HashMap::new();
    let mut out = Vec::new();
    let mut tested = 0;
    // Each search counts once; a width that stops early counts the searches
    // it skips, so that the counts reach the total.
    let mut searched = 0;
    for _ in 0..nmotifs {
        if progress::cancelled(progress) {
            break;
        }
        searched += 1;
        progress::add(progress, 1);
        let (seeds, candidates) = enumerate_seeds(
            &pos_words,
            work.len() as u64,
            &neg_words,
            neg.len() as u64,
            lf,
            &mut memo,
        );
        tested += candidates;
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
            let eval = evaluate(&work, neg, &m, rc, lf);
            if eval.ln_pvalue > stop_pval.ln() {
                continue;
            }
            // Masking only removes words, so those a sequence loses are the
            // ones it held before and not after.
            let lost: Vec<u64> = work
                .par_iter_mut()
                .enumerate()
                .flat_map_iter(|(si, s)| {
                    if !(0..s.len()).any(|p| eval.covered.contains(si, p)) {
                        return Vec::new();
                    }
                    let room = (s.len() + 1).saturating_sub(w);
                    let mut before = vec![0; room];
                    let n = distinct_words(s, w, rc, &mut before);
                    before.truncate(n);
                    for (p, b) in s.iter_mut().enumerate() {
                        if eval.covered.contains(si, p) {
                            *b = 4;
                        }
                    }
                    let mut after = vec![0; room];
                    let n = distinct_words(s, w, rc, &mut after);
                    after.truncate(n);
                    before.retain(|x| after.binary_search(x).is_err());
                    before
                })
                .collect();
            for word in lost {
                let i = pos_words.0.binary_search(&word).unwrap();
                pos_words.1[i] -= 1;
            }
            let order = out.len();
            out.push(Found {
                motif: m,
                eval,
                width: w,
                order,
                candidates,
            });
            accepted = true;
            break;
        }
        if !accepted {
            break;
        }
    }
    progress::add(progress, (nmotifs - searched) as u64);
    (out, tested)
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
    /// Narrowest motif to search for; 6 by default.
    pub min_width: usize,
    /// Widest motif to search for, at most 30; 15 by default.
    pub max_width: usize,
    /// Most motifs searched for at each width; 10 by default.
    pub nmotifs: usize,
    /// P-value for a site to count as a hit; 1e-4 by default.
    pub hit_pvalue: f64,
    /// A motif is kept if its hits are enriched with at most this P-value;
    /// 1e-3 by default.
    pub stop_pvalue: f64,
    /// Motifs with q-values above this are dropped; 0.05 by default.
    pub qvalue: f64,
    /// Fraction of shared site positions above which the less significant
    /// of two motifs is dropped; 0.5 by default.
    pub dedup_overlap: f64,
    /// Search both strands, as by default. Motifs found without it are for
    /// the forward strand only, which [`crate::scan`] then respects.
    pub rc: bool,
    /// The size of the k-mers kept when shuffling the sequences into
    /// negatives, when none are given; 2 by default.
    pub shuffle_k: usize,
    /// How many shuffles of the sequences make up the negatives, when none
    /// are given; 3 by default. More shuffles steady the background, so
    /// that motifs in few of the sequences are found whatever the seed.
    pub bkg_shuffles: usize,
    /// Seed for shuffling the sequences into negatives.
    pub seed: u64,
    /// Pseudocount added to the motifs' counts when scoring sites, which the
    /// motifs found keep as their pseudocount.
    pub pseudocount: u32,
    /// Threads to use, one width at a time; 0 means all available.
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
            qvalue: 0.05,
            dedup_overlap: 0.5,
            rc: true,
            shuffle_k: 2,
            bkg_shuffles: 3,
            seed: 1,
            pseudocount: 1,
            nthreads: 1,
        }
    }
}

/// A motif found by [`motif_finder`], with its discovery statistics.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct FoundMotif {
    /// The motif, a PPM named `motif_1`, `motif_2` and so on, with its q-value
    /// as its E-value.
    pub motif: Motif,
    /// 1-based rank by P-value.
    pub rank: usize,
    /// Positive sequences with at least one hit.
    pub seqs_pos: u64,
    /// Negative sequences with at least one hit.
    pub seqs_neg: u64,
    /// Hits in positive sequences.
    pub sites_pos: u64,
    /// Hits in negative sequences.
    pub sites_neg: u64,
    /// The number of positive sequences.
    pub n_pos: usize,
    /// The number of negative sequences.
    pub n_neg: usize,
    /// Words tested in the search that found the motif, which its P-value is
    /// corrected for.
    pub candidates: u64,
    /// Fisher's exact test P-value of the motif's hits in positive against
    /// negative sequences, multiplied by `candidates` and capped at 1.
    pub pvalue: f64,
    /// Base 10 log of `pvalue`, which stays finite when `pvalue` is too small
    /// for a double.
    pub log10_pvalue: f64,
    /// Benjamini-Hochberg q-value, counting every word tested at every width
    /// as a test.
    pub qvalue: f64,
    /// Base 10 log of `qvalue`, which stays finite when `qvalue` is too small
    /// for a double.
    pub log10_qvalue: f64,
}

/// The background used by [`motif_finder`]: A, C, G and T or U frequencies over
/// all sequences, floored at 0.001.
pub(crate) fn finder_background<S: AsRef<str>>(pos: &[S], neg: &[S]) -> [f64; 4] {
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

/// Find motifs enriched in DNA or RNA sequences relative to negative sequences,
/// which by default are the positives shuffled keeping k-mers of size `shuffle_k`. Motifs
/// come back as PPMs named `motif_1`, `motif_2` and so on, in order of P-value,
/// with the q-value stored as their E-value.
///
/// Each motif is the best of many candidate words, so its P-value is
/// corrected for the number of words tested in the search that found it, and
/// q-values count every word tested at every width. This accounts for most of
/// the search, though not for the refinement of each motif from its word, so
/// the P-values remain somewhat optimistic.
///
/// Shuffled negatives suit sequences without much structure, but genomic
/// sequences hold many loose repeats, such as GA-rich stretches, that
/// shuffling breaks up, so such repeats come back as motifs enriched over the
/// shuffles even in regions with nothing to find. For genomic sequences,
/// negatives drawn from the genome and matched for composition, as from
/// [`crate::background::match_bkg`], give far fewer false motifs.
///
/// # References
///
/// - Bailey TL (2021). STREME: accurate and versatile sequence motif discovery.
///   *Bioinformatics* 37(18):2834–2840.
///   <https://doi.org/10.1093/bioinformatics/btab203>
/// - Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
///   practical and powerful approach to multiple testing. *Journal of the Royal
///   Statistical Society: Series B* 57(1):289–300.
///   <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
pub fn motif_finder<S: AsRef<str> + Sync>(
    sequences: &[S],
    negatives: Option<&[S]>,
    opts: &FinderOptions,
) -> Result<Vec<FoundMotif>> {
    let progress = progress::current();
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
            shuffled = shuffled_backgrounds(
                sequences,
                &ShuffleOptions {
                    k: opts.shuffle_k,
                    method: ShuffleMethod::Euler,
                    seed: opts.seed,
                    window: None,
                    nthreads: opts.nthreads,
                },
                opts.bkg_shuffles,
            )?;
            shuffled.iter().map(|s| s.as_str()).collect()
        }
    };
    if neg_strs.is_empty() {
        return arg_err("no negative sequences");
    }
    let pos_strs: Vec<&str> = sequences.iter().map(|s| s.as_ref()).collect();
    let bkg = finder_background(&pos_strs, &neg_strs);
    // The sequences as letter indices, 0 to 3 for A, C, G and T or U and 4
    // for anything else, as the search steps read them.
    let encode = |s: &&str| s.bytes().map(|b| CHAR2IDX[b as usize]).collect();
    let pos: Vec<Vec<u8>> = pos_strs.iter().map(encode).collect();
    let neg: Vec<Vec<u8>> = neg_strs.iter().map(encode).collect();

    let lf = LogFact::new(pos.len() + neg.len());
    let widths: Vec<usize> = (opts.min_width..=opts.max_width).collect();
    progress::begin(&progress, "searching", (widths.len() * opts.nmotifs) as u64);
    let per_width: Vec<(Vec<Found>, u64)> = with_threads(opts.nthreads, || {
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
                    &lf,
                    &progress,
                )
            })
            .collect()
    })?;
    progress::check(&progress)?;
    let total_tests: u64 = per_width.iter().map(|(_, t)| t).sum();
    let mut found: Vec<Found> = per_width.into_iter().flat_map(|(f, _)| f).collect();
    found.sort_by_key(|f| (f.width, f.order));
    if found.is_empty() {
        return Ok(Vec::new());
    }
    // Re-evaluate on the unmasked sequences so coverage is comparable.
    progress::begin(&progress, "evaluating", found.len() as u64);
    for f in found.iter_mut() {
        progress::check(&progress)?;
        f.eval = evaluate(&pos, &neg, &f.motif, opts.rc, &lf);
        progress::add(&progress, 1);
    }
    let corrected = |f: &Found| (f.eval.ln_pvalue + (f.candidates as f64).ln()).min(0.0);
    // Of two motifs covering the same sites, the one more enriched keeps them.
    found.sort_by(|a, b| a.eval.ln_pvalue.total_cmp(&b.eval.ln_pvalue));
    let mut dropped = vec![false; found.len()];
    let counts: Vec<u64> = found.iter().map(|f| f.eval.covered.count()).collect();
    for i in 0..found.len() {
        if dropped[i] {
            continue;
        }
        for j in i + 1..found.len() {
            if dropped[j] {
                continue;
            }
            let inter = found[i].eval.covered.shared(&found[j].eval.covered);
            let smin = counts[i].min(counts[j]);
            if smin > 0 && inter as f64 / smin as f64 > opts.dedup_overlap {
                dropped[j] = true;
            }
        }
    }
    let mut survivors: Vec<Found> = found
        .into_iter()
        .zip(dropped)
        .filter(|(_, d)| !d)
        .map(|(f, _)| f)
        .collect();
    survivors.sort_by(|a, b| corrected(a).total_cmp(&corrected(b)));
    let ln_q = crate::stats::ln_p_adjust_bh(
        &survivors
            .iter()
            .map(|f| f.eval.ln_pvalue)
            .collect::<Vec<_>>(),
        total_tests.max(survivors.len() as u64) as f64,
    );

    let mut out = Vec::new();
    for (f, ln_qv) in survivors.iter().zip(ln_q) {
        let qv = ln_qv.exp();
        let ln_pv = corrected(f);
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
            name: format!("motif_{}", out.len() + 1),
            nsites: Some(f.motif.nsites as f64),
            pseudocount: opts.pseudocount as f64,
            bkg: Some(bkg.to_vec()),
            pval: Some(ln_pv.exp()),
            eval: Some(qv),
            // A motif found on one strand is only scanned on that strand.
            strand: if opts.rc { Strand::Both } else { Strand::Plus },
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
            candidates: f.candidates,
            pvalue: ln_pv.exp(),
            log10_pvalue: ln_pv / std::f64::consts::LN_10,
            qvalue: qv,
            log10_qvalue: ln_qv / std::f64::consts::LN_10,
        });
    }
    Ok(out)
}

/// A pair of k-mers at a range of spacings, with the number of positive and
/// negative sequences holding it and its enrichment.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SpacedPair {
    /// The first k-mer, two bits per letter in ACGT order, first letter
    /// highest.
    pub left: u64,
    /// The second k-mer.
    pub right: u64,
    /// The fewest letters between them.
    pub gap: usize,
    /// The number of spacings, from `gap` up, that the pair may have.
    pub width: usize,
    /// Positive sequences holding the pair at one of these spacings.
    pub positives: u32,
    /// Negative sequences holding it.
    pub negatives: u32,
    /// The natural log of the one-sided Fisher's exact test P-value of the
    /// pair being more common in the positive sequences.
    pub ln_pvalue: f64,
}

/// Count every pair of `k`-mers separated by `min_gap` to `max_gap` letters
/// in each set of sequences, once per sequence, both at each spacing and
/// within windows of 2, 4, 8 and 16 neighbouring spacings, so that a pair
/// whose spacing varies is counted as a whole, and test the pairs found in
/// at least `min_count` positive sequences for enrichment with Fisher's
/// exact test, most enriched first. Windows wider than 2 start every half
/// width. With `rc`, a pair and its reverse complement, the reverse
/// complement of the second k-mer then of the first at the same spacing,
/// count as one, given as the smaller of the two. Windows with letters
/// other than A, C, G and T are skipped.
///
/// With `mismatch`, each pair also counts for every pair one letter away
/// from it, so that a pair is counted in each sequence holding it with at
/// most one mismatch.
///
/// Pairs are first counted over all spacings together, in the positive
/// sequences only, and only those found in at least `min_count` of them are
/// then counted spacing by spacing, which keeps memory to the pairs that can
/// pass. If more than 65,536 pairs are that common, only the 65,536 most
/// enriched over all spacings together are counted spacing by spacing.
#[allow(clippy::too_many_arguments)]
pub fn spaced_pairs<S: AsRef<str> + Sync>(
    positives: &[S],
    negatives: &[S],
    k: usize,
    min_gap: usize,
    max_gap: usize,
    rc: bool,
    min_count: u32,
    mismatch: bool,
    nthreads: usize,
) -> Result<Vec<SpacedPair>> {
    spaced_pairs_with(
        positives, negatives, k, min_gap, max_gap, rc, min_count, mismatch, nthreads, None,
    )
}

/// The window counts of the pairs `kept`, taken from those of `before`, or
/// `None` if a pair kept is not among them. Both lists are sorted.
fn subset_counts(
    kept: &[usize],
    before: &[usize],
    counts: &[u32],
    n_windows: usize,
) -> Option<Vec<u32>> {
    if before.is_empty() || kept.len() > before.len() {
        return None;
    }
    let mut out = Vec::with_capacity(kept.len() * n_windows);
    let mut j = 0;
    for &p in kept {
        while j < before.len() && before[j] < p {
            j += 1;
        }
        if j == before.len() || before[j] != p {
            return None;
        }
        out.extend_from_slice(&counts[j * n_windows..(j + 1) * n_windows]);
    }
    Some(out)
}

/// The negatives' counts from an earlier call of [`spaced_pairs_with`] with
/// the same negatives and options, which a later call can reuse: the pairs
/// counted window by window, sorted, with their counts, and the counts of
/// every pair over all spacings when they were needed. The P-values of the
/// tables tested are kept as well, for calls with as many positives.
#[derive(Default)]
pub(crate) struct NegativePairs {
    kept: Vec<usize>,
    by_window: Vec<u32>,
    /// Counts over all spacings, of the pairs in the bit set alone.
    any: Option<(Vec<u64>, Vec<u32>)>,
    sizes: (u64, u64),
    ln_p: HashMap<(u64, u64), f64>,
}

/// Presence counts over a set of slots, each slot counted at most once per
/// sequence through a bit set cleared after each one.
struct Counter {
    counts: Vec<u32>,
    seen: Vec<u64>,
    touched: Vec<usize>,
}

impl Counter {
    #[inline]
    fn add(&mut self, p: usize) {
        let (word, bit) = (p / 64, 1u64 << (p % 64));
        if self.seen[word] & bit == 0 {
            self.seen[word] |= bit;
            self.touched.push(p);
            self.counts[p] += 1;
        }
    }

    fn next_sequence(&mut self) {
        for &p in &self.touched {
            self.seen[p / 64] = 0;
        }
        self.touched.clear();
    }
}

/// Count, over `len` slots, the sequences holding each pair of k-mers
/// separated by min_gap to max_gap letters, and with `mismatch` each pair
/// one letter away from one held, both given as one 4k-bit number in the
/// reading `rc` picks. `slots` turns a pair and the index of its spacing
/// into the slots it counts in. Each pair is handled without calls through
/// closures, which matters for the 37 pairs a mismatch makes of each.
#[allow(clippy::too_many_arguments)]
fn count_pairs<S, F>(
    seqs: &[S],
    k: usize,
    min_gap: usize,
    max_gap: usize,
    rc: bool,
    mismatch: bool,
    threads: usize,
    len: usize,
    slots: &F,
) -> Result<Vec<u32>>
where
    S: AsRef<str> + Sync,
    F: Fn(usize, usize, &mut Counter) + Sync,
{
    // The shifts that change one letter of a pair and the same letter of
    // its reverse complement, which sits at the other end.
    let shifts: Vec<(u32, u32)> = (0..2 * k)
        .map(|i| (2 * i as u32, 2 * (2 * k - 1 - i) as u32))
        .collect();
    let chunk = seqs.len().div_ceil(threads).max(1);
    let parts: Vec<Vec<u32>> = with_threads(threads, || {
        seqs.par_chunks(chunk)
            .map(|part| {
                let mut c = Counter {
                    counts: vec![0u32; len],
                    seen: vec![0u64; len.div_ceil(64)],
                    touched: Vec::new(),
                };
                let mut codes: Vec<Option<(u64, u64)>> = Vec::new();
                for seq in part {
                    let bytes = seq.as_ref().as_bytes();
                    if bytes.len() < 2 * k + min_gap {
                        continue;
                    }
                    // Each k-mer with its reverse complement.
                    codes.clear();
                    codes.extend(
                        (0..=bytes.len() - k)
                            .map(|i| kmer_encode(bytes, k, i).map(|c| (c, rc_kmer(c, k)))),
                    );
                    for i in 0..codes.len() {
                        let Some((a, ra)) = codes[i] else { continue };
                        for g in min_gap..=max_gap {
                            let j = i + k + g;
                            if j >= codes.len() {
                                break;
                            }
                            let Some((b, rb)) = codes[j] else { continue };
                            // The pair as one 2k-mer, whose reverse complement
                            // is the reverse complement of the second k-mer
                            // then of the first.
                            let pair = (a << (2 * k)) | b;
                            let rev = (rb << (2 * k)) | ra;
                            let gi = g - min_gap;
                            let canon = |p: u64, r: u64| (if rc { p.min(r) } else { p }) as usize;
                            slots(canon(pair, rev), gi, &mut c);
                            if mismatch {
                                for &(sh, rsh) in &shifts {
                                    for x in 1..4u64 {
                                        slots(
                                            canon(pair ^ (x << sh), rev ^ (x << rsh)),
                                            gi,
                                            &mut c,
                                        );
                                    }
                                }
                            }
                        }
                    }
                    c.next_sequence();
                }
                c.counts
            })
            .collect()
    })?;
    let mut total = vec![0u32; len];
    for part in parts {
        for (t, v) in total.iter_mut().zip(part) {
            *t += v;
        }
    }
    Ok(total)
}

/// [`spaced_pairs`], reusing the negatives' counts in `cache` when the pairs
/// to count are among those counted before, and leaving this call's counts
/// there otherwise. Positives that only lose sites between calls, as when
/// the sites of motifs found are masked, keep fewer pairs at each call, so
/// the negatives are counted once.
#[allow(clippy::too_many_arguments)]
pub(crate) fn spaced_pairs_with<S: AsRef<str> + Sync>(
    positives: &[S],
    negatives: &[S],
    k: usize,
    min_gap: usize,
    max_gap: usize,
    rc: bool,
    min_count: u32,
    mismatch: bool,
    nthreads: usize,
    mut cache: Option<&mut NegativePairs>,
) -> Result<Vec<SpacedPair>> {
    if !(2..=6).contains(&k) {
        return arg_err("k must be from 2 to 6");
    }
    if min_gap > max_gap || max_gap > 100 {
        return arg_err("gaps must satisfy min_gap <= max_gap <= 100");
    }
    let min_count = min_count.max(1);
    let n_gaps = max_gap - min_gap + 1;
    // Windows of spacings, as (first spacing, width), each spacing's
    // windows listed with it.
    let mut windows: Vec<(usize, usize)> = Vec::new();
    for w in [1usize, 2, 4, 8, 16] {
        if w > n_gaps {
            break;
        }
        let step = if w <= 2 { 1 } else { w / 2 };
        let mut st = 0;
        loop {
            windows.push((st, w));
            if st + w >= n_gaps {
                break;
            }
            st = (st + step).min(n_gaps - w);
        }
    }
    let n_windows = windows.len();
    let of_gap: Vec<Vec<u32>> = (0..n_gaps)
        .map(|g| {
            windows
                .iter()
                .enumerate()
                .filter(|(_, (st, w))| *st <= g && g < st + w)
                .map(|(i, _)| i as u32)
                .collect()
        })
        .collect();
    let n_kmers = 1usize << (2 * k);
    let n_pairs = n_kmers * n_kmers;
    let threads = nthreads.max(1);
    let (n_pos, n_neg) = (positives.len() as u64, negatives.len() as u64);
    let lf = LogFact::new((n_pos + n_neg) as usize);
    // Many pairs share the same counts, so each table is tested once.
    let mut memo = match cache.as_deref_mut() {
        Some(c) if c.sizes == (n_pos, n_neg) => std::mem::take(&mut c.ln_p),
        _ => HashMap::new(),
    };
    let mut ln_p = |c: u64, d: u64| {
        *memo
            .entry((c, d))
            .or_insert_with(|| ln_fisher_greater(&lf, c, n_pos - c, d, n_neg - d))
    };
    // Pairs at any spacing, in the positives.
    let mut index = count_pairs(
        positives,
        k,
        min_gap,
        max_gap,
        rc,
        mismatch,
        threads,
        n_pairs,
        &|pair, _, c: &mut Counter| c.add(pair),
    )?;
    let n_kept = index.iter().filter(|&&c| c >= min_count).count();
    let mut kept: Vec<usize> = Vec::with_capacity(n_kept);
    kept.extend((0..n_pairs).filter(|&p| index[p] >= min_count));
    if kept.len() > MAX_SPACED_PAIRS {
        // Too many to count spacing by spacing: keep those most enriched
        // over all spacings together. Only these pairs are counted in the
        // negatives, which a bit per pair rules in or out cheaply.
        let mut wanted = vec![0u64; n_pairs.div_ceil(64)];
        for &p in &kept {
            wanted[p / 64] |= 1 << (p % 64);
        }
        let count_wanted = |bits: Vec<u64>| -> Result<(Vec<u64>, Vec<u32>)> {
            let counts = count_pairs(
                negatives,
                k,
                min_gap,
                max_gap,
                rc,
                mismatch,
                threads,
                n_pairs,
                &|pair, _, c: &mut Counter| {
                    if bits[pair / 64] & (1 << (pair % 64)) != 0 {
                        c.add(pair)
                    }
                },
            )?;
            Ok((bits, counts))
        };
        let fresh;
        let any_neg: &Vec<u32> = match cache.as_deref_mut() {
            Some(c) => {
                // Counts from an earlier call serve if they cover every
                // pair wanted now.
                let covered = c.any.as_ref().is_some_and(|(bits, _)| {
                    kept.iter().all(|&p| bits[p / 64] & (1 << (p % 64)) != 0)
                });
                if !covered {
                    c.any = Some(count_wanted(wanted)?);
                }
                &c.any.as_ref().unwrap().1
            }
            None => {
                fresh = count_wanted(wanted)?.1;
                &fresh
            }
        };
        let mut scored: Vec<(f64, usize)> = Vec::with_capacity(kept.len());
        scored.extend(kept.iter().map(|&p| {
            let (c, d) = (index[p] as u64, any_neg[p] as u64);
            (ln_p(c, d), p)
        }));
        scored.select_nth_unstable_by(MAX_SPACED_PAIRS - 1, |a, b| a.0.total_cmp(&b.0));
        kept.clear();
        kept.extend(scored[..MAX_SPACED_PAIRS].iter().map(|&(_, p)| p));
        kept.sort_unstable();
    }
    index.fill(u32::MAX);
    // A bit per pair, small enough to stay in cache, rules out most pairs
    // before the index is read.
    let mut is_kept = vec![0u64; n_pairs.div_ceil(64)];
    for (i, &p) in kept.iter().enumerate() {
        index[p] = i as u32;
        is_kept[p / 64] |= 1 << (p % 64);
    }
    // The pairs kept, window by window.
    let by_window = |pair: usize, gi: usize, c: &mut Counter| {
        if is_kept[pair / 64] & (1 << (pair % 64)) == 0 {
            return;
        }
        let i = index[pair];
        if i != u32::MAX {
            let base = i as usize * n_windows;
            for &win in &of_gap[gi] {
                c.add(base + win as usize);
            }
        }
    };
    let len = kept.len() * n_windows;
    let pos = count_pairs(
        positives, k, min_gap, max_gap, rc, mismatch, threads, len, &by_window,
    )?;
    // A pair's count in each window doesn't depend on the other pairs kept,
    // so counts from an earlier call serve for any of its pairs.
    let reused = cache
        .as_deref()
        .and_then(|c| subset_counts(&kept, &c.kept, &c.by_window, n_windows));
    let neg = match reused {
        Some(neg) => neg,
        None => {
            let neg = count_pairs(
                negatives, k, min_gap, max_gap, rc, mismatch, threads, len, &by_window,
            )?;
            if let Some(c) = cache.as_deref_mut() {
                c.kept = kept.clone();
                c.by_window = neg.clone();
            }
            neg
        }
    };
    // Sized up front: growing a vector this large by reallocation leaves
    // memory behind with some allocators.
    let mut out: Vec<SpacedPair> =
        Vec::with_capacity(pos.iter().filter(|&&c| c >= min_count).count());
    out.extend(
        pos.iter()
            .enumerate()
            .filter(|&(_, &c)| c >= min_count)
            .map(|(p, &c)| {
                let (st, w) = windows[p % n_windows];
                let pair = kept[p / n_windows];
                let d = neg[p] as u64;
                SpacedPair {
                    left: (pair / n_kmers) as u64,
                    right: (pair % n_kmers) as u64,
                    gap: min_gap + st,
                    width: w,
                    positives: c,
                    negatives: d as u32,
                    ln_pvalue: ln_p(c as u64, d),
                }
            }),
    );
    out.sort_by(|a, b| a.ln_pvalue.total_cmp(&b.ln_pvalue));
    if let Some(c) = cache {
        c.sizes = (n_pos, n_neg);
        c.ln_p = memo;
    }
    Ok(out)
}
