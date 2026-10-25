//! Shuffling sequences and motifs with [`shuffle_sequences`] and
//! [`shuffle_motifs`], and generating random sequences with
//! [`create_sequences`].

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use rayon::prelude::*;

use crate::convert::convert_type;
use crate::error::{Result, arg_err};
use crate::kmer::{WindowOptions, windows};
use crate::matrix::Matrix;
use crate::motif::{Motif, MotifType};
use crate::pvalue::with_threads;

/// A generator for item `i` of a seeded batch, independent of threading.
pub(crate) fn rng_for(seed: u64, i: usize) -> StdRng {
    StdRng::seed_from_u64(
        seed ^ (i as u64)
            .wrapping_add(1)
            .wrapping_mul(0x9E37_79B9_7F4A_7C15),
    )
}

/// Draw an index with probability proportional to `weights`; `None` if they
/// are all zero.
pub(crate) fn weighted_index(rng: &mut StdRng, weights: &[f64]) -> Option<usize> {
    let total: f64 = weights.iter().sum();
    if total.is_nan() || total <= 0.0 {
        return None;
    }
    let mut r = rng.random::<f64>() * total;
    for (i, &w) in weights.iter().enumerate() {
        if r < w {
            return Some(i);
        }
        r -= w;
    }
    weights.iter().rposition(|&w| w > 0.0)
}

pub(crate) fn shuffle_in_place<T>(rng: &mut StdRng, v: &mut [T]) {
    for i in (1..v.len()).rev() {
        let j = rng.random_range(0..=i);
        v.swap(i, j);
    }
}

/// How sequences are shuffled when `k > 1`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ShuffleMethod {
    /// Keep the exact k-mer counts with Altschul and Erickson's Eulerian walk.
    #[default]
    Euler,
    /// Generate a new sequence from the sequence's own Markov model of order
    /// `k - 1`; k-mer counts are only kept on average.
    Markov,
    /// Shuffle consecutive blocks of `k` letters; any leftover letters stay at
    /// the end.
    Linear,
}

impl ShuffleMethod {
    /// The name [`str::parse`] reads, such as `"euler"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            ShuffleMethod::Euler => "euler",
            ShuffleMethod::Markov => "markov",
            ShuffleMethod::Linear => "linear",
        }
    }
}

impl std::str::FromStr for ShuffleMethod {
    type Err = crate::error::Error;

    /// Parse `"euler"`, `"markov"` or `"linear"`.
    fn from_str(s: &str) -> Result<ShuffleMethod> {
        match s {
            "euler" => Ok(ShuffleMethod::Euler),
            "markov" => Ok(ShuffleMethod::Markov),
            "linear" => Ok(ShuffleMethod::Linear),
            _ => arg_err(format!(
                "unknown shuffling method '{s}' (use euler, markov or linear)"
            )),
        }
    }
}

impl std::fmt::Display for ShuffleMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Options for [`shuffle_sequences`].
#[derive(Debug, Clone, PartialEq)]
pub struct ShuffleOptions {
    /// Size of the k-mers to preserve.
    pub k: usize,
    /// The shuffling method; the Euler walk by default.
    pub method: ShuffleMethod,
    /// Seed for the shuffle. Results depend on it but not on the number of
    /// threads.
    pub seed: u64,
    /// Shuffle within sliding windows instead of whole sequences.
    pub window: Option<WindowOptions>,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for ShuffleOptions {
    fn default() -> Self {
        ShuffleOptions {
            k: 1,
            method: ShuffleMethod::Euler,
            seed: 1,
            window: None,
            nthreads: 1,
        }
    }
}

/// A sequence as indices into its own sorted alphabet.
fn encode<T: Copy + Ord>(seq: &[T]) -> (Vec<T>, Vec<u32>) {
    let mut alph = seq.to_vec();
    alph.sort_unstable();
    alph.dedup();
    let ints = seq
        .iter()
        .map(|c| alph.binary_search(c).unwrap() as u32)
        .collect();
    (alph, ints)
}

fn kmer_counts(ints: &[u32], k: usize, n: usize) -> Vec<usize> {
    let mut counts = vec![0; n.pow(k as u32)];
    for w in ints.windows(k) {
        counts[w.iter().fold(0, |acc, &i| acc * n + i as usize)] += 1;
    }
    counts
}

fn shuffle_k1<T: Copy>(seq: &[T], rng: &mut StdRng) -> Vec<T> {
    let mut out = seq.to_vec();
    shuffle_in_place(rng, &mut out);
    out
}

fn shuffle_linear<T: Copy>(seq: &[T], k: usize, rng: &mut StdRng) -> Vec<T> {
    let nblocks = seq.len() / k;
    let mut idx: Vec<usize> = (0..nblocks).map(|i| i * k).collect();
    shuffle_in_place(rng, &mut idx);
    let mut out = Vec::with_capacity(seq.len());
    for i in idx {
        out.extend_from_slice(&seq[i..i + k]);
    }
    out.extend_from_slice(&seq[nblocks * k..]);
    out
}

fn shuffle_markov<T: Copy + Ord>(seq: &[T], k: usize, rng: &mut StdRng) -> Vec<T> {
    if seq.len() < k {
        return seq.to_vec();
    }
    let (alph, ints) = encode(seq);
    let n = alph.len();
    let counts = kmer_counts(&ints, k, n);
    let letters = kmer_counts(&ints, 1, n);
    let weights: Vec<f64> = counts.iter().map(|&c| c as f64).collect();
    let first = weighted_index(rng, &weights).unwrap_or(0);
    let mut out: Vec<u32> = (0..k)
        .map(|p| ((first / n.pow((k - 1 - p) as u32)) % n) as u32)
        .collect();
    let mlets = n.pow(k as u32 - 1);
    let fallback: Vec<f64> = letters.iter().map(|&c| c as f64).collect();
    while out.len() < seq.len() {
        let ctx = out[out.len() - (k - 1)..]
            .iter()
            .fold(0, |acc, &i| acc * n + i as usize)
            % mlets;
        let row: Vec<f64> = counts[ctx * n..(ctx + 1) * n]
            .iter()
            .map(|&c| c as f64)
            .collect();
        let next = weighted_index(rng, &row)
            .or_else(|| weighted_index(rng, &fallback))
            .unwrap_or(0);
        out.push(next as u32);
    }
    out.truncate(seq.len());
    out.into_iter().map(|i| alph[i as usize]).collect()
}

fn shuffle_euler<T: Copy + Ord>(seq: &[T], k: usize, rng: &mut StdRng) -> Vec<T> {
    let (alph, ints) = encode(seq);
    let n = alph.len();
    if n == 1 || seq.len() <= k {
        return seq.to_vec();
    }
    let mlets = n.pow(k as u32 - 1);
    let mmlets = n.pow(k as u32 - 2);
    let counts = kmer_counts(&ints, k, n);
    let mut edges: Vec<Vec<usize>> = (0..mlets)
        .map(|v| counts[v * n..(v + 1) * n].to_vec())
        .collect();
    let next_vertex = |v: usize, letter: usize| (v % mmlets) * n + letter;
    let last = ints[ints.len() - (k - 1)..]
        .iter()
        .fold(0, |acc, &i| acc * n + i as usize);

    // A random spanning tree of last exits, pointing towards the last vertex
    // (Wilson's algorithm).
    let empty: Vec<bool> = edges.iter().map(|e| e.iter().all(|&c| c == 0)).collect();
    let mut in_tree: Vec<bool> = empty.clone();
    in_tree[last] = true;
    let mut exit = vec![0usize; mlets];
    for i in 0..mlets {
        let mut u = i;
        while !in_tree[u] {
            let w: Vec<f64> = edges[u].iter().map(|&c| c as f64).collect();
            exit[u] = weighted_index(rng, &w).unwrap_or(0);
            u = next_vertex(u, exit[u]);
        }
        let mut u = i;
        while !in_tree[u] {
            in_tree[u] = true;
            u = next_vertex(u, exit[u]);
        }
    }

    // Each vertex's remaining edges in random order, with its last exit last.
    let mut lists: Vec<Vec<u32>> = vec![Vec::new(); mlets];
    for v in 0..mlets {
        if empty[v] {
            continue;
        }
        if v != last {
            edges[v][exit[v]] -= 1;
        }
        let mut l: Vec<u32> = edges[v]
            .iter()
            .enumerate()
            .flat_map(|(letter, &c)| std::iter::repeat_n(letter as u32, c))
            .collect();
        shuffle_in_place(rng, &mut l);
        if v != last {
            l.push(exit[v] as u32);
        }
        lists[v] = l;
    }

    let mut out: Vec<u32> = ints[..k - 1].to_vec();
    let mut used = vec![0usize; mlets];
    while out.len() < ints.len() {
        let v = out[out.len() - (k - 1)..]
            .iter()
            .fold(0, |acc, &i| acc * n + i as usize);
        out.push(lists[v][used[v]]);
        used[v] += 1;
    }
    out.into_iter().map(|i| alph[i as usize]).collect()
}

fn shuffle_one<T: Copy + Ord>(
    seq: &[T],
    k: usize,
    method: ShuffleMethod,
    rng: &mut StdRng,
) -> Vec<T> {
    if k == 1 || seq.len() < 2 {
        return shuffle_k1(seq, rng);
    }
    match method {
        ShuffleMethod::Euler => shuffle_euler(seq, k, rng),
        ShuffleMethod::Markov => shuffle_markov(seq, k, rng),
        ShuffleMethod::Linear => shuffle_linear(seq, k, rng),
    }
}

/// Shuffle sequences while preserving their k-mer composition, exactly with the
/// Euler and linear methods and on average with the Markov method. Each
/// sequence is shuffled with its own alphabet.
///
/// # Example
///
/// ```
/// use universalmotif::shuffle::{ShuffleOptions, shuffle_sequences};
///
/// let opts = ShuffleOptions { k: 2, seed: 1, ..Default::default() };
/// let shuffled = shuffle_sequences(&["ACGTTGCAAGGCTTAC"], &opts)?;
/// let mut letters: Vec<char> = shuffled[0].chars().collect();
/// letters.sort();
/// assert_eq!(letters.iter().collect::<String>(), "AAAACCCCGGGGTTTT");
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn shuffle_sequences<S: AsRef<str> + Sync>(
    seqs: &[S],
    opts: &ShuffleOptions,
) -> Result<Vec<String>> {
    if opts.k < 1 {
        return arg_err("k must be at least 1");
    }
    let mut plans: Vec<Vec<(usize, usize)>> = Vec::with_capacity(seqs.len());
    for s in seqs {
        let s = s.as_ref();
        let len = if s.is_ascii() {
            s.len()
        } else {
            s.chars().count()
        };
        match opts.window {
            None => {
                if len < opts.k {
                    return arg_err(format!(
                        "k ({}) must not exceed the shortest sequence length ({len})",
                        opts.k
                    ));
                }
                plans.push(vec![(0, len)]);
            }
            Some(w) => {
                if w.size.is_nan() || w.size <= 0.0 || w.overlap < 0.0 {
                    return arg_err(
                        "the window size must be positive and the overlap non-negative",
                    );
                }
                let resolve = |v: f64| {
                    if v < 1.0 {
                        (v * len as f64) as usize
                    } else {
                        v as usize
                    }
                };
                let (size, overlap) = (resolve(w.size), resolve(w.overlap));
                if size <= overlap {
                    return arg_err("the window overlap must be smaller than the window size");
                }
                if size < opts.k {
                    return arg_err("the window size must be at least k");
                }
                plans.push(windows(len, size, overlap, true));
            }
        }
    }
    // ASCII sequences (the usual case) are shuffled as bytes; anything else
    // letter by letter. The random draws are the same either way.
    fn shuffle_windows<T: Copy + Ord>(
        mut out: Vec<T>,
        plan: &[(usize, usize)],
        opts: &ShuffleOptions,
        rng: &mut StdRng,
    ) -> Vec<T> {
        for &(start, end) in plan {
            let part = shuffle_one(&out[start..end], opts.k, opts.method, rng);
            out[start..end].copy_from_slice(&part);
        }
        out
    }
    with_threads(opts.nthreads, || {
        seqs.par_iter()
            .zip(&plans)
            .enumerate()
            .map(|(i, (s, plan))| {
                let s = s.as_ref();
                let mut rng = rng_for(opts.seed, i);
                if s.is_ascii() {
                    let out = shuffle_windows(s.as_bytes().to_vec(), plan, opts, &mut rng);
                    String::from_utf8(out).expect("shuffled ASCII is valid UTF-8")
                } else {
                    let chars: Vec<char> = s.chars().collect();
                    shuffle_windows(chars, plan, opts, &mut rng)
                        .into_iter()
                        .collect()
                }
            })
            .collect()
    })
}

/// Random sequences of `alphabet` letters.
///
/// `freqs` gives the probabilities of every k-mer of the alphabet in
/// [`crate::kmer::kmers`] order, and its length decides k. Sequences are drawn
/// from the Markov model of order `k - 1` they define. By default letters are
/// uniform.
pub fn create_sequences(
    alphabet: &[char],
    n: usize,
    len: usize,
    freqs: Option<&[f64]>,
    seed: u64,
    nthreads: usize,
) -> Result<Vec<String>> {
    let a = alphabet.len();
    if a < 1 {
        return arg_err("the alphabet needs at least one letter");
    }
    let uniform = vec![1.0 / a as f64; a];
    let freqs = freqs.unwrap_or(&uniform);
    let mut k = 1usize;
    while a.pow(k as u32) < freqs.len() {
        k += 1;
    }
    if a.pow(k as u32) != freqs.len() {
        return arg_err("the number of frequencies must be a power of the alphabet size");
    }
    if freqs.iter().any(|f| !f.is_finite() || *f < 0.0) {
        return arg_err("frequencies must be finite and non-negative");
    }
    if !freqs.iter().any(|&f| f > 0.0) {
        return arg_err("at least one frequency must be above zero");
    }
    with_threads(nthreads, || {
        (0..n)
            .into_par_iter()
            .map(|i| {
                let mut rng = rng_for(seed, i);
                let mut out = String::with_capacity(len);
                if k == 1 {
                    for _ in 0..len {
                        out.push(alphabet[weighted_index(&mut rng, freqs).unwrap_or(0)]);
                    }
                } else {
                    // The first k letters, then each letter given the last
                    // k - 1, whose code is kept as it rolls along.
                    let first = weighted_index(&mut rng, freqs).unwrap_or(0);
                    for p in 0..k.min(len) {
                        out.push(alphabet[(first / a.pow((k - 1 - p) as u32)) % a]);
                    }
                    let mlets = a.pow(k as u32 - 1);
                    let mut ctx = first % mlets;
                    for _ in k..len {
                        let row = &freqs[ctx * a..(ctx + 1) * a];
                        let x = weighted_index(&mut rng, row).unwrap_or(0);
                        out.push(alphabet[x]);
                        ctx = (ctx * a + x) % mlets;
                    }
                }
                out
            })
            .collect()
    })
}

/// Shuffle blocks of `k` items, keeping any leftover items at the two ends.
fn shuffle_blocks_keep_ends(n: usize, k: usize, rng: &mut StdRng) -> Vec<usize> {
    let idx: Vec<usize> = (0..n).collect();
    let rem = n % k;
    let (left, right) = match rem {
        0 => (0, 0),
        1 => (1, 0),
        _ => {
            let l = ((rem as f64 + 0.1) / 2.0).round() as usize;
            (l, rem - l)
        }
    };
    let middle = &idx[left..n - right];
    let mut blocks: Vec<&[usize]> = middle.chunks(k).collect();
    shuffle_in_place(rng, &mut blocks);
    let mut out = idx[..left].to_vec();
    for b in blocks {
        out.extend_from_slice(b);
    }
    out.extend_from_slice(&idx[n - right..]);
    out
}

/// Shuffle the columns of motifs among themselves in blocks of `k` columns,
/// producing new motifs of the same widths and types with `" [shuffled]"` added
/// to their names.
pub fn shuffle_motifs(motifs: &[Motif], k: usize, seed: u64) -> Result<Vec<Motif>> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs given");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if k < 1 {
        return arg_err("k must be at least 1");
    }
    let ppms: Vec<Motif> = motifs
        .iter()
        .map(|m| convert_type(m, MotifType::Ppm, None, false))
        .collect::<Result<_>>()?;
    let cols: Vec<&[f64]> = ppms.iter().flat_map(|m| m.matrix.columns()).collect();
    let mut rng = rng_for(seed, 0);
    let order = if k == 1 {
        let mut o: Vec<usize> = (0..cols.len()).collect();
        shuffle_in_place(&mut rng, &mut o);
        o
    } else {
        shuffle_blocks_keep_ends(cols.len(), k, &mut rng)
    };
    let mut out = Vec::with_capacity(motifs.len());
    let mut next = 0;
    for (orig, ppm) in motifs.iter().zip(&ppms) {
        let w = ppm.width();
        let new: Vec<&[f64]> = order[next..next + w].iter().map(|&i| cols[i]).collect();
        next += w;
        let m = Motif::builder(Matrix::from_columns(&new)?, ppm.alphabet.clone())
            .name(format!("{} [shuffled]", ppm.name))
            .bkg(Some(ppm.bkg.clone()))
            .bkgsites(ppm.bkgsites)
            .nsites(ppm.nsites)
            .strand(ppm.strand)
            .pseudocount(ppm.pseudocount)
            .build()?;
        out.push(
            if orig.motif_type != MotifType::Ppm && m.motif_type != orig.motif_type {
                convert_type(&m, orig.motif_type, None, false)?
            } else {
                m
            },
        );
    }
    Ok(out)
}
