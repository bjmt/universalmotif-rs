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
///
/// # References
///
/// - Altschul SF, Erickson BW (1985). Significance of nucleotide sequence
///   alignments: a method for random sequence permutation that preserves
///   dinucleotide and codon usage. *Molecular Biology and Evolution*
///   2(6):526–538. <https://doi.org/10.1093/oxfordjournals.molbev.a040370>
/// - Jiang M, Anderson J, Gillespie J, Mayne M (2008). uShuffle: a useful tool
///   for shuffling biological sequences while preserving the k-let counts. *BMC
///   Bioinformatics* 9:192. <https://doi.org/10.1186/1471-2105-9-192>
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
    /// Size of the k-mers to preserve; 2 by default, since the dinucleotide
    /// make-up of DNA, CpG depletion above all, changes how often motifs
    /// match, so backgrounds should keep it.
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
            k: 2,
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

/// A generator for copy `j` of item `i` of a seeded batch, which does not
/// depend on how many copies are made or on threading.
fn rng_for_copy(seed: u64, i: u64, j: u64) -> StdRng {
    // SplitMix64's finalizer, applied in turn to the item and the copy.
    let mix = |mut z: u64| {
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    };
    let a = mix(seed ^ i.wrapping_add(1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    StdRng::seed_from_u64(mix(a ^ j
        .wrapping_add(1)
        .wrapping_mul(0xD6E8_FEB8_6659_FD93)))
}

/// Options for [`shuffle_codes`].
#[derive(Debug, Clone, PartialEq)]
pub struct CodeShuffleOptions {
    /// Size of the k-mers to preserve, 2 by default.
    pub k: usize,
    /// The shuffling method; the Euler walk by default.
    pub method: ShuffleMethod,
    /// Seed for the shuffles.
    pub seed: u64,
    /// The index of the first sequence among all those shuffled with this
    /// seed, so that shuffling a large set in batches gives the same
    /// shuffles as shuffling it at once.
    pub offset: u64,
    /// A code for unknown letters, such as N, which stay where they are; the
    /// letters between them are shuffled stretch by stretch.
    pub unknown: Option<u8>,
    /// Threads to use, at least 1.
    pub nthreads: usize,
}

impl Default for CodeShuffleOptions {
    fn default() -> Self {
        CodeShuffleOptions {
            k: 2,
            method: ShuffleMethod::Euler,
            seed: 1,
            offset: 0,
            unknown: None,
            nthreads: 1,
        }
    }
}

/// Several shuffles of each of a set of sequences of equal length, given as
/// letter codes, such as the references that DeepLIFT and DeepSHAP compare
/// each input with.
///
/// `codes` holds `codes.len() / length` sequences one after another, and the
/// result holds `n` shuffles of each in turn, `sequences * n * length`
/// codes. With the Euler method, each shuffle keeps the exact counts of every
/// k-mer of its sequence, along with its first `k - 1` letters and its last
/// letter, and is drawn uniformly from all the sequences that do. Positions
/// holding `opts.unknown` stay where they are, and each stretch of letters
/// between them is shuffled on its own, so that k-mers touching them are
/// kept too. Copy `j` of sequence `i` depends only on `opts.seed`,
/// `opts.offset + i` and `j`, not on `n` or the number of threads.
///
/// # Example
///
/// ```
/// use universalmotif::shuffle::{CodeShuffleOptions, shuffle_codes};
///
/// // ACGTTGCA and AACCNNGT, as codes 0 to 3, with 4 for N.
/// let codes = [0, 1, 2, 3, 3, 2, 1, 0, 0, 0, 1, 1, 4, 4, 2, 3];
/// let opts = CodeShuffleOptions { unknown: Some(4), ..Default::default() };
/// let out = shuffle_codes(&codes, 8, 3, &opts)?;
/// assert_eq!(out.len(), 2 * 3 * 8);
/// // The Ns of the second sequence stay where they are.
/// for copy in out[24..].chunks(8) {
///     assert_eq!(&copy[4..6], &[4, 4]);
/// }
/// # Ok::<(), universalmotif::Error>(())
/// ```
///
/// # References
///
/// - Shrikumar A, Greenside P, Kundaje A (2017). Learning important features
///   through propagating activation differences. *Proceedings of the 34th
///   International Conference on Machine Learning*, PMLR 70:3145–3153.
///   <https://proceedings.mlr.press/v70/shrikumar17a.html>
/// - Lundberg SM, Lee SI (2017). A unified approach to interpreting model
///   predictions. *Advances in Neural Information Processing Systems*
///   30:4765–4774. <https://arxiv.org/abs/1705.07874>
pub fn shuffle_codes(
    codes: &[u8],
    length: usize,
    n: usize,
    opts: &CodeShuffleOptions,
) -> Result<Vec<u8>> {
    if opts.k < 1 {
        return arg_err("k must be at least 1");
    }
    if length == 0 {
        return if codes.is_empty() {
            Ok(Vec::new())
        } else {
            arg_err("sequences of length 0 cannot hold any codes")
        };
    }
    if !codes.len().is_multiple_of(length) {
        return arg_err(format!(
            "{} codes do not divide into sequences of length {length}",
            codes.len()
        ));
    }
    let shuffle_seq = |seq: &[u8], rng: &mut StdRng| -> Vec<u8> {
        let mut out = seq.to_vec();
        let mut start = 0;
        while start < seq.len() {
            if opts.unknown == Some(seq[start]) {
                start += 1;
                continue;
            }
            let end = seq[start..]
                .iter()
                .position(|&c| opts.unknown == Some(c))
                .map_or(seq.len(), |p| start + p);
            let part = &seq[start..end];
            if part.len() >= opts.k.max(2) {
                out[start..end].copy_from_slice(&shuffle_one(part, opts.k, opts.method, rng));
            }
            start = end;
        }
        out
    };
    let parts: Vec<Vec<u8>> = with_threads(opts.nthreads, || {
        codes
            .par_chunks(length)
            .enumerate()
            .map(|(i, seq)| {
                let mut out = Vec::with_capacity(n * length);
                for j in 0..n {
                    let mut rng = rng_for_copy(opts.seed, opts.offset + i as u64, j as u64);
                    out.extend(shuffle_seq(seq, &mut rng));
                }
                out
            })
            .collect()
    })?;
    Ok(parts.concat())
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
            .motif_type(Some(MotifType::Ppm))
            .name(format!("{} [shuffled]", ppm.name))
            .bkg(Some(ppm.bkg.clone()))
            .bkgsites(ppm.bkgsites)
            .nsites(ppm.nsites)
            .strand(ppm.strand)
            .pseudocount(ppm.pseudocount)
            .build()?;
        out.push(if orig.motif_type != MotifType::Ppm {
            convert_type(&m, orig.motif_type, None, false)?
        } else {
            m
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(seq: &[u8], k: usize) -> std::collections::BTreeMap<Vec<u8>, usize> {
        let mut m = std::collections::BTreeMap::new();
        for w in seq.windows(k) {
            *m.entry(w.to_vec()).or_insert(0) += 1;
        }
        m
    }

    #[test]
    fn code_shuffles_keep_kmers_and_unknowns() {
        let seqs: Vec<u8> = b"ACGTTGCAAGGCTTACNNGATTACAGGATTCCAGN"
            .iter()
            .map(|c| b"ACGTN".iter().position(|x| x == c).unwrap() as u8)
            .collect();
        let length = seqs.len();
        for k in [1, 2, 3] {
            let opts = CodeShuffleOptions {
                k,
                unknown: Some(4),
                ..Default::default()
            };
            let out = shuffle_codes(&seqs, length, 20, &opts).unwrap();
            let mut differ = 0;
            for copy in out.chunks(length) {
                assert_eq!(counts(copy, k), counts(&seqs, k));
                for (a, b) in copy.iter().zip(&seqs) {
                    assert_eq!(*a == 4, *b == 4);
                }
                if k > 1 {
                    // Each stretch keeps its first k - 1 letters and its last.
                    assert_eq!(copy[..k - 1], seqs[..k - 1]);
                    assert_eq!(copy[15], seqs[15]);
                    assert_eq!(copy[18..18 + k - 1], seqs[18..18 + k - 1]);
                }
                differ += usize::from(copy != seqs.as_slice());
            }
            assert!(differ > 0);
        }
    }

    #[test]
    fn code_shuffles_are_reproducible_in_batches() {
        let seqs: Vec<u8> = (0..6 * 40).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
        let opts = CodeShuffleOptions {
            seed: 7,
            nthreads: 2,
            ..Default::default()
        };
        let all = shuffle_codes(&seqs, 40, 5, &opts).unwrap();
        let first = shuffle_codes(&seqs[..2 * 40], 40, 5, &opts).unwrap();
        let rest = shuffle_codes(
            &seqs[2 * 40..],
            40,
            5,
            &CodeShuffleOptions {
                offset: 2,
                nthreads: 1,
                ..opts.clone()
            },
        )
        .unwrap();
        assert_eq!([first, rest].concat(), all);
        // More copies extend the earlier ones.
        let more = shuffle_codes(&seqs[..40], 40, 8, &opts).unwrap();
        assert_eq!(more[..5 * 40], all[..5 * 40]);
        let other = shuffle_codes(&seqs, 40, 5, &CodeShuffleOptions { seed: 8, ..opts }).unwrap();
        assert_ne!(other, all);
    }

    #[test]
    fn code_shuffles_are_uniform() {
        // ACGTACGATCGA has 15 distinct shuffles keeping its dinucleotides.
        let seq: Vec<u8> = b"ACGTACGATCGA"
            .iter()
            .map(|c| b"ACGT".iter().position(|x| x == c).unwrap() as u8)
            .collect();
        let n = 6000;
        let out = shuffle_codes(&seq, seq.len(), n, &CodeShuffleOptions::default()).unwrap();
        let mut seen = std::collections::HashMap::new();
        for copy in out.chunks(seq.len()) {
            *seen.entry(copy.to_vec()).or_insert(0usize) += 1;
        }
        let expected = n as f64 / seen.len() as f64;
        assert_eq!(seen.len(), 15);
        for &c in seen.values() {
            assert!(
                (c as f64 - expected).abs() < 5.0 * expected.sqrt(),
                "{c} vs {expected}"
            );
        }
    }

    #[test]
    fn code_shuffle_checks() {
        let opts = CodeShuffleOptions::default();
        assert!(shuffle_codes(&[0, 1, 2], 2, 1, &opts).is_err());
        assert!(
            shuffle_codes(
                &[0, 1],
                2,
                1,
                &CodeShuffleOptions {
                    k: 0,
                    ..opts.clone()
                }
            )
            .is_err()
        );
        assert!(shuffle_codes(&[], 0, 3, &opts).unwrap().is_empty());
        assert!(
            shuffle_codes(&[0, 1, 2, 3], 4, 0, &opts)
                .unwrap()
                .is_empty()
        );
    }
}
