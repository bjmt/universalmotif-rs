//! Co-occurring motif pairs, tested by [`motif_coocc`].
//!
//! Every pair of motifs is tested for sharing sequences more often than
//! expected, with a one-sided Fisher's exact test on the 2x2 table of
//! per-sequence presence. Optionally, the distances between the two motifs'
//! hits in shared sequences are summarized as well; they do not feed the
//! test.

use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::pvalue::with_threads;
use crate::stats::{median, p_adjust_bh, phyper};

/// A motif hit, as far as [`motif_coocc`] is concerned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CooccHit {
    /// Index of the motif.
    pub motif: usize,
    /// Index of the sequence.
    pub sequence: usize,
    /// Start position, counting from 0, used only for distances.
    pub start: i64,
}

/// Options for [`motif_coocc`].
#[derive(Debug, Clone, PartialEq)]
pub struct CooccOptions {
    /// Report how many shared sequences have the two motifs' hits within this
    /// distance of each other, comparing start positions, and the median of the
    /// nearest distances.
    pub max_distance: Option<u64>,
    /// Pairs sharing fewer sequences than this are not tested.
    pub min_coocc: u64,
    /// Added to every cell of the 2x2 table.
    pub pseudocount: u64,
    /// Also pair each motif with itself. Such a pair is not tested: `both`
    /// counts the sequences with hits at two or more positions and `a_only`
    /// those with a hit at only one, and distances are between a motif's own
    /// hits.
    pub self_pairs: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for CooccOptions {
    fn default() -> Self {
        CooccOptions {
            max_distance: None,
            min_coocc: 1,
            pseudocount: 0,
            self_pairs: false,
            nthreads: 1,
        }
    }
}

/// The co-occurrence of two motifs.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct CooccPair {
    /// The first motif of the pair.
    pub a: usize,
    /// The second motif of the pair.
    pub b: usize,
    /// Sequences with a hit from `a` but not `b`.
    pub a_only: u64,
    /// Sequences with a hit from `b` but not `a`.
    pub b_only: u64,
    /// Sequences with hits from both motifs.
    pub both: u64,
    /// Sequences with hits from neither motif.
    pub neither: u64,
    /// Sample odds ratio of the table with pseudocounts; `None` if untested.
    pub odds_ratio: Option<f64>,
    /// One-sided Fisher's exact P-value; `None` if untested.
    pub pvalue: Option<f64>,
    /// Benjamini-Hochberg q-value across the tested pairs.
    pub qvalue: Option<f64>,
    /// Shared sequences whose nearest pair of hits is within `max_distance`,
    /// when that is given.
    pub both_clustered: Option<u64>,
    /// Median of those nearest distances.
    pub median_distance: Option<f64>,
}

/// One motif's hits: start positions sorted by sequence and start, with
/// the run of positions belonging to each sequence, and a bitset of the
/// sequences hit.
struct MotifHits {
    starts: Vec<i64>,
    /// (sequence, first, end) for each sequence hit, in sequence order.
    runs: Vec<(usize, usize, usize)>,
    present: Vec<u64>,
    /// Sequences with hits before each word of `present`, which is the index
    /// in `runs` of the first sequence in that word.
    rank: Vec<u32>,
}

impl MotifHits {
    fn new(mut hits: Vec<(usize, i64)>, words: usize) -> MotifHits {
        hits.sort_unstable();
        // Hits on both strands at the same position, as palindromes give,
        // count once.
        hits.dedup();
        let mut runs: Vec<(usize, usize, usize)> = Vec::new();
        let mut present = vec![0u64; words];
        for (k, &(seq, _)) in hits.iter().enumerate() {
            match runs.last_mut() {
                Some(r) if r.0 == seq => r.2 = k + 1,
                _ => {
                    runs.push((seq, k, k + 1));
                    present[seq / 64] |= 1 << (seq % 64);
                }
            }
        }
        let mut rank = Vec::with_capacity(words);
        let mut before = 0u32;
        for w in &present {
            rank.push(before);
            before += w.count_ones();
        }
        MotifHits {
            starts: hits.into_iter().map(|h| h.1).collect(),
            runs,
            present,
            rank,
        }
    }

    fn n_sequences(&self) -> u64 {
        self.runs.len() as u64
    }
}

/// The smallest distance between start positions of the two motifs' hits in
/// each sequence they share, keeping those of at most `max`. For a motif
/// paired with itself, it is the distance between two of its hits.
fn nearest_distances(a: &MotifHits, b: &MotifHits, same: bool, max: i64) -> Vec<f64> {
    let mut out = Vec::new();
    if same {
        for &(_, lo, hi) in &a.runs {
            if let Some(d) = a.starts[lo..hi].windows(2).map(|w| w[1] - w[0]).min()
                && d <= max
            {
                out.push(d as f64);
            }
        }
        return out;
    }
    // Only the sequences both motifs hit are visited, in order. Each one's
    // run is found by counting the sequences with hits before it.
    fn run_of(m: &MotifHits, w: usize, below: u64) -> &[i64] {
        let (_, lo, hi) = m.runs[(m.rank[w] + (m.present[w] & below).count_ones()) as usize];
        &m.starts[lo..hi]
    }
    for (w, (&pa, &pb)) in a.present.iter().zip(&b.present).enumerate() {
        let mut shared = pa & pb;
        while shared != 0 {
            let bit = shared.trailing_zeros();
            shared &= shared - 1;
            let below = (1u64 << bit) - 1;
            let (x, y) = (run_of(a, w, below), run_of(b, w, below));
            // Merge the two sorted lists, tracking the closest pair.
            let (mut p, mut q, mut best) = (0, 0, i64::MAX);
            while p < x.len() && q < y.len() && best > 0 {
                let d = x[p] - y[q];
                best = best.min(d.abs());
                let step = (d < 0) as usize;
                p += step;
                q += 1 - step;
            }
            if best <= max {
                out.push(best as f64);
            }
        }
    }
    out
}

/// Test motif pairs for co-occurrence across `n_sequences` sequences, from
/// a precomputed hit table. Pairs come
/// sorted by q-value then P-value, with untested pairs last.
pub fn motif_coocc(
    hits: &[CooccHit],
    n_motifs: usize,
    n_sequences: usize,
    opts: &CooccOptions,
) -> Result<Vec<CooccPair>> {
    if n_motifs < 2 && !opts.self_pairs {
        return arg_err("need at least 2 motifs (or self pairs)");
    }
    if let Some(h) = hits
        .iter()
        .find(|h| h.motif >= n_motifs || h.sequence >= n_sequences)
    {
        return arg_err(format!(
            "a hit refers to motif {} and sequence {}, but there are {n_motifs} motifs and \
             {n_sequences} sequences",
            h.motif, h.sequence
        ));
    }
    let words = n_sequences.div_ceil(64);
    let mut grouped: Vec<Vec<(usize, i64)>> = vec![Vec::new(); n_motifs];
    for h in hits {
        grouped[h.motif].push((h.sequence, h.start));
    }
    let motifs: Vec<MotifHits> = grouped
        .into_iter()
        .map(|g| MotifHits::new(g, words))
        .collect();

    let mut pairs = Vec::new();
    for i in 0..n_motifs {
        let from = if opts.self_pairs { i } else { i + 1 };
        for j in from..n_motifs {
            pairs.push((i, j));
        }
    }
    let pc = opts.pseudocount;
    // Distances are whole numbers, so comparing them with the cutoff as an
    // integer keeps the same ones as comparing as floats.
    let max_d = |m: u64| i64::try_from(m).unwrap_or(i64::MAX);
    let mut out: Vec<CooccPair> = with_threads(opts.nthreads, || {
        pairs
            .par_iter()
            .map(|&(i, j)| {
                let (mi, mj) = (&motifs[i], &motifs[j]);
                if i == j {
                    let multi = mi.runs.iter().filter(|r| r.2 - r.1 >= 2).count() as u64;
                    let mut p = CooccPair {
                        a: i,
                        b: j,
                        a_only: mi.n_sequences() - multi,
                        b_only: 0,
                        both: multi,
                        neither: n_sequences as u64 - mi.n_sequences(),
                        odds_ratio: None,
                        pvalue: None,
                        qvalue: None,
                        both_clustered: None,
                        median_distance: None,
                    };
                    if let Some(maxd) = opts.max_distance {
                        let close = nearest_distances(mi, mi, true, max_d(maxd));
                        p.both_clustered = Some(close.len() as u64);
                        p.median_distance = median(&close);
                    }
                    return p;
                }
                let both = mi
                    .present
                    .iter()
                    .zip(&mj.present)
                    .map(|(x, y)| (x & y).count_ones() as u64)
                    .sum::<u64>();
                let a_only = mi.n_sequences() - both;
                let b_only = mj.n_sequences() - both;
                let neither = n_sequences as u64 - a_only - b_only - both;
                let tested = both >= opts.min_coocc;
                let (a, b, c, d) = (both + pc, a_only + pc, b_only + pc, neither + pc);
                let mut p = CooccPair {
                    a: i,
                    b: j,
                    a_only,
                    b_only,
                    both,
                    neither,
                    odds_ratio: tested.then(|| (a as f64 * d as f64) / (b as f64 * c as f64)),
                    pvalue: tested.then(|| phyper(a as i64 - 1, a + b, c + d, a + c, false)),
                    qvalue: None,
                    both_clustered: None,
                    median_distance: None,
                };
                if let Some(maxd) = opts.max_distance {
                    let close = if both > 0 {
                        nearest_distances(mi, mj, false, max_d(maxd))
                    } else {
                        Vec::new()
                    };
                    p.both_clustered = Some(close.len() as u64);
                    p.median_distance = median(&close);
                }
                p
            })
            .collect()
    })?;

    let tested: Vec<usize> = (0..out.len())
        .filter(|&k| out[k].pvalue.is_some())
        .collect();
    let q = p_adjust_bh(
        &tested
            .iter()
            .map(|&k| out[k].pvalue.unwrap())
            .collect::<Vec<_>>(),
    );
    for (&k, q) in tested.iter().zip(q) {
        out[k].qvalue = Some(q);
    }
    // Stable sort: tested pairs by q-value then P-value, untested pairs last.
    out.sort_by(|x, y| match (x.qvalue, y.qvalue) {
        (Some(a), Some(b)) => a
            .total_cmp(&b)
            .then(x.pvalue.unwrap().total_cmp(&y.pvalue.unwrap())),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });
    Ok(out)
}
