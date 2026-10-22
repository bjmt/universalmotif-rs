//! Co-occurring motif pairs, tested by [`motif_coocc`].
//!
//! Every pair of motifs is tested for sharing sequences more often than
//! expected, with a one-sided Fisher's exact test on the 2x2 table of
//! per-sequence presence. Optionally, the distances between the two motifs'
//! hits in shared sequences are summarised as well; they do not feed the
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
    /// Start position, used only for distances.
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
    /// Also test each motif against itself.
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
pub struct CooccPair {
    pub a: usize,
    pub b: usize,
    /// Sequences with a hit from `a` but not `b`.
    pub a_only: u64,
    pub b_only: u64,
    pub both: u64,
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
}

impl MotifHits {
    fn new(mut hits: Vec<(usize, i64)>, words: usize) -> MotifHits {
        hits.sort_unstable();
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
        MotifHits {
            starts: hits.into_iter().map(|h| h.1).collect(),
            runs,
            present,
        }
    }

    fn n_sequences(&self) -> u64 {
        self.runs.len() as u64
    }
}

/// The smallest distance between start positions of the two motifs' hits in
/// each sequence they share. For a motif paired with itself, it is the distance
/// between two of its hits.
fn nearest_distances(a: &MotifHits, b: &MotifHits, same: bool) -> Vec<f64> {
    let mut out = Vec::new();
    if same {
        for &(_, lo, hi) in &a.runs {
            if let Some(d) = a.starts[lo..hi].windows(2).map(|w| w[1] - w[0]).min() {
                out.push(d as f64);
            }
        }
        return out;
    }
    let (mut i, mut j) = (0, 0);
    while i < a.runs.len() && j < b.runs.len() {
        let (sa, la, ha) = a.runs[i];
        let (sb, lb, hb) = b.runs[j];
        if sa < sb {
            i += 1;
        } else if sb < sa {
            j += 1;
        } else {
            // Merge the two sorted lists, tracking the closest pair.
            let (x, y) = (&a.starts[la..ha], &b.starts[lb..hb]);
            let (mut p, mut q, mut best) = (0, 0, i64::MAX);
            while p < x.len() && q < y.len() {
                best = best.min((x[p] - y[q]).abs());
                if x[p] < y[q] {
                    p += 1;
                } else {
                    q += 1;
                }
            }
            out.push(best as f64);
            i += 1;
            j += 1;
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
    let mut out: Vec<CooccPair> = with_threads(opts.nthreads, || {
        pairs
            .par_iter()
            .map(|&(i, j)| {
                let (mi, mj) = (&motifs[i], &motifs[j]);
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
                    let mut close = if both > 0 {
                        nearest_distances(mi, mj, i == j)
                    } else {
                        Vec::new()
                    };
                    close.retain(|&d| d <= maxd as f64);
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
