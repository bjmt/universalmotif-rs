//! Co-occurring motif pairs (`motif_coocc()`).
//!
//! Every pair of motifs is tested for sharing sequences more often than
//! expected, with a one-sided Fisher's exact test on the 2x2 table of
//! per-sequence presence. Optionally, the distances between the two motifs'
//! hits in shared sequences are summarised as well; they do not feed the
//! test.

use std::collections::BTreeSet;

use crate::error::{Result, arg_err};
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
    /// Report how many shared sequences have the two motifs' hits within
    /// this distance of each other (comparing start positions), and the
    /// median of the nearest distances.
    pub max_distance: Option<u64>,
    /// Pairs sharing fewer sequences than this are not tested.
    pub min_coocc: u64,
    /// Added to every cell of the 2x2 table.
    pub pseudocount: u64,
    /// Also test each motif against itself.
    pub self_pairs: bool,
}

impl Default for CooccOptions {
    fn default() -> Self {
        CooccOptions {
            max_distance: None,
            min_coocc: 1,
            pseudocount: 0,
            self_pairs: false,
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
    /// Shared sequences whose nearest pair of hits is within
    /// `max_distance` (only with `max_distance`).
    pub both_clustered: Option<u64>,
    /// Median of those nearest distances.
    pub median_distance: Option<f64>,
}

/// Nearest distance between hits of two motifs (or between two hits of one
/// motif) in each shared sequence.
fn nearest_distances(sa: &[(usize, i64)], sb: &[(usize, i64)], same: bool) -> Vec<f64> {
    // Both lists are sorted by sequence, then start.
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < sa.len() {
        let seq = sa[i].0;
        let ia = sa[i..].iter().take_while(|h| h.0 == seq).count();
        let a = &sa[i..i + ia];
        i += ia;
        if same {
            if let Some(d) = a.windows(2).map(|w| w[1].1 - w[0].1).min() {
                out.push(d as f64);
            }
            continue;
        }
        while j < sb.len() && sb[j].0 < seq {
            j += 1;
        }
        let jb = sb[j..].iter().take_while(|h| h.0 == seq).count();
        let b = &sb[j..j + jb];
        if b.is_empty() {
            continue;
        }
        let d = a
            .iter()
            .map(|&(_, x)| {
                let k = b.partition_point(|&(_, y)| y <= x);
                let below = if k > 0 { x - b[k - 1].1 } else { i64::MAX };
                let above = if k < b.len() { b[k].1 - x } else { i64::MAX };
                below.min(above)
            })
            .min()
            .unwrap();
        out.push(d as f64);
    }
    out
}

/// Test motif pairs for co-occurrence across `n_sequences` sequences, as
/// universalmotif's `motif_coocc()` with a precomputed hit table. Pairs come
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
    let mut presence: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n_motifs];
    let mut by_motif: Vec<Vec<(usize, i64)>> = vec![Vec::new(); n_motifs];
    for h in hits {
        presence[h.motif].insert(h.sequence);
        by_motif[h.motif].push((h.sequence, h.start));
    }
    by_motif.iter_mut().for_each(|v| v.sort_unstable());

    let mut pairs = Vec::new();
    for i in 0..n_motifs {
        let from = if opts.self_pairs { i } else { i + 1 };
        for j in from..n_motifs {
            pairs.push((i, j));
        }
    }
    let pc = opts.pseudocount;
    let mut out: Vec<CooccPair> = pairs
        .iter()
        .map(|&(i, j)| {
            let both = presence[i].intersection(&presence[j]).count() as u64;
            let a_only = presence[i].len() as u64 - both;
            let b_only = presence[j].len() as u64 - both;
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
                let mut close = Vec::new();
                if both > 0 {
                    let sa: Vec<(usize, i64)> = by_motif[i]
                        .iter()
                        .copied()
                        .filter(|h| presence[j].contains(&h.0))
                        .collect();
                    let sb: Vec<(usize, i64)> = by_motif[j]
                        .iter()
                        .copied()
                        .filter(|h| presence[i].contains(&h.0))
                        .collect();
                    close = nearest_distances(&sa, &sb, i == j);
                    close.retain(|&d| d <= maxd as f64);
                }
                p.both_clustered = Some(close.len() as u64);
                p.median_distance = median(&close);
            }
            p
        })
        .collect();

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
