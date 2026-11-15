//! Removing overlapping hits with [`dedup`].
//!
//! Hits are grouped by sequence, motif and strand, though motif and strand can
//! be ignored. Within a group, overlapping hits form clusters; in each cluster
//! hits are taken in order of priority, with ties going to the earlier hit, and
//! kept unless they overlap a hit already kept.

use std::collections::BTreeMap;

use crate::motif::Strand;

/// A hit as seen by [`dedup`]. Positions are half-open, and a higher
/// priority is better.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    /// Index of the sequence.
    pub sequence: usize,
    /// Index of the motif.
    pub motif: usize,
    /// The strand of the hit, with [`Strand::Both`] for hits on neither.
    /// Hits only overlap others on the same strand, unless strands are
    /// ignored.
    pub strand: Strand,
    /// First position, counting from 0.
    pub start: i64,
    /// Position just past the hit, so that it covers `start..end`.
    pub end: i64,
    /// How good the hit is; higher is better.
    pub priority: f64,
}

/// Which hits to keep, in input order.
pub fn dedup(hits: &[Interval], ignore_strand: bool, ignore_motif: bool) -> Vec<bool> {
    let n = hits.len();
    let mut keep = vec![false; n];
    let key = |i: usize| {
        let h = &hits[i];
        (
            h.sequence,
            if ignore_motif { 0 } else { h.motif },
            if ignore_strand { 0 } else { h.strand as u8 },
        )
    };
    let mut ord: Vec<usize> = (0..n).collect();
    ord.sort_by_key(|&i| (key(i), hits[i].start, i));

    let mut flush = |cluster: &[usize]| {
        if let [only] = cluster {
            keep[*only] = true;
            return;
        }
        let mut by_priority = cluster.to_vec();
        by_priority.sort_by(|&a, &b| {
            hits[b]
                .priority
                .total_cmp(&hits[a].priority)
                .then(a.cmp(&b))
        });
        // Kept hits never overlap each other, so their ends rise with their
        // starts, and a hit overlaps one of them exactly when the last to
        // start before it ends runs past its start.
        let mut kept: BTreeMap<i64, i64> = BTreeMap::new();
        for i in by_priority {
            let (s, e) = (hits[i].start, hits[i].end);
            if !kept.range(..e).next_back().is_some_and(|(_, &ke)| ke > s) {
                keep[i] = true;
                let end = kept.entry(s).or_insert(e);
                *end = (*end).max(e);
            }
        }
    };

    let mut g = 0;
    while g < n {
        let mut g_end = g + 1;
        while g_end < n && key(ord[g_end]) == key(ord[g]) {
            g_end += 1;
        }
        let mut lo = g;
        let mut max_end = hits[ord[g]].end;
        for k in g + 1..g_end {
            let h = &hits[ord[k]];
            if h.start >= max_end {
                flush(&ord[lo..k]);
                lo = k;
                max_end = h.end;
            } else if h.end > max_end {
                max_end = h.end;
            }
        }
        flush(&ord[lo..g_end]);
        g = g_end;
    }
    keep
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The plain check of every kept hit, which `dedup` must agree with.
    fn naive(hits: &[Interval]) -> Vec<bool> {
        let mut ord: Vec<usize> = (0..hits.len()).collect();
        ord.sort_by(|&a, &b| {
            hits[b]
                .priority
                .total_cmp(&hits[a].priority)
                .then(a.cmp(&b))
        });
        let mut keep = vec![false; hits.len()];
        let mut kept: Vec<(i64, i64)> = Vec::new();
        for i in ord {
            let (s, e) = (hits[i].start, hits[i].end);
            if !kept.iter().any(|&(ks, ke)| s < ke && ks < e) {
                keep[i] = true;
                kept.push((s, e));
            }
        }
        keep
    }

    #[test]
    fn agrees_with_checking_every_kept_hit() {
        let mut x: u64 = 7;
        let mut next = |m: u64| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) % m
        };
        for _ in 0..200 {
            let hits: Vec<Interval> = (0..60)
                .map(|_| {
                    let start = next(200) as i64;
                    Interval {
                        sequence: 0,
                        motif: 0,
                        strand: Strand::Plus,
                        start,
                        end: start + 1 + next(12) as i64,
                        priority: next(5) as f64,
                    }
                })
                .collect();
            assert_eq!(dedup(&hits, true, true), naive(&hits));
        }
    }
}
