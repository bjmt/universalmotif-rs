//! Removing overlapping hits (`dedup_hits()`).
//!
//! Hits are grouped by sequence, motif and strand (motif and strand can be
//! ignored). Within a group, overlapping hits form clusters; in each cluster
//! hits are taken in order of priority (ties go to the earlier hit) and kept
//! unless they overlap a hit already kept. This follows universalmotif's
//! `src/dedup_hits.cpp`.

/// A hit as seen by [`dedup`]. Coordinates are inclusive, and a higher
/// priority is better.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub sequence: usize,
    pub motif: usize,
    pub strand: i32,
    pub start: i64,
    pub end: i64,
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
            if ignore_strand { 0 } else { h.strand },
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
        let mut kept: Vec<(i64, i64)> = Vec::new();
        for i in by_priority {
            let (s, e) = (hits[i].start, hits[i].end);
            if !kept.iter().any(|&(ks, ke)| s <= ke && ks <= e) {
                keep[i] = true;
                kept.push((s, e));
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
            if h.start > max_end {
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
