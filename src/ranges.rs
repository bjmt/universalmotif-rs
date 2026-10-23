//! Genomic ranges and the few operations on them that motif proximity testing
//! needs: reducing, intersecting, overlaps, flanks and distances.
//!
//! Ranges are 1-based and inclusive, on chromosomes identified by index. A
//! range with `end < start` is empty.

use crate::error::{Result, arg_err};

/// Strand of a genomic range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GStrand {
    /// The forward strand, written `+`.
    Plus,
    /// The reverse strand, written `-`.
    Minus,
    /// No strand, written `*`.
    Any,
}

impl GStrand {
    /// Parse `"+"` or `"-"`; `"*"`, `"+-"` and an empty string mean no strand.
    pub fn parse(s: &str) -> Result<GStrand> {
        match s {
            "+" => Ok(GStrand::Plus),
            "-" => Ok(GStrand::Minus),
            "*" | "" | "+-" => Ok(GStrand::Any),
            _ => arg_err(format!("unknown strand '{s}'; expected '+', '-' or '*'")),
        }
    }
}

/// A genomic range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GRange {
    /// The chromosome, by index.
    pub chrom: usize,
    /// The first position, counting from 1.
    pub start: i64,
    /// The last position, counting from 1 and inclusive.
    pub end: i64,
    /// The strand; no strand by default.
    pub strand: GStrand,
}

impl GRange {
    /// A range on no strand.
    pub fn new(chrom: usize, start: i64, end: i64) -> GRange {
        GRange {
            chrom,
            start,
            end,
            strand: GStrand::Any,
        }
    }

    /// The number of positions in the range, or 0 if it is empty.
    pub fn width(&self) -> i64 {
        (self.end - self.start + 1).max(0)
    }

    /// Whether the range is empty, which it is when `end < start`.
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }

    fn overlaps(&self, other: &GRange) -> bool {
        self.chrom == other.chrom && self.start <= other.end && other.start <= self.end
    }
}

/// Merge ranges whose gap is smaller than `min_gap`, ignoring strand: with
/// `min_gap = 1` adjacent ranges merge, with 0 only overlapping ones do. Empty
/// ranges are dropped; the result is sorted by chromosome and start.
pub fn reduce(ranges: &[GRange], min_gap: i64) -> Vec<GRange> {
    let mut v: Vec<GRange> = ranges.iter().filter(|r| !r.is_empty()).copied().collect();
    v.sort_by_key(|r| (r.chrom, r.start, r.end));
    let mut out: Vec<GRange> = Vec::with_capacity(v.len());
    for r in v {
        match out.last_mut() {
            Some(last) if last.chrom == r.chrom && r.start - last.end - 1 < min_gap => {
                last.end = last.end.max(r.end);
            }
            _ => out.push(GRange::new(r.chrom, r.start, r.end)),
        }
    }
    out
}

/// Intersection of two reduced range sets, as made by [`reduce`].
pub fn intersect(a: &[GRange], b: &[GRange]) -> Vec<GRange> {
    let (mut i, mut j) = (0, 0);
    let mut out = Vec::new();
    while i < a.len() && j < b.len() {
        let (x, y) = (&a[i], &b[j]);
        if x.chrom != y.chrom {
            if x.chrom < y.chrom {
                i += 1;
            } else {
                j += 1;
            }
            continue;
        }
        let (s, e) = (x.start.max(y.start), x.end.min(y.end));
        if s <= e {
            out.push(GRange::new(x.chrom, s, e));
        }
        if x.end < y.end {
            i += 1;
        } else {
            j += 1;
        }
    }
    out
}

/// Total width of non-overlapping ranges.
pub fn total_width(ranges: &[GRange]) -> f64 {
    ranges.iter().map(|r| r.width() as f64).sum()
}

/// Whether `q` overlaps any range in a reduced set.
pub fn overlaps_any(q: &GRange, reduced: &[GRange]) -> bool {
    if q.is_empty() {
        return false;
    }
    // First range on q's chromosome ending at or after q.start.
    let k = reduced.partition_point(|r| (r.chrom, r.end) < (q.chrom, q.start));
    reduced.get(k).is_some_and(|r| r.overlaps(q))
}

/// How many of `queries` overlap the reduced ranges, counted in one pass but
/// the same as counting with [`overlaps_any`]. The queries must be sorted by
/// chromosome and start.
pub(crate) fn count_overlapping(queries: &[GRange], reduced: &[GRange]) -> usize {
    let mut k = 0;
    let mut n = 0;
    for q in queries.iter().filter(|q| !q.is_empty()) {
        while reduced
            .get(k)
            .is_some_and(|r| (r.chrom, r.end) < (q.chrom, q.start))
        {
            k += 1;
        }
        if reduced.get(k).is_some_and(|r| r.overlaps(q)) {
            n += 1;
        }
    }
    n
}

/// Which side of an anchor [`zone`] covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Both sides and the anchor itself.
    Undirected,
    /// The `d` positions before the anchor on its strand.
    Upstream,
    /// The `d` positions after the anchor on its strand.
    Downstream,
}

impl Orientation {
    /// Parse `"undirected"`, `"upstream"` or `"downstream"`.
    pub fn parse(s: &str) -> Result<Orientation> {
        match s {
            "undirected" => Ok(Orientation::Undirected),
            "upstream" => Ok(Orientation::Upstream),
            "downstream" => Ok(Orientation::Downstream),
            _ => arg_err(format!(
                "unknown orientation '{s}'; expected 'undirected', 'upstream' or 'downstream'"
            )),
        }
    }
}

/// The positions within `d` of an anchor, on both sides of it or, for
/// directional zones, on one side only, clipped to the chromosome when its
/// length is known. Unstranded anchors count as `+`.
pub fn zone(anchor: &GRange, d: i64, orientation: Orientation, chrom_len: Option<i64>) -> GRange {
    let minus = anchor.strand == GStrand::Minus;
    let (start, end) = match orientation {
        Orientation::Undirected => (anchor.start - d, anchor.end + d),
        Orientation::Upstream if !minus => (anchor.start - d, anchor.start - 1),
        Orientation::Upstream => (anchor.end + 1, anchor.end + d),
        Orientation::Downstream if !minus => (anchor.end + 1, anchor.end + d),
        Orientation::Downstream => (anchor.start - d, anchor.start - 1),
    };
    let mut z = GRange::new(anchor.chrom, start, end);
    if let Some(len) = chrom_len {
        // Clip to [1, len]; a range wholly outside becomes empty.
        if z.end < 1 {
            z.start = 1;
            z.end = 0;
        } else if z.start > len {
            z.start = len + 1;
            z.end = len;
        } else {
            z.start = z.start.max(1);
            z.end = z.end.min(len);
        }
    }
    z
}

/// The gap between two ranges on the same chromosome: 0 when they overlap
/// or touch.
pub fn distance(a: &GRange, b: &GRange) -> i64 {
    (b.start - a.end - 1).max(a.start - b.end - 1).max(0)
}

/// Anchors sorted for nearest-anchor lookups.
pub(crate) struct AnchorIndex {
    /// (range, original index), sorted by chromosome and start.
    sorted: Vec<(GRange, usize)>,
    /// Running maximum of ends within each chromosome.
    max_end: Vec<i64>,
}

impl AnchorIndex {
    pub(crate) fn new(anchors: &[GRange]) -> AnchorIndex {
        let mut sorted: Vec<(GRange, usize)> = anchors.iter().copied().zip(0..).collect();
        sorted.sort_by_key(|(r, i)| (r.chrom, r.start, *i));
        let mut max_end = Vec::with_capacity(sorted.len());
        for (k, (r, _)) in sorted.iter().enumerate() {
            let prev = if k > 0 && sorted[k - 1].0.chrom == r.chrom {
                max_end[k - 1]
            } else {
                i64::MIN
            };
            max_end.push(prev.max(r.end));
        }
        AnchorIndex { sorted, max_end }
    }

    /// The original index of the anchor nearest to `q` on its chromosome,
    /// with the distance; ties go to the lowest index.
    pub(crate) fn nearest(&self, q: &GRange) -> Option<(usize, i64)> {
        let lo = self.sorted.partition_point(|(a, _)| a.chrom < q.chrom);
        let hi = self.sorted.partition_point(|(a, _)| a.chrom <= q.chrom);
        if lo == hi {
            return None;
        }
        // Anchors from p on start after q; the first of them is the closest
        // on that side. Before p, the furthest-reaching end is the closest.
        let p = lo + self.sorted[lo..hi].partition_point(|(a, _)| a.start <= q.end);
        let mut best = i64::MAX;
        if p < hi {
            best = distance(q, &self.sorted[p].0);
        }
        if p > lo {
            best = best.min((q.start - self.max_end[p - 1] - 1).max(0));
        }
        // Collect every anchor at that distance and keep the lowest index.
        let mut pick: Option<usize> = None;
        let mut k = p;
        while k < hi && self.sorted[k].0.start - q.end - 1 <= best {
            if distance(q, &self.sorted[k].0) == best {
                pick = Some(pick.map_or(self.sorted[k].1, |x| x.min(self.sorted[k].1)));
            }
            k += 1;
        }
        let mut k = p;
        while k > lo && self.max_end[k - 1] >= q.start - best - 1 {
            k -= 1;
            if distance(q, &self.sorted[k].0) == best {
                pick = Some(pick.map_or(self.sorted[k].1, |x| x.min(self.sorted[k].1)));
            }
        }
        pick.map(|i| (i, best))
    }

    /// The anchor nearest to position `pos` on chromosome `chrom`, with the
    /// number of positions from `pos` to the anchor. That number is 0 inside
    /// the anchor, so a zone widened by `d` contains `pos` exactly when it is
    /// at most `d`. Ties go to the lowest index.
    pub(crate) fn nearest_position(&self, chrom: usize, pos: i64) -> Option<(usize, i64)> {
        let q = GRange::new(chrom, pos, pos);
        let (i, gap) = self.nearest(&q)?;
        if gap > 0 {
            return Some((i, gap + 1));
        }
        // Gap 0 covers anchors containing pos and anchors next to it; those
        // containing it are closer.
        let lo = self.sorted.partition_point(|(a, _)| a.chrom < chrom);
        let p = lo + self.sorted[lo..].partition_point(|(a, _)| a.chrom == chrom && a.start <= pos);
        let mut inside: Option<usize> = None;
        let mut k = p;
        while k > lo && self.max_end[k - 1] >= pos {
            k -= 1;
            let (a, idx) = &self.sorted[k];
            if a.start <= pos && pos <= a.end {
                inside = Some(inside.map_or(*idx, |x| x.min(*idx)));
            }
        }
        Some(match inside {
            Some(idx) => (idx, 0),
            None => (i, 1),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(c: usize, s: i64, e: i64) -> GRange {
        GRange::new(c, s, e)
    }

    #[test]
    fn reduce_and_intersect() {
        let x = [r(0, 7, 9), r(0, 1, 4), r(0, 3, 6), r(1, 2, 2)];
        assert_eq!(reduce(&x, 0), vec![r(0, 1, 6), r(0, 7, 9), r(1, 2, 2)]);
        assert_eq!(reduce(&x, 1), vec![r(0, 1, 9), r(1, 2, 2)]);
        let u = [r(0, 5, 20), r(1, 1, 10)];
        assert_eq!(intersect(&reduce(&x, 1), &u), vec![r(0, 5, 9), r(1, 2, 2)]);
        assert!(overlaps_any(&r(0, 9, 12), &reduce(&x, 1)));
        assert!(!overlaps_any(&r(0, 10, 12), &reduce(&x, 1)));
        assert!(!overlaps_any(&r(2, 1, 100), &reduce(&x, 1)));
    }

    #[test]
    fn counting_matches_overlaps_any() {
        let z = reduce(&[r(0, 5, 9), r(0, 20, 30), r(1, 3, 4), r(3, 10, 12)], 1);
        let mut q: Vec<GRange> = (0..4)
            .flat_map(|c| (0..40).step_by(3).map(move |s| r(c, s, s + (s % 4))))
            .collect();
        q.push(r(0, 8, 7));
        q.sort_by_key(|g| (g.chrom, g.start, g.end));
        let brute = q.iter().filter(|g| overlaps_any(g, &z)).count();
        assert_eq!(count_overlapping(&q, &z), brute);
        assert!(brute > 0);
    }

    #[test]
    fn zones_match_genomicranges() {
        let mut a = r(0, 5, 5);
        assert_eq!(zone(&a, 3, Orientation::Upstream, None), r(0, 2, 4));
        a.strand = GStrand::Minus;
        assert_eq!(zone(&a, 3, Orientation::Upstream, None), r(0, 6, 8));
        assert_eq!(zone(&a, 3, Orientation::Downstream, None), r(0, 2, 4));
        let t = zone(&r(0, 2, 2), 5, Orientation::Undirected, Some(10));
        assert_eq!((t.start, t.end), (1, 7));
        let e = zone(&r(0, 1, 1), 3, Orientation::Upstream, Some(10));
        assert_eq!((e.start, e.end, e.width()), (1, 0, 0));
        assert_eq!(distance(&r(0, 1, 5), &r(0, 6, 10)), 0);
        assert_eq!(distance(&r(0, 1, 5), &r(0, 7, 10)), 1);
        assert_eq!(distance(&r(0, 1, 5), &r(0, 3, 10)), 0);
    }

    #[test]
    fn nearest_anchor() {
        let anchors = [
            r(0, 30, 30),
            r(0, 10, 10),
            r(0, 50, 80),
            r(1, 5, 5),
            r(0, 12, 40),
        ];
        let idx = AnchorIndex::new(&anchors);
        // Brute force with ties to the lowest index.
        for q in [
            r(0, 20, 20),
            r(0, 1, 1),
            r(0, 45, 46),
            r(0, 100, 101),
            r(1, 9, 9),
            r(0, 31, 33),
        ] {
            let want = anchors
                .iter()
                .enumerate()
                .filter(|(_, a)| a.chrom == q.chrom)
                .map(|(i, a)| (distance(&q, a), i))
                .min()
                .map(|(d, i)| (i, d));
            assert_eq!(idx.nearest(&q), want, "{q:?}");
        }
        assert_eq!(idx.nearest(&r(2, 1, 1)), None);
        // Positions: inside an anchor beats next to one.
        let idx = AnchorIndex::new(&[r(0, 11, 11), r(0, 5, 10)]);
        assert_eq!(idx.nearest_position(0, 10), Some((1, 0)));
        assert_eq!(idx.nearest_position(0, 12), Some((0, 1)));
        assert_eq!(idx.nearest_position(0, 2), Some((1, 3)));
    }
}
