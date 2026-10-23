//! Motif enrichment near genomic anchors, tested by [`motif_proximity`].
//!
//! The genome-wide counterpart of [`crate::peaks`]: are a motif's hits closer
//! to a set of anchor coordinates, such as summits, TSSs or SNPs, than chance
//! would put them? Each hit is placed at its middle position. For each
//! candidate distance `d`, the near zone is every scanned position within `d`
//! of an anchor; under the null that hits are uniform over the scanned
//! universe, the number in the zone is binomial. The best distance is
//! Bonferroni-corrected for the distances tried. A permutation null with random
//! anchors and a one-sided Kolmogorov-Smirnov test of the distances to the
//! nearest anchor are also available.

use std::collections::{BTreeSet, HashMap};

use rand::RngExt;
use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::pvalue::with_threads;
use crate::ranges::{
    AnchorIndex, GRange, GStrand, Orientation, count_overlapping, intersect, overlaps_any, reduce,
    total_width, zone,
};
use crate::shuffle::rng_for;
use crate::stats::{p_adjust_bh, pbinom_upper};

/// The null model of [`motif_proximity`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProximityMethod {
    /// Hits uniform over the universe.
    Binomial,
    /// Anchors placed at random in the universe.
    Permutation,
    /// One-sided Kolmogorov-Smirnov test of the distances.
    Ks,
}

impl ProximityMethod {
    /// Parse `"binomial"`, `"permutation"` or `"ks"`.
    pub fn parse(s: &str) -> Result<ProximityMethod> {
        match s {
            "binomial" => Ok(ProximityMethod::Binomial),
            "permutation" => Ok(ProximityMethod::Permutation),
            "ks" => Ok(ProximityMethod::Ks),
            _ => arg_err(format!(
                "unknown method '{s}'; expected 'binomial', 'permutation' or 'ks'"
            )),
        }
    }
}

/// What [`motif_proximity`] counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProximityCount {
    /// Hits in the near zone, where overlapping hits of a motif count once.
    Hits,
    /// Anchors with a hit nearby.
    Anchors,
}

impl ProximityCount {
    /// Parse `"hits"` or `"anchors"`.
    pub fn parse(s: &str) -> Result<ProximityCount> {
        match s {
            "hits" => Ok(ProximityCount::Hits),
            "anchors" => Ok(ProximityCount::Anchors),
            _ => arg_err(format!("unknown count '{s}'; expected 'hits' or 'anchors'")),
        }
    }
}

/// Options for [`motif_proximity`].
#[derive(Debug, Clone, PartialEq)]
pub struct ProximityOptions {
    /// The null model; binomial by default.
    pub method: ProximityMethod,
    /// Whether to count hits or anchors; hits by default.
    pub count: ProximityCount,
    /// Which side of stranded anchors the near zone covers; both sides by
    /// default.
    pub orientation: Orientation,
    /// Distances to test; by default a log-spaced grid up to half the
    /// smallest universe segment.
    pub dist_thresholds: Option<Vec<i64>>,
    /// Random placements of the anchors with the permutation null; 1000 by
    /// default.
    pub n_perm: usize,
    /// Keep permuted anchors on their own chromosome.
    pub perm_per_chrom: bool,
    /// Collapse overlapping hits regardless of strand.
    pub ignore_strand: bool,
    /// Seed for the permutation null.
    pub seed: u64,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for ProximityOptions {
    fn default() -> Self {
        ProximityOptions {
            method: ProximityMethod::Binomial,
            count: ProximityCount::Hits,
            orientation: Orientation::Undirected,
            dist_thresholds: None,
            n_perm: 1000,
            perm_per_chrom: true,
            ignore_strand: true,
            seed: 1,
            nthreads: 1,
        }
    }
}

/// A motif hit for [`motif_proximity`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProximityHit {
    /// Motif identifier, of any value; results come in increasing order of it.
    pub motif: usize,
    /// Where the hit is.
    pub range: GRange,
    /// The hit's score; overlapping hits keep the best.
    pub score: f64,
}

/// Proximity of one motif's hits to the anchors.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Proximity {
    /// The motif, as given in the hits.
    pub motif: usize,
    /// Hits tested, or anchors when counting anchors.
    pub nhits: usize,
    /// The most significant distance.
    pub best_dist: Option<i64>,
    /// Hits, or anchors, within the best distance.
    pub hits_near: Option<u64>,
    /// Hits, or anchors, beyond the best distance.
    pub hits_far: Option<u64>,
    /// The fraction of hits, or anchors, within the best distance.
    pub frac_near: Option<f64>,
    /// The fraction expected within the best distance under the null.
    pub expected_near: Option<f64>,
    /// The observed fraction divided by the expected one.
    pub enrichment: Option<f64>,
    /// The base 2 logarithm of the enrichment.
    pub log2_enrichment: Option<f64>,
    /// The P-value of the test. With the binomial null it is that of the best
    /// distance, Bonferroni-corrected for the distances tried.
    pub pvalue: f64,
    /// Benjamini-Hochberg q-value across motifs.
    pub qvalue: f64,
    /// Kolmogorov-Smirnov statistic, for the KS method only.
    pub ks_d: Option<f64>,
    /// Distance from each hit's middle to its nearest anchor, which is 0 inside
    /// the anchor and negative upstream of it in the anchor's orientation;
    /// `None` when no anchor shares its chromosome.
    pub dist_to_nearest: Vec<Option<f64>>,
}

struct Universe {
    segments: Vec<GRange>,
    total: f64,
    lengths: Vec<Option<i64>>,
}

impl Universe {
    fn len_of(&self, chrom: usize) -> Option<i64> {
        self.lengths.get(chrom).copied().flatten()
    }
}

/// The middle position of a range, taking the lower one for even widths.
fn middle(r: &GRange) -> GRange {
    let m = (r.start + r.end).div_euclid(2);
    GRange {
        start: m,
        end: m,
        ..*r
    }
}

/// Of hits that overlap one another on the same strand, or on either strand
/// with `ignore_strand`, only the best-scoring one is kept.
fn collapse(hits: &[ProximityHit], ignore_strand: bool) -> Vec<ProximityHit> {
    let mut keep = vec![false; hits.len()];
    let mut by_group: HashMap<(usize, Option<GStrand>), Vec<usize>> = HashMap::new();
    for (i, h) in hits.iter().enumerate() {
        let strand = (!ignore_strand).then_some(h.range.strand);
        by_group.entry((h.motif, strand)).or_default().push(i);
    }
    for idx in by_group.values() {
        let mut order = idx.clone();
        order.sort_by_key(|&i| {
            (
                hits[i].range.chrom,
                hits[i].range.start,
                hits[i].range.end,
                i,
            )
        });
        let mut k = 0;
        while k < order.len() {
            // A cluster of overlapping ranges.
            let first = &hits[order[k]].range;
            let (chrom, mut end) = (first.chrom, first.end);
            let mut j = k + 1;
            while j < order.len()
                && hits[order[j]].range.chrom == chrom
                && hits[order[j]].range.start <= end
            {
                end = end.max(hits[order[j]].range.end);
                j += 1;
            }
            let best = order[k..j]
                .iter()
                .copied()
                .max_by(|&a, &b| hits[a].score.total_cmp(&hits[b].score).then(b.cmp(&a)))
                .unwrap();
            keep[best] = true;
            k = j;
        }
    }
    hits.iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(h, _)| *h)
        .collect()
}

fn default_thresholds(u: &Universe) -> Vec<i64> {
    let min_w = u.segments.iter().map(|s| s.width()).min().unwrap_or(0);
    let cap = 10.max(min_w / 2);
    let grid: Vec<i64> = [
        10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000,
    ]
    .into_iter()
    .filter(|&g| g <= cap)
    .collect();
    if grid.is_empty() { vec![cap] } else { grid }
}

/// Fraction of the universe within `d` of an anchor.
fn zone_fraction(
    anchors: &[GRange],
    u: &Universe,
    d: i64,
    orientation: Orientation,
) -> (Vec<GRange>, f64) {
    let zones: Vec<GRange> = anchors
        .iter()
        .map(|a| zone(a, d, orientation, u.len_of(a.chrom)))
        .collect();
    let z = intersect(&reduce(&zones, 1), &u.segments);
    let f = total_width(&z) / u.total;
    (z, f)
}

#[derive(Clone)]
struct Scan {
    z: f64,
    p: f64,
    d: Option<i64>,
    k: u64,
    n: u64,
    expected: Option<f64>,
    frac: Option<f64>,
    enrichment: Option<f64>,
    n_tested: usize,
}

/// The best-scoring distance for anchors and one set of hits sorted by
/// chromosome and start.
fn scan(
    hits: &[GRange],
    anchors: &[GRange],
    u: &Universe,
    dists: &[i64],
    count: ProximityCount,
    orientation: Orientation,
) -> Scan {
    let n = match count {
        ProximityCount::Hits => hits.len(),
        ProximityCount::Anchors => anchors.len(),
    } as u64;
    let mut best = Scan {
        z: f64::NEG_INFINITY,
        p: 1.0,
        d: None,
        k: 0,
        n,
        expected: None,
        frac: None,
        enrichment: None,
        n_tested: 0,
    };
    if n == 0 || hits.is_empty() {
        return best;
    }
    let lambda = hits.len() as f64 / u.total;
    let hit_union = reduce(hits, 1);
    // Zones are reduced (sorted) for every distance; starting from sorted
    // anchors makes those sorts cheap. Anchor counts keep the given order,
    // since their zone fractions are summed in it.
    let mut sorted_anchors = anchors.to_vec();
    if count == ProximityCount::Hits {
        sorted_anchors.sort_by_key(|a| (a.chrom, a.start, a.end));
    }
    let anchors = &sorted_anchors[..];
    let mut n_tested = 0;
    for &d in dists {
        let (f, k) = match count {
            ProximityCount::Hits => {
                let (z, f) = zone_fraction(anchors, u, d, orientation);
                if !(f.is_finite() && f > 0.0 && f < 1.0) {
                    continue;
                }
                (f, count_overlapping(hits, &z) as u64)
            }
            ProximityCount::Anchors => {
                let zones: Vec<GRange> = anchors
                    .iter()
                    .map(|a| zone(a, d, orientation, u.len_of(a.chrom)))
                    .collect();
                let f = zones
                    .iter()
                    .map(|z| 1.0 - (-lambda * z.width() as f64).exp())
                    .sum::<f64>()
                    / zones.len() as f64;
                if !(f.is_finite() && f > 0.0 && f < 1.0) {
                    continue;
                }
                (
                    f,
                    zones.iter().filter(|z| overlaps_any(z, &hit_union)).count() as u64,
                )
            }
        };
        n_tested += 1;
        let nf = n as f64;
        let expected = nf * f;
        let p = pbinom_upper(k, n, f);
        let z = (k as f64 - expected) / (nf * f * (1.0 - f)).sqrt();
        if p < best.p || (p == best.p && z.is_finite() && z > best.z) {
            best = Scan {
                z,
                p,
                d: Some(d),
                k,
                n,
                expected: Some(expected),
                frac: Some(k as f64 / nf),
                enrichment: (expected > 0.0).then(|| k as f64 / expected),
                n_tested,
            };
        }
    }
    best.n_tested = n_tested;
    best
}

/// Signed distance from each hit to its nearest anchor.
fn signed_distances(hits: &[GRange], anchors: &[GRange], index: &AnchorIndex) -> Vec<Option<f64>> {
    hits.iter()
        .map(|h| {
            let (i, d) = index.nearest_position(h.chrom, h.start)?;
            let a = &anchors[i];
            let raw = h.start as f64 - (a.start + a.end) as f64 / 2.0;
            let off = if raw < 0.0 { -1.0 } else { 1.0 };
            let strand = if a.strand == GStrand::Minus {
                -1.0
            } else {
                1.0
            };
            Some(d as f64 * off * strand)
        })
        .collect()
}

/// The null distribution function of the distance to the nearest anchor,
/// tabulated on a log-spaced grid and interpolated linearly.
struct NullCdf {
    d: Vec<f64>,
    f: Vec<f64>,
}

impl NullCdf {
    fn new(anchors: &[GRange], u: &Universe, orientation: Orientation) -> NullCdf {
        let dmax = u.segments.iter().map(|s| s.width()).max().unwrap_or(1) as f64;
        let top = dmax.max(2.0).ln();
        let mut grid = vec![0i64];
        for i in 0..200 {
            let x = (top * i as f64 / 199.0).exp().round_ties_even() as i64;
            if !grid.contains(&x) {
                grid.push(x);
            }
        }
        let f = grid
            .iter()
            .map(|&d| {
                if d <= 0 {
                    0.0
                } else {
                    zone_fraction(anchors, u, d, orientation).1.clamp(0.0, 1.0)
                }
            })
            .collect();
        NullCdf {
            d: grid.iter().map(|&x| x as f64).collect(),
            f,
        }
    }

    fn at(&self, x: f64) -> f64 {
        if x < self.d[0] {
            return 0.0;
        }
        if x > *self.d.last().unwrap() {
            return 1.0;
        }
        let k = self.d.partition_point(|&g| g <= x);
        if k == 0 {
            return self.f[0];
        }
        let i = k - 1;
        if i + 1 >= self.d.len() || self.d[i] == x {
            return self.f[i];
        }
        self.f[i] + (self.f[i + 1] - self.f[i]) * (x - self.d[i]) / (self.d[i + 1] - self.d[i])
    }
}

/// One-sided Kolmogorov-Smirnov test for the "greater" alternative, using the
/// exact distribution for fewer than 100 values without ties and the asymptotic
/// one otherwise.
fn ks_greater(x: &[f64], cdf: impl Fn(f64) -> f64) -> (f64, f64) {
    let mut v = x.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len() as f64;
    let d = v
        .iter()
        .enumerate()
        .map(|(i, &xi)| 1.0 / n - (cdf(xi) - i as f64 / n))
        .fold(f64::NEG_INFINITY, f64::max);
    let ties = v.windows(2).any(|w| w[0] == w[1]);
    let p = if d <= 0.0 {
        1.0
    } else if d > 1.0 {
        0.0
    } else if v.len() < 100 && !ties {
        let lchoose = |n: f64, k: f64| {
            libm::lgamma(n + 1.0) - libm::lgamma(k + 1.0) - libm::lgamma(n - k + 1.0)
        };
        let jmax = (n * (1.0 - d)).floor() as i64;
        d * (0..=jmax)
            .map(|j| {
                let j = j as f64;
                (lchoose(n, j) + (n - j) * (1.0 - d - j / n).ln() + (j - 1.0) * (d + j / n).ln())
                    .exp()
            })
            .sum::<f64>()
    } else {
        (-2.0 * n * d * d).exp()
    };
    (d, p.clamp(0.0, 1.0))
}

/// Random anchors with the same widths and strands, placed uniformly at
/// legal positions in the universe.
fn permute_anchors(
    anchors: &[GRange],
    u: &Universe,
    per_chrom: bool,
    rng: &mut rand::rngs::StdRng,
) -> Result<Vec<GRange>> {
    let mut out = Vec::with_capacity(anchors.len());
    for a in anchors {
        let w = a.width();
        let segs: Vec<&GRange> = u
            .segments
            .iter()
            .filter(|s| !per_chrom || s.chrom == a.chrom)
            .collect();
        let legal: Vec<f64> = segs
            .iter()
            .map(|s| (s.width() - w + 1).max(0) as f64)
            .collect();
        let total: f64 = legal.iter().sum();
        if total <= 0.0 {
            return arg_err(format!(
                "could not place an anchor of width {w} in the universe; the universe is too small"
            ));
        }
        let mut r = rng.random::<f64>() * total;
        let mut pick = legal.len() - 1;
        for (i, &l) in legal.iter().enumerate() {
            if r < l {
                pick = i;
                break;
            }
            r -= l;
        }
        let off = (rng.random::<f64>() * legal[pick]).floor() as i64;
        let start = segs[pick].start + off;
        out.push(GRange {
            chrom: segs[pick].chrom,
            start,
            end: start + w - 1,
            strand: a.strand,
        });
    }
    Ok(out)
}

/// Test each motif's hits for proximity to `anchors`. `universe` is the scanned
/// regions, which form the null support, and `chrom_lengths` the chromosome
/// lengths where known, to which zones are clipped. Only chromosomes carrying
/// anchors are used. Returns one result per motif with hits, in increasing
/// order of motif, unfiltered.
pub fn motif_proximity(
    hits: &[ProximityHit],
    anchors: &[GRange],
    universe: &[GRange],
    chrom_lengths: &[Option<i64>],
    opts: &ProximityOptions,
) -> Result<Vec<Proximity>> {
    if opts.method == ProximityMethod::Ks && opts.count == ProximityCount::Anchors {
        return arg_err(
            "the KS test is not available when counting anchors; use the permutation test",
        );
    }
    if opts.method == ProximityMethod::Ks && opts.orientation != Orientation::Undirected {
        return arg_err("the KS test is not available for directional zones");
    }
    if opts.n_perm < 1 {
        return arg_err("n_perm must be at least 1");
    }
    if opts.orientation != Orientation::Undirected
        && !anchors.iter().any(|a| a.strand != GStrand::Any)
    {
        return arg_err("directional zones need stranded anchors");
    }
    let anchor_chroms: BTreeSet<usize> = anchors.iter().map(|a| a.chrom).collect();
    let uni_chroms: BTreeSet<usize> = universe.iter().map(|r| r.chrom).collect();
    if anchor_chroms.intersection(&uni_chroms).next().is_none() {
        return arg_err("the anchors and the universe share no chromosomes");
    }
    let segments = reduce(
        &universe
            .iter()
            .filter(|r| anchor_chroms.contains(&r.chrom))
            .copied()
            .collect::<Vec<_>>(),
        1,
    );
    let u = Universe {
        total: total_width(&segments),
        segments,
        lengths: chrom_lengths.to_vec(),
    };
    let anchors: Vec<GRange> = anchors
        .iter()
        .filter(|a| uni_chroms.contains(&a.chrom))
        .copied()
        .collect();
    // Hits are placed at their middle positions, so that "near" means the
    // same thing for a hit as for a position of the universe.
    let mut hits: Vec<ProximityHit> = hits
        .iter()
        .filter(|h| overlaps_any(&middle(&h.range), &u.segments))
        .copied()
        .collect();
    if hits.is_empty() || anchors.is_empty() || u.total <= 0.0 {
        return Ok(Vec::new());
    }
    if opts.count == ProximityCount::Hits {
        hits = collapse(&hits, opts.ignore_strand);
    }
    for h in hits.iter_mut() {
        h.range = middle(&h.range);
    }
    let dists: Vec<i64> = match &opts.dist_thresholds {
        None => default_thresholds(&u),
        Some(d) => {
            if d.iter().any(|&x| x <= 0) {
                return arg_err("distance thresholds must be positive");
            }
            let mut d = d.clone();
            d.sort_unstable();
            d.dedup();
            d
        }
    };
    let null_cdf =
        (opts.method == ProximityMethod::Ks).then(|| NullCdf::new(&anchors, &u, opts.orientation));
    let index = AnchorIndex::new(&anchors);

    let motif_ids: Vec<usize> = hits
        .iter()
        .map(|h| h.motif)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut out: Vec<Proximity> = with_threads(opts.nthreads, || {
        motif_ids
            .par_iter()
            .map(|&m| -> Result<Proximity> {
                let hm: Vec<GRange> = hits
                    .iter()
                    .filter(|h| h.motif == m)
                    .map(|h| h.range)
                    .collect();
                let mut sorted = hm.clone();
                sorted.sort_by_key(|h| (h.chrom, h.start, h.end));
                let mut res = Proximity {
                    motif: m,
                    nhits: match opts.count {
                        ProximityCount::Hits => hm.len(),
                        ProximityCount::Anchors => anchors.len(),
                    },
                    best_dist: None,
                    hits_near: None,
                    hits_far: None,
                    frac_near: None,
                    expected_near: None,
                    enrichment: None,
                    log2_enrichment: None,
                    pvalue: 1.0,
                    qvalue: 1.0,
                    ks_d: None,
                    dist_to_nearest: signed_distances(&hm, &anchors, &index),
                };
                let fill = |res: &mut Proximity, sc: &Scan, p: f64| {
                    res.best_dist = sc.d;
                    res.hits_near = Some(sc.k);
                    res.hits_far = Some(sc.n - sc.k);
                    res.frac_near = sc.frac;
                    res.expected_near = sc.expected;
                    res.enrichment = sc.enrichment;
                    res.log2_enrichment = sc.enrichment.filter(|&e| e > 0.0).map(f64::log2);
                    res.pvalue = p;
                };
                match opts.method {
                    ProximityMethod::Binomial => {
                        let sc = scan(&sorted, &anchors, &u, &dists, opts.count, opts.orientation);
                        let p = (sc.p * sc.n_tested.max(1) as f64).min(1.0);
                        fill(&mut res, &sc, p);
                    }
                    ProximityMethod::Ks => {
                        let d: Vec<f64> = res
                            .dist_to_nearest
                            .iter()
                            .flatten()
                            .map(|x| x.abs())
                            .collect();
                        if !d.is_empty() {
                            let cdf = null_cdf.as_ref().unwrap();
                            let (stat, p) = ks_greater(&d, |x| cdf.at(x));
                            res.ks_d = Some(stat);
                            res.pvalue = p;
                        }
                    }
                    ProximityMethod::Permutation => {
                        let obs = scan(&sorted, &anchors, &u, &dists, opts.count, opts.orientation);
                        let mut rng = rng_for(opts.seed, m);
                        let mut above = 0usize;
                        for _ in 0..opts.n_perm {
                            let ra = permute_anchors(&anchors, &u, opts.perm_per_chrom, &mut rng)?;
                            let z = scan(&sorted, &ra, &u, &dists, opts.count, opts.orientation).z;
                            if z >= obs.z {
                                above += 1;
                            }
                        }
                        let p = (1 + above) as f64 / (1 + opts.n_perm) as f64;
                        fill(&mut res, &obs, p);
                    }
                }
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()
    })??;
    let q = p_adjust_bh(&out.iter().map(|p| p.pvalue).collect::<Vec<_>>());
    for (p, q) in out.iter_mut().zip(q) {
        p.qvalue = q;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ks_matches_r() {
        // ks.test(c(0.1, 0.25, 0.3, 0.42, 0.9), "punif", alternative = "greater")
        let (d, p) = ks_greater(&[0.1, 0.25, 0.3, 0.42, 0.9], |x| x.clamp(0.0, 1.0));
        assert!((d - 0.38).abs() < 1e-12, "{d}");
        assert!((p - 0.1837559872).abs() < 1e-9, "{p}");
    }
}
