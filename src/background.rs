//! Composition-matched background sequences, drawn by [`match_bkg`].
//!
//! HOMER-style matching: the universe of candidate sequences is binned by GC
//! fraction into equal-width bins over [0, 1], by length into quantile bins and
//! by any extra covariates, and each target draws background sequences from its
//! own bin, widening the search ring by ring when the bin runs short.

use std::collections::HashMap;

use rand::RngExt;
use rand::rngs::StdRng;

use crate::error::{Result, arg_err};
use crate::ranges::{GRange, overlaps_any, reduce};
use crate::shuffle::rng_for;

/// Fraction of G and C in a sequence, counting either case.
pub fn gc_fraction(seq: &str) -> f64 {
    let gc = seq
        .bytes()
        .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
        .count();
    gc as f64 / seq.len().max(1) as f64
}

/// How a matching axis is binned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinType {
    /// Breaks at quantiles of the universe values.
    Quantile,
    /// Equal-width breaks over the range of all values.
    EqualWidth,
}

impl BinType {
    /// The name [`str::parse`] reads, such as `"quantile"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            BinType::Quantile => "quantile",
            BinType::EqualWidth => "equalwidth",
        }
    }
}

impl std::str::FromStr for BinType {
    type Err = crate::error::Error;

    /// Parse `"quantile"` or `"equalwidth"`.
    fn from_str(s: &str) -> Result<BinType> {
        match s {
            "quantile" => Ok(BinType::Quantile),
            "equalwidth" => Ok(BinType::EqualWidth),
            _ => arg_err(format!(
                "unknown bin type '{s}'; expected 'quantile' or 'equalwidth'"
            )),
        }
    }
}

impl std::fmt::Display for BinType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// An extra per-sequence value to match on.
#[derive(Debug, Clone, PartialEq)]
pub struct Covariate {
    /// The covariate's name, used in error messages.
    pub name: String,
    /// The value for each target.
    pub targets: Vec<f64>,
    /// The value for each universe sequence.
    pub universe: Vec<f64>,
    /// The number of bins.
    pub n_bins: usize,
    /// How the bins are made.
    pub bin_type: BinType,
}

/// Per-target and per-universe 0-based bins on one axis.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Axis {
    /// The bin of each target.
    pub bin_targets: Vec<usize>,
    /// The bin of each universe sequence.
    pub bin_universe: Vec<usize>,
    /// The number of bins.
    pub n_bins: usize,
}

/// `n` evenly spaced values from `from` to `to`, the last one set to `to`
/// exactly.
fn r_seq(from: f64, to: f64, n: usize) -> Vec<f64> {
    if n == 1 {
        return vec![from];
    }
    if from == to {
        return vec![from; n];
    }
    let n1 = (n - 1) as f64;
    let mut v: Vec<f64> = (0..n - 1)
        .map(|i| from + i as f64 * ((to - from) / n1))
        .collect();
    v.push(to);
    v
}

/// Sample quantiles of `x` at `probs`, interpolating linearly between the
/// sorted values so that the quantile at `p` sits at position `1 + (n - 1) p`,
/// which is Hyndman and Fan's type 7.
fn quantile7(x: &[f64], probs: &[f64]) -> Vec<f64> {
    let mut v = x.to_vec();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    probs
        .iter()
        .map(|&p| {
            let index = 1.0 + (n.saturating_sub(1)) as f64 * p;
            let (lo, hi) = (index.floor(), index.ceil());
            let qs = v[lo as usize - 1];
            let xh = v[hi as usize - 1];
            if index > lo && xh != qs {
                let h = index - lo;
                (1.0 - h) * qs + h * xh
            } else {
                qs
            }
        })
        .collect()
}

/// The 0-based bin of `x` among the increasing `breaks`: bin `k` runs from
/// `breaks[k]` up to but not including `breaks[k + 1]`, except that the last
/// bin also holds its upper break. Values outside the breaks go to the first
/// or last bin.
fn bin_of(x: f64, breaks: &[f64]) -> usize {
    let nb = breaks.len() - 1;
    let mut k = breaks.partition_point(|&b| b <= x);
    if x == breaks[nb] {
        k = nb;
    }
    k.clamp(1, nb) - 1
}

/// Bin targets and universe on one axis. Equal-width breaks span `range`, which
/// defaults to all the values, while quantile breaks come from the universe. A
/// constant axis collapses to one bin.
#[doc(hidden)]
pub fn match_axis(
    targets: &[f64],
    universe: &[f64],
    n_bins: usize,
    bin_type: BinType,
    range: Option<(f64, f64)>,
) -> Axis {
    let all = targets.iter().chain(universe);
    let lo = all.clone().copied().fold(f64::INFINITY, f64::min);
    let hi = all.copied().fold(f64::NEG_INFINITY, f64::max);
    let breaks = match bin_type {
        BinType::EqualWidth => {
            let (a, b) = range.unwrap_or((lo, hi));
            let br = r_seq(a, b, n_bins + 1);
            let mut u = br.clone();
            u.dedup();
            if u.len() < 2 {
                vec![lo - 1.0, hi + 1.0]
            } else {
                br
            }
        }
        BinType::Quantile => {
            let probs = r_seq(0.0, 1.0, n_bins + 1);
            let mut br = quantile7(universe, &probs);
            br.dedup();
            if br.len() < 2 {
                vec![lo - 1.0, hi + 1.0]
            } else {
                let last = br.len() - 1;
                br[0] = br[0].min(lo) - 1.0;
                br[last] = br[last].max(hi) + 1.0;
                br
            }
        }
    };
    Axis {
        bin_targets: targets.iter().map(|&x| bin_of(x, &breaks)).collect(),
        bin_universe: universe.iter().map(|&x| bin_of(x, &breaks)).collect(),
        n_bins: breaks.len() - 1,
    }
}

/// Options for [`match_bkg`].
#[derive(Debug, Clone, PartialEq)]
pub struct MatchOptions {
    /// Universe sequences to draw for each target; 1 by default.
    pub n_per_target: usize,
    /// Equal-width GC content bins between 0 and 1; 20 by default.
    pub n_bins_gc: usize,
    /// Quantile bins of sequence length; 10 by default.
    pub n_bins_length: usize,
    /// Draw each universe sequence at most once across all targets. Either
    /// way, the draws for one target are distinct, coming from nearby bins
    /// when its own bin holds too few and then with replacement when those
    /// run out too.
    pub unique: bool,
    /// Seed for the random draws.
    pub seed: u64,
}

impl Default for MatchOptions {
    fn default() -> Self {
        MatchOptions {
            n_per_target: 1,
            n_bins_gc: 20,
            n_bins_length: 10,
            unique: true,
            seed: 1,
        }
    }
}

/// Result of [`match_bkg`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Matched {
    /// The universe sequence drawn for each target, `n_per_target` per
    /// target in target order.
    pub universe_index: Vec<usize>,
    /// Targets that needed drawing with replacement from nearby bins.
    pub fallback_targets: usize,
    /// The GC fraction of each target.
    pub gc_targets: Vec<f64>,
    /// The GC fraction of each universe sequence.
    pub gc_universe: Vec<f64>,
}

/// Offsets with absolute values summing to `ring`, forming a Manhattan ring.
fn ring_offsets(ring: i64, n_axes: usize) -> Vec<Vec<i64>> {
    fn rec(axis: usize, n: usize, remaining: i64, acc: &mut Vec<i64>, out: &mut Vec<Vec<i64>>) {
        if axis + 1 == n {
            for v in if remaining == 0 {
                vec![0]
            } else {
                vec![remaining, -remaining]
            } {
                acc.push(v);
                out.push(acc.clone());
                acc.pop();
            }
            return;
        }
        for m in 0..=remaining {
            for v in if m == 0 { vec![0] } else { vec![m, -m] } {
                acc.push(v);
                rec(axis + 1, n, remaining - m, acc, out);
                acc.pop();
            }
        }
    }
    let mut out = Vec::new();
    rec(0, n_axes, ring, &mut Vec::new(), &mut out);
    out
}

/// `k` distinct positions in `0..n`, by Floyd's algorithm.
fn distinct_positions(rng: &mut StdRng, n: usize, k: usize) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::with_capacity(k);
    for j in n - k..n {
        let t = rng.random_range(0..=j);
        out.push(if out.contains(&t) { j } else { t });
    }
    out
}

/// The cell and the index within it of position `pos` in the concatenation
/// of `pools[c]` over `cells`.
fn locate(cells: &[usize], pools: &[Vec<usize>], mut pos: usize) -> (usize, usize) {
    for &c in cells {
        if pos < pools[c].len() {
            return (c, pos);
        }
        pos -= pools[c].len();
    }
    unreachable!("position past the end of the pool")
}

/// Draw `opts.n_per_target` universe sequences for each target with the
/// same GC fraction, length and covariate bins. When a target's bin runs
/// short, up to `max(2, axes)`
/// rings of neighbouring bins are searched; after that the rest are drawn
/// with replacement from the bins already searched.
pub fn match_bkg<S: AsRef<str>, U: AsRef<str>>(
    targets: &[S],
    universe: &[U],
    covariates: &[Covariate],
    opts: &MatchOptions,
) -> Result<Matched> {
    if targets.is_empty() {
        return arg_err("no target sequences given");
    }
    if universe.is_empty() {
        return arg_err("the universe is empty");
    }
    if universe.iter().any(|s| s.as_ref().is_empty()) {
        return arg_err("the universe must not contain empty sequences");
    }
    if opts.n_per_target < 1 {
        return arg_err("n_per_target must be at least 1");
    }
    if opts.n_bins_gc < 2 || opts.n_bins_length < 2 || covariates.iter().any(|c| c.n_bins < 2) {
        return arg_err("every axis needs at least 2 bins");
    }
    for c in covariates {
        if c.targets.len() != targets.len() || c.universe.len() != universe.len() {
            return arg_err(format!(
                "covariate '{}' needs one value per target and per universe sequence",
                c.name
            ));
        }
        if c.targets.iter().chain(&c.universe).any(|v| !v.is_finite()) {
            return arg_err(format!(
                "covariate '{}' has missing or infinite values",
                c.name
            ));
        }
    }
    let gc_t: Vec<f64> = targets.iter().map(|s| gc_fraction(s.as_ref())).collect();
    let gc_u: Vec<f64> = universe.iter().map(|s| gc_fraction(s.as_ref())).collect();
    let len_t: Vec<f64> = targets.iter().map(|s| s.as_ref().len() as f64).collect();
    let len_u: Vec<f64> = universe.iter().map(|s| s.as_ref().len() as f64).collect();
    let mut axes = vec![
        match_axis(
            &gc_t,
            &gc_u,
            opts.n_bins_gc,
            BinType::EqualWidth,
            Some((0.0, 1.0)),
        ),
        match_axis(&len_t, &len_u, opts.n_bins_length, BinType::Quantile, None),
    ];
    for c in covariates {
        axes.push(match_axis(
            &c.targets,
            &c.universe,
            c.n_bins,
            c.bin_type,
            None,
        ));
    }
    let n_axes = axes.len();
    let key = |coords: &[usize]| -> usize {
        coords
            .iter()
            .zip(&axes)
            .fold(0, |k, (&c, a)| k * a.n_bins + c)
    };
    // Cells are disjoint, so the cells of a ring never share sequences.
    // `free` holds each cell's sequences not yet drawn (all of them when
    // draws need not be unique), so a draw only touches a ring's few cells.
    let mut cell_ids: HashMap<usize, usize> = HashMap::new();
    let mut all: Vec<Vec<usize>> = Vec::new();
    for i in 0..universe.len() {
        let coords: Vec<usize> = axes.iter().map(|a| a.bin_universe[i]).collect();
        let id = *cell_ids.entry(key(&coords)).or_insert_with(|| {
            all.push(Vec::new());
            all.len() - 1
        });
        all[id].push(i);
    }
    let mut free = all.clone();
    let max_ring = 2.max(n_axes as i64);
    let rings: Vec<Vec<Vec<i64>>> = (0..=max_ring).map(|r| ring_offsets(r, n_axes)).collect();

    let mut rng = rng_for(opts.seed, 0);
    let mut out = Vec::with_capacity(targets.len() * opts.n_per_target);
    let mut fallbacks = 0;
    for t in 0..targets.len() {
        let home: Vec<i64> = axes.iter().map(|a| a.bin_targets[t] as i64).collect();
        let mut picked: Vec<usize> = Vec::new();
        let mut searched: Vec<usize> = Vec::new();
        for offsets in &rings {
            let ring: Vec<usize> = offsets
                .iter()
                .filter_map(|off| {
                    let coords: Option<Vec<usize>> = home
                        .iter()
                        .zip(off)
                        .zip(&axes)
                        .map(|((&h, &o), a)| {
                            let c = h + o;
                            (c >= 0 && (c as usize) < a.n_bins).then_some(c as usize)
                        })
                        .collect();
                    coords.and_then(|c| cell_ids.get(&key(&c)).copied())
                })
                .collect();
            searched.extend(&ring);
            let total: usize = ring.iter().map(|&c| free[c].len()).sum();
            let take = total.min(opts.n_per_target - picked.len());
            let mut chosen: Vec<(usize, usize)> = distinct_positions(&mut rng, total, take)
                .into_iter()
                .map(|p| locate(&ring, &free, p))
                .collect();
            picked.extend(chosen.iter().map(|&(c, k)| free[c][k]));
            if opts.unique {
                // Highest index first, so earlier ones stay put.
                chosen.sort_unstable_by(|a, b| b.cmp(a));
                for (c, k) in chosen {
                    free[c].swap_remove(k);
                }
            }
            if picked.len() >= opts.n_per_target {
                break;
            }
        }
        let need = opts.n_per_target - picked.len();
        if need > 0 {
            // Draw the rest with replacement from every searched cell,
            // drawn already or not.
            fallbacks += 1;
            let total: usize = searched.iter().map(|&c| all[c].len()).sum();
            for _ in 0..need {
                picked.push(if total == 0 {
                    rng.random_range(0..universe.len())
                } else {
                    let (c, k) = locate(&searched, &all, rng.random_range(0..total));
                    all[c][k]
                });
            }
        }
        out.extend(picked);
    }
    Ok(Matched {
        universe_index: out,
        fallback_targets: fallbacks,
        gc_targets: gc_t,
        gc_universe: gc_u,
    })
}

/// Random genomic windows with widths drawn from `target_widths`, for a
/// background universe for [`match_bkg`].
/// Chromosomes are chosen in proportion to their lengths; windows that do
/// not fit are dropped, as are windows overlapping `exclude`. Returns the
/// first `n_candidates` windows.
pub fn genome_windows(
    chrom_lengths: &[i64],
    target_widths: &[i64],
    n_candidates: usize,
    exclude: &[GRange],
    seed: u64,
) -> Result<Vec<GRange>> {
    if target_widths.is_empty() || target_widths.iter().any(|&w| w < 1) {
        return arg_err("target sequences must not be empty");
    }
    if n_candidates < target_widths.len() {
        return arg_err("n_candidates must be at least the number of targets");
    }
    let total: f64 = chrom_lengths.iter().map(|&l| l.max(0) as f64).sum();
    if total <= 0.0 {
        return arg_err("no usable chromosomes");
    }
    let n_try = if exclude.is_empty() {
        n_candidates
    } else {
        2 * n_candidates
    };
    let exclude = reduce(exclude, 1);
    let mut rng = rng_for(seed, 0);
    let mut out = Vec::new();
    for _ in 0..n_try {
        let mut r = rng.random::<f64>() * total;
        let mut chrom = chrom_lengths.len() - 1;
        for (i, &l) in chrom_lengths.iter().enumerate() {
            if r < l.max(0) as f64 {
                chrom = i;
                break;
            }
            r -= l.max(0) as f64;
        }
        let w = target_widths[rng.random_range(0..target_widths.len())];
        let n_starts = chrom_lengths[chrom] - w + 1;
        if n_starts < 1 {
            continue;
        }
        let start = (rng.random::<f64>() * n_starts as f64).floor() as i64;
        let g = GRange::new(chrom, start, start + w);
        if !overlaps_any(&g, &exclude) {
            out.push(g);
        }
    }
    if out.len() < target_widths.len() {
        return arg_err(format!(
            "could not sample enough windows outside the excluded regions (got {}, need {})",
            out.len(),
            target_widths.len()
        ));
    }
    out.truncate(n_candidates);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_positions_cover_the_range_evenly() {
        let mut rng = rng_for(3, 0);
        let mut seen = [0usize; 10];
        for _ in 0..20_000 {
            let mut p = distinct_positions(&mut rng, 10, 3);
            for &x in &p {
                seen[x] += 1;
            }
            p.sort_unstable();
            p.dedup();
            assert_eq!(p.len(), 3);
        }
        // Each position is drawn in 3 of every 10 draws.
        assert!(
            seen.iter().all(|&n| (5_600..6_400).contains(&n)),
            "{seen:?}"
        );
        assert_eq!(distinct_positions(&mut rng, 4, 4).len(), 4);
        assert!(distinct_positions(&mut rng, 0, 0).is_empty());
        let pools = vec![vec![7, 8], vec![], vec![9]];
        assert_eq!(locate(&[0, 1, 2], &pools, 2), (2, 0));
    }

    #[test]
    fn ring_offsets_cover_the_ring() {
        for (ring, n) in [(0, 2), (1, 2), (2, 2), (2, 3), (3, 4)] {
            let offs = ring_offsets(ring, n);
            assert!(
                offs.iter()
                    .all(|o| o.len() == n && o.iter().map(|x| x.abs()).sum::<i64>() == ring)
            );
            let uniq: std::collections::HashSet<_> = offs.iter().collect();
            assert_eq!(
                uniq.len(),
                offs.len(),
                "duplicates for ring {ring}, {n} axes"
            );
        }
        // 4 * ring cells on a 2-D ring.
        assert_eq!(ring_offsets(3, 2).len(), 12);
    }

    #[test]
    fn quantiles_match_r() {
        // quantile(c(5, 1, 9, 3, 7, 2), seq(0, 1, length.out = 5))
        let q = quantile7(&[5.0, 1.0, 9.0, 3.0, 7.0, 2.0], &r_seq(0.0, 1.0, 5));
        let want = [1.0, 2.25, 4.0, 6.5, 9.0];
        for (a, b) in q.iter().zip(want) {
            assert!((a - b).abs() < 1e-12, "{q:?}");
        }
    }
}
