//! Hierarchical clustering.
//!
//! [`hclust`] is a port of R's `hclust()` (F. Murtagh's `HCLUST` and
//! `HCASS2` Fortran routines in R's `stats` package), so merges, heights and
//! leaf order match R exactly, ties included.

use crate::error::{Result, arg_err};

/// Agglomeration method, named as in R's `hclust()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    WardD,
    Single,
    #[default]
    Complete,
    Average,
    Mcquitty,
    Median,
    Centroid,
    WardD2,
}

impl Method {
    pub fn parse(s: &str) -> Result<Method> {
        Ok(match s {
            "ward.D" => Method::WardD,
            "single" => Method::Single,
            "complete" => Method::Complete,
            "average" => Method::Average,
            "mcquitty" => Method::Mcquitty,
            "median" => Method::Median,
            "centroid" => Method::Centroid,
            "ward.D2" => Method::WardD2,
            _ => {
                return arg_err(format!(
                    "unknown clustering method '{s}' (use ward.D, single, complete, average, \
                     mcquitty, median, centroid or ward.D2)"
                ));
            }
        })
    }

    fn code(self) -> u8 {
        match self {
            Method::WardD => 1,
            Method::Single => 2,
            Method::Complete => 3,
            Method::Average => 4,
            Method::Mcquitty => 5,
            Method::Median => 6,
            Method::Centroid => 7,
            Method::WardD2 => 8,
        }
    }
}

/// A clustering tree in R's `hclust` form.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    /// Merge steps, as in R: a negative value `-i` is leaf `i` (1-based),
    /// a positive value `j` is the cluster formed at step `j`.
    pub merge: Vec<[i64; 2]>,
    /// Height of each merge.
    pub height: Vec<f64>,
    /// Leaves (0-based) in dendrogram order.
    pub order: Vec<usize>,
}

/// Cluster `n` objects from a distance matrix. Only the lower triangle
/// (`dist[i][j]` with `i > j`) is used, as R's `as.dist()` does.
pub fn hclust(dist: &[Vec<f64>], method: Method) -> Result<Tree> {
    let n = dist.len();
    if n < 2 {
        return arg_err("need at least 2 objects to cluster");
    }
    if dist.iter().any(|r| r.len() != n) {
        return arg_err("the distance matrix must be square");
    }
    // R's dist vector: lower triangle, column by column, 1-based.
    let mut diss = vec![0.0f64];
    for j in 0..n {
        for row in dist.iter().skip(j + 1) {
            let d = row[j];
            if !d.is_finite() {
                return arg_err("distances must be finite");
            }
            diss.push(d);
        }
    }
    let (ia, ib, crit) = hclust_fortran(n, &mut diss, method.code());
    let (merge, order) = hcass2(n, &ia, &ib);
    Ok(Tree {
        merge,
        height: crit[1..n].to_vec(),
        order,
    })
}

/// Index of (i, j), i < j, both 1-based, in the 1-based dist vector.
fn ioffst(n: usize, i: usize, j: usize) -> usize {
    j + (i - 1) * n - (i * (i + 1)) / 2
}

#[allow(clippy::needless_range_loop)]
fn hclust_fortran(n: usize, diss: &mut [f64], iopt: u8) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    const INF: f64 = 1e300;
    let mut ia = vec![0usize; n + 1];
    let mut ib = vec![0usize; n + 1];
    let mut crit = vec![0.0f64; n + 1];
    let mut membr = vec![1.0f64; n + 1];
    let mut nn = vec![0usize; n + 1];
    let mut disnn = vec![0.0f64; n + 1];
    let mut flag = vec![true; n + 1];
    let (mut im, mut jj, mut jm) = (0usize, 0usize, 0usize);
    let mut ncl = n;

    if iopt == 8 {
        for d in diss.iter_mut().skip(1) {
            *d *= *d;
        }
    }

    for i in 1..n {
        let mut dmin = INF;
        for j in i + 1..=n {
            let ind = ioffst(n, i, j);
            if dmin > diss[ind] {
                dmin = diss[ind];
                jm = j;
            }
        }
        nn[i] = jm;
        disnn[i] = dmin;
    }

    loop {
        let mut dmin = INF;
        for i in 1..n {
            if flag[i] && disnn[i] < dmin {
                dmin = disnn[i];
                im = i;
                jm = nn[i];
            }
        }
        ncl -= 1;

        let i2 = im.min(jm);
        let j2 = im.max(jm);
        ia[n - ncl] = i2;
        ib[n - ncl] = j2;
        let is_ward = iopt == 1 || iopt == 8;
        if iopt == 8 {
            dmin = dmin.sqrt();
        }
        crit[n - ncl] = dmin;
        flag[j2] = false;

        dmin = INF;
        for k in 1..=n {
            if flag[k] && k != i2 {
                let ind1 = if i2 < k {
                    ioffst(n, i2, k)
                } else {
                    ioffst(n, k, i2)
                };
                let ind2 = if j2 < k {
                    ioffst(n, j2, k)
                } else {
                    ioffst(n, k, j2)
                };
                let d12 = diss[ioffst(n, i2, j2)];
                let (mi, mj, mk) = (membr[i2], membr[j2], membr[k]);
                if is_ward {
                    diss[ind1] = (mi + mk) * diss[ind1] + (mj + mk) * diss[ind2] - mk * d12;
                    diss[ind1] /= mi + mj + mk;
                } else {
                    match iopt {
                        2 => diss[ind1] = diss[ind1].min(diss[ind2]),
                        3 => diss[ind1] = diss[ind1].max(diss[ind2]),
                        4 => diss[ind1] = (mi * diss[ind1] + mj * diss[ind2]) / (mi + mj),
                        5 => diss[ind1] = (diss[ind1] + diss[ind2]) / 2.0,
                        6 => diss[ind1] = ((diss[ind1] + diss[ind2]) - d12 / 2.0) / 2.0,
                        7 => {
                            diss[ind1] = (mi * diss[ind1] + mj * diss[ind2]
                                - mi * mj * d12 / (mi + mj))
                                / (mi + mj)
                        }
                        _ => unreachable!(),
                    }
                }
                if i2 < k {
                    if diss[ind1] < dmin {
                        dmin = diss[ind1];
                        jj = k;
                    }
                } else if diss[ind1] < disnn[k] {
                    disnn[k] = diss[ind1];
                    nn[k] = i2;
                }
            }
        }
        membr[i2] += membr[j2];
        disnn[i2] = dmin;
        nn[i2] = jj;

        for i in 1..n {
            if flag[i] && (nn[i] == i2 || nn[i] == j2) {
                let mut dmin = INF;
                for j in i + 1..=n {
                    if flag[j] {
                        let ind = ioffst(n, i, j);
                        if diss[ind] < dmin {
                            dmin = diss[ind];
                            jj = j;
                        }
                    }
                }
                nn[i] = jj;
                disnn[i] = dmin;
            }
        }

        if ncl <= 1 {
            break;
        }
    }
    (ia, ib, crit)
}

// The loop bound changes only right before `break`, as in the Fortran.
#[allow(clippy::needless_range_loop, clippy::mut_range_bound)]
fn hcass2(n: usize, ia: &[usize], ib: &[usize]) -> (Vec<[i64; 2]>, Vec<usize>) {
    let mut iia: Vec<i64> = ia.iter().map(|&x| x as i64).collect();
    let mut iib: Vec<i64> = ib.iter().map(|&x| x as i64).collect();
    for i in 1..=n.saturating_sub(2) {
        let k = ia[i].min(ib[i]);
        for j in i + 1..n {
            if ia[j] == k {
                iia[j] = -(i as i64);
            }
            if ib[j] == k {
                iib[j] = -(i as i64);
            }
        }
    }
    for i in 1..n {
        iia[i] = -iia[i];
        iib[i] = -iib[i];
    }
    for i in 1..n {
        if iia[i] > 0 && iib[i] < 0 {
            std::mem::swap(&mut iia[i], &mut iib[i]);
        }
        if iia[i] > 0 && iib[i] > 0 {
            let (k1, k2) = (iia[i].min(iib[i]), iia[i].max(iib[i]));
            iia[i] = k1;
            iib[i] = k2;
        }
    }

    let mut order = vec![0i64; n + 1];
    order[1] = iia[n - 1];
    order[2] = iib[n - 1];
    let mut loc = 2;
    for i in (1..=n.saturating_sub(2)).rev() {
        for j in 1..=loc {
            if order[j] == i as i64 {
                order[j] = iia[i];
                loc += 1;
                if j == loc - 1 {
                    order[loc] = iib[i];
                } else {
                    for k in (j + 2..=loc).rev() {
                        order[k] = order[k - 1];
                    }
                    order[j + 1] = iib[i];
                }
                break;
            }
        }
    }
    let merge = (1..n).map(|i| [iia[i], iib[i]]).collect();
    let order = (1..=n).map(|i| (-order[i] - 1) as usize).collect();
    (merge, order)
}

/// C's `%.10g`, as used by ape's `write.tree()`.
fn format_g10(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let sci = format!("{x:.9e}");
    let (_, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let strip = |s: String| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    };
    if (-4..10).contains(&exp) {
        strip(format!("{x:.*}", (9 - exp) as usize))
    } else {
        let (mant, _) = sci.split_once('e').unwrap();
        format!(
            "{}e{}{:02}",
            strip(mant.to_string()),
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    }
}

/// Clean a tip label the way ape's `checkLabel()` does.
fn newick_label(s: &str) -> String {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace() || c == '(');
    let s = s.trim_end_matches(|c: char| c.is_ascii_whitespace() || c == ')');
    s.chars()
        .map(|c| match c {
            c if c.is_ascii_whitespace() => '_',
            ',' | ':' | ';' | '(' | ')' => '-',
            c => c,
        })
        .collect()
}

impl Tree {
    /// The tree in Newick format, as `ape::write.tree(ape::as.phylo(tree))`
    /// writes it: branch lengths are half the height differences, so leaves
    /// sit at half the merge height from each other.
    pub fn newick(&self, labels: &[String]) -> String {
        let mut strs: Vec<String> = Vec::with_capacity(self.merge.len());
        for (i, m) in self.merge.iter().enumerate() {
            let h = self.height[i];
            let child = |c: i64| -> String {
                if c < 0 {
                    let label = labels
                        .get((-c - 1) as usize)
                        .map(|l| newick_label(l))
                        .unwrap_or_else(|| (-c).to_string());
                    format!("{label}:{}", format_g10(h / 2.0))
                } else {
                    let y = (c - 1) as usize;
                    format!("{}:{}", strs[y], format_g10((h - self.height[y]) / 2.0))
                }
            };
            let s = format!("({},{})", child(m[0]), child(m[1]));
            strs.push(s);
        }
        format!("{};", strs.last().cloned().unwrap_or_default())
    }

    /// The tree as a SciPy linkage matrix: one row per merge of
    /// `[a, b, height, size]`, where leaves are `0..n` and the cluster from
    /// merge `i` (0-based) is `n + i`.
    pub fn linkage(&self) -> Vec<[f64; 4]> {
        let n = self.merge.len() + 1;
        let mut sizes: Vec<f64> = Vec::with_capacity(self.merge.len());
        let id = |c: i64| -> f64 {
            if c < 0 {
                (-c - 1) as f64
            } else {
                (n as i64 + c - 1) as f64
            }
        };
        let mut out = Vec::with_capacity(self.merge.len());
        for (i, m) in self.merge.iter().enumerate() {
            let size = |c: i64, sizes: &[f64]| if c < 0 { 1.0 } else { sizes[(c - 1) as usize] };
            let s = size(m[0], &sizes) + size(m[1], &sizes);
            sizes.push(s);
            out.push([id(m[0]), id(m[1]), self.height[i], s]);
        }
        out
    }

    /// Cut the tree into `k` groups, as R's `cutree(tree, k)`. Groups are
    /// numbered from 0 in the order their first leaf appears in the input.
    pub fn cutree(&self, k: usize) -> Result<Vec<usize>> {
        let n = self.merge.len() + 1;
        if k < 1 || k > n {
            return arg_err(format!("k must be between 1 and {n}"));
        }
        // The step (1-based) whose cluster each leaf currently belongs to,
        // or 0 while it is a singleton.
        let mut step = vec![0usize; n];
        for (s, m) in self.merge.iter().take(n - k).enumerate() {
            let s = s + 1;
            for &c in m {
                if c < 0 {
                    step[(-c - 1) as usize] = s;
                } else {
                    let old = c as usize;
                    step.iter_mut().filter(|x| **x == old).for_each(|x| *x = s);
                }
            }
        }
        let mut ids = vec![usize::MAX; n];
        let mut next = 0;
        Ok(step
            .iter()
            .map(|&s| {
                if s == 0 {
                    next += 1;
                    next - 1
                } else {
                    if ids[s - 1] == usize::MAX {
                        ids[s - 1] = next;
                        next += 1;
                    }
                    ids[s - 1]
                }
            })
            .collect())
    }

    /// Cut the tree at height `h`, as R's `cutree(tree, h = h)`: merges at
    /// heights up to `h` are kept. The heights must be increasing, which is
    /// not guaranteed for the centroid and median methods.
    pub fn cutree_height(&self, h: f64) -> Result<Vec<usize>> {
        if self.height.windows(2).any(|w| w[1] < w[0]) {
            return arg_err("the tree heights are not sorted increasingly");
        }
        let n = self.merge.len() + 1;
        let first_above = self
            .height
            .iter()
            .position(|&x| x > h)
            .unwrap_or(self.height.len());
        self.cutree(n - first_above)
    }
}

#[cfg(test)]
mod tests {
    use super::format_g10 as g;

    #[test]
    fn percent_g() {
        assert_eq!(g(0.5), "0.5");
        assert_eq!(g(1.0 / 3.0), "0.3333333333");
        assert_eq!(g(123.456), "123.456");
        assert_eq!(g(0.00012345678901), "0.000123456789");
        assert_eq!(g(1.5e-7), "1.5e-07");
        assert_eq!(g(12345678901.0), "1.23456789e+10");
    }
}

/// Options for [`motif_tree`].
#[derive(Debug, Clone, PartialEq)]
pub struct TreeOptions {
    /// Fewest overlapping columns when comparing motifs.
    pub min_overlap: usize,
    /// Also try reverse complements when comparing motifs.
    pub rc: bool,
    pub method: Method,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for TreeOptions {
    fn default() -> Self {
        TreeOptions {
            min_overlap: 6,
            rc: true,
            method: Method::Complete,
            nthreads: 1,
        }
    }
}

/// Cluster DNA or RNA motifs into a tree, as `motif_tree_lite()` does:
/// motifs are compared with [`crate::compare::compare_matrix`], mean column
/// correlations become distances as `(1 - score) / 2`, and the distances are
/// clustered with [`hclust`].
pub fn motif_tree(motifs: &[crate::motif::Motif], opts: &TreeOptions) -> Result<Tree> {
    let copts = crate::compare::CompareOptions {
        min_overlap: opts.min_overlap,
        rc: opts.rc,
        nthreads: opts.nthreads,
        ..Default::default()
    };
    let scores = crate::compare::compare_matrix(motifs, &copts, false)?.scores;
    if scores.iter().flatten().any(|s| s.is_nan()) {
        return arg_err("some motifs could not be compared; check the inputs or lower min_overlap");
    }
    let dist: Vec<Vec<f64>> = scores
        .iter()
        .map(|r| r.iter().map(|s| (1.0 - s) / 2.0).collect())
        .collect();
    hclust(&dist, opts.method)
}
