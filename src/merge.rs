//! Merging motifs (`merge_motifs_lite()` and `merge_similar_lite()`).

use crate::alphabet::Alphabet;
use crate::compare::{CompareOptions, NullModel, align, compare_matrix};
use crate::convert::convert_type;
use crate::create::{CreateOptions, from_matrix};
use crate::error::{Result, arg_err};
use crate::matrix::{Matrix, which_max};
use crate::motif::{Motif, MotifType, Strand};

/// Options for [`merge_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct MergeOptions {
    /// Fewest overlapping columns for a motif to be aligned to the anchor;
    /// motifs that cannot be aligned are left out.
    pub min_overlap: usize,
    /// Allow motifs to be reverse complemented to fit the anchor.
    pub rc: bool,
    /// Name of the merged motif; by default `merged_` followed by the names
    /// of the motifs used, joined by `+`.
    pub name: Option<String>,
    /// Weight each motif by its `nsites` when averaging.
    pub weighted: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for MergeOptions {
    fn default() -> Self {
        MergeOptions {
            min_overlap: 5,
            rc: true,
            name: None,
            weighted: false,
            nthreads: 1,
        }
    }
}

fn shared_alphabet(motifs: &[Motif]) -> Result<Alphabet> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs given");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    Ok(first.alphabet.clone())
}

/// Merge motifs into one.
///
/// The motif with the highest total IC is the anchor. Every other motif is
/// aligned to it (see [`crate::compare::align`]), reverse complemented if
/// that fits better (DNA and RNA only), and the aligned probabilities are
/// averaged position by position. Positions covered by no motif are uniform. The result is a PPM
/// whose `nsites` is the sum of the motifs' (1 for motifs without one), with
/// the average background and the anchor's pseudocount.
pub fn merge_motifs(motifs: &[Motif], opts: &MergeOptions) -> Result<Motif> {
    let alphabet = shared_alphabet(motifs)?;
    let nl = alphabet.len();
    if motifs.len() == 1 {
        let mut p = convert_type(&motifs[0], MotifType::Ppm, None, false)?.into_parts();
        if let Some(n) = &opts.name {
            p.name = n.clone();
        }
        return Motif::from_parts(p);
    }
    let ppms: Vec<Motif> = motifs
        .iter()
        .map(|m| convert_type(m, MotifType::Ppm, None, false))
        .collect::<Result<_>>()?;
    let nsites: Vec<f64> = ppms
        .iter()
        .map(|m| m.nsites.filter(|n| *n > 0.0).unwrap_or(1.0))
        .collect();
    let ic: Vec<f64> = ppms.iter().map(|m| m.icscore).collect();
    let anchor = which_max(&ic);

    let pairs: Vec<(usize, usize)> = (0..ppms.len()).map(|i| (anchor, i)).collect();
    let copts = CompareOptions {
        min_overlap: opts.min_overlap,
        rc: opts.rc,
        nthreads: opts.nthreads,
        ..Default::default()
    };
    let al = align(&ppms, &pairs, &copts)?;
    let w_anchor = ppms[anchor].width() as i64;

    let mut used: Vec<(usize, Matrix, i64)> = Vec::new();
    for (i, m) in ppms.iter().enumerate() {
        if i == anchor {
            used.push((i, m.matrix.clone(), 0));
            continue;
        }
        let a = &al[i];
        let len = a.overlap as i64;
        if a.overlap < opts.min_overlap {
            continue;
        }
        let t_start = (-a.offset).max(0);
        let w = m.width() as i64;
        if a.strand == Strand::Minus {
            let q_start = w_anchor - a.query_start as i64 - len;
            let t_start_oriented = w - t_start - len;
            used.push((i, m.matrix.reversed(), q_start - t_start_oriented));
        } else {
            used.push((i, m.matrix.clone(), a.query_start as i64 - t_start));
        }
    }

    let lo = used.iter().map(|u| u.2).min().unwrap();
    let hi = used.iter().map(|u| u.2 + u.1.ncol() as i64).max().unwrap();
    let width = (hi - lo) as usize;
    let weights: Vec<f64> = used
        .iter()
        .map(|u| if opts.weighted { nsites[u.0] } else { 1.0 })
        .collect();
    let mut cols = Vec::with_capacity(width);
    for p in 0..width as i64 {
        let frame = lo + p;
        let mut sum = vec![0.0; nl];
        let mut wsum = 0.0;
        for (u, &w) in used.iter().zip(&weights) {
            let local = frame - u.2;
            if local < 0 || local >= u.1.ncol() as i64 {
                continue;
            }
            for (s, v) in sum.iter_mut().zip(u.1.col(local as usize)) {
                *s += w * v;
            }
            wsum += w;
        }
        cols.push(if wsum == 0.0 {
            vec![1.0 / nl as f64; nl]
        } else {
            sum.iter().map(|s| s / wsum).collect()
        });
    }

    let total_nsites: f64 = used.iter().map(|u| nsites[u.0]).sum();
    let mut bkg = vec![0.0; nl];
    for u in &used {
        for (b, v) in bkg.iter_mut().zip(&ppms[u.0].bkg) {
            *b += v;
        }
    }
    let bkg: Vec<f64> = bkg.iter().map(|b| b / used.len() as f64).collect();
    let name = opts.name.clone().unwrap_or_else(|| {
        let names: Vec<&str> = used.iter().map(|u| ppms[u.0].name.as_str()).collect();
        format!("merged_{}", names.join("+"))
    });
    let copts = CreateOptions {
        motif_type: MotifType::Ppm,
        name,
        nsites: Some(total_nsites),
        bkg: Some(bkg),
        pseudocount: ppms[anchor].pseudocount,
        ..Default::default()
    };
    from_matrix(Matrix::from_columns(&cols)?, None, Some(alphabet), &copts)
}

/// Options for [`merge_similar`].
#[derive(Debug, Clone, PartialEq)]
pub struct MergeSimilarOptions {
    /// Motifs are linked when their comparison q-value (in either direction)
    /// is at most this.
    pub qvalue: f64,
    pub null: NullModel,
    pub min_overlap: usize,
    pub rc: bool,
    /// Background (A, C, G, T) for the parametric null (DNA and RNA only).
    pub bkg: [f64; 4],
    pub weighted: bool,
    pub nthreads: usize,
}

impl Default for MergeSimilarOptions {
    fn default() -> Self {
        MergeSimilarOptions {
            qvalue: 0.01,
            null: NullModel::Empirical,
            min_overlap: 5,
            rc: true,
            bkg: [0.25; 4],
            weighted: false,
            nthreads: 1,
        }
    }
}

/// Cluster similar motifs and merge each cluster.
///
/// Motifs are linked when their comparison q-value is within
/// `opts.qvalue` in either direction, and clusters are the connected groups.
/// Returns the cluster of each input motif (numbered by first appearance)
/// and one motif per cluster: single motifs as they were, larger clusters
/// merged with [`merge_motifs`].
#[allow(clippy::needless_range_loop)]
pub fn merge_similar(
    motifs: &[Motif],
    opts: &MergeSimilarOptions,
) -> Result<(Vec<usize>, Vec<Motif>)> {
    if motifs.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    shared_alphabet(motifs)?;
    if !(opts.qvalue > 0.0 && opts.qvalue <= 1.0) {
        return arg_err("qvalue must be in (0, 1]");
    }
    let n = motifs.len();
    if n == 1 {
        return Ok((vec![0], motifs.to_vec()));
    }
    let copts = CompareOptions {
        min_overlap: opts.min_overlap,
        rc: opts.rc,
        null: opts.null,
        bkg: opts.bkg,
        nthreads: opts.nthreads,
    };
    let q = compare_matrix(motifs, &copts, true)?.qvalues;

    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        let mut root = x;
        while parent[root] != root {
            root = parent[root];
        }
        while parent[x] != root {
            let next = parent[x];
            parent[x] = root;
            x = next;
        }
        root
    }
    for i in 0..n - 1 {
        for j in i + 1..n {
            let qij = q[i][j].min(q[j][i]);
            if !qij.is_nan() && qij <= opts.qvalue {
                let (ra, rb) = (find(&mut parent, i), find(&mut parent, j));
                if ra != rb {
                    parent[rb] = ra;
                }
            }
        }
    }
    let mut seen: Vec<usize> = Vec::new();
    let clusters: Vec<usize> = (0..n)
        .map(|i| {
            let r = find(&mut parent, i);
            match seen.iter().position(|&s| s == r) {
                Some(p) => p,
                None => {
                    seen.push(r);
                    seen.len() - 1
                }
            }
        })
        .collect();

    let mopts = MergeOptions {
        min_overlap: opts.min_overlap,
        rc: opts.rc,
        weighted: opts.weighted,
        nthreads: opts.nthreads,
        name: None,
    };
    let merged = (0..seen.len())
        .map(|c| {
            let members: Vec<Motif> = (0..n)
                .filter(|&i| clusters[i] == c)
                .map(|i| motifs[i].clone())
                .collect();
            if members.len() == 1 {
                Ok(members.into_iter().next().unwrap())
            } else {
                merge_motifs(&members, &mopts)
            }
        })
        .collect::<Result<_>>()?;
    Ok((clusters, merged))
}
