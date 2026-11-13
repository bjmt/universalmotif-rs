//! Motif hits in the contribution scores a sequence model gives its inputs,
//! with [`scan_contributions`].
//!
//! After TF-MoDISco has gathered a model's contribution scores into patterns,
//! the next step is usually to find where each pattern occurs, by matching
//! its contribution weight matrix (CWM) against the scores of every input.
//! Each window is scored by two measures, as FiNeMo and TF-MoDISco's own hit
//! caller score their hits: its similarity, the cosine between the window's
//! scores and the CWM, which says whether the contributions have the shape of
//! the pattern, and its importance, the sum of the window's absolute scores,
//! which says whether there is much contribution there at all. A window is a
//! hit when both reach the motif's thresholds.
//!
//! Every motif is scored on its own, in a single pass, so motifs that share a
//! site all have hits there; FiNeMo instead fits all the motifs to each input
//! at once, so that they compete for the same scores.

use rayon::prelude::*;

use crate::error::{Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

/// How a motif's hits are called: windows need at least this similarity and
/// importance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContribThreshold {
    /// The least cosine between a window's scores and the CWM, from -1 to 1.
    pub min_similarity: f64,
    /// The least importance of a window, the sum of its absolute scores, as a
    /// fraction of the sum of the CWM's absolute values. A CWM from
    /// TF-MoDISco is the average of its seqlets, so 1 is about the importance
    /// of a typical seqlet.
    pub min_importance: f64,
}

impl Default for ContribThreshold {
    fn default() -> Self {
        ContribThreshold {
            min_similarity: 0.7,
            min_importance: 0.5,
        }
    }
}

/// Options for [`scan_contributions`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ContribScanOptions {
    /// Also match the reverse complement of each CWM.
    pub rc: bool,
    /// Threads to use, at least 1.
    pub nthreads: usize,
}

impl Default for ContribScanOptions {
    fn default() -> Self {
        ContribScanOptions {
            rc: true,
            nthreads: 1,
        }
    }
}

/// A window of a sequence's scores matching a CWM.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ContribHit {
    /// The index of the motif.
    pub motif: usize,
    /// The index of the sequence.
    pub sequence: usize,
    /// The first position of the window, counted from 0.
    pub start: usize,
    /// The position after the last of the window.
    pub end: usize,
    /// `Minus` when the window matches the CWM's reverse complement.
    pub strand: Strand,
    /// The dot product of the window's scores with the CWM scaled to unit
    /// norm, which is the similarity times the Euclidean norm of the window,
    /// so that strong and well-shaped windows score highest.
    pub score: f64,
    /// The cosine between the window's scores and the CWM.
    pub similarity: f64,
    /// The sum of the window's absolute scores.
    pub importance: f64,
}

/// A CWM ready for scanning: its columns scaled to unit norm, as `width * 4`
/// values position by position, for each strand searched.
struct Filter {
    width: usize,
    strands: Vec<(Strand, Vec<f32>)>,
    min_similarity: f64,
    min_importance: f64,
}

fn filters(cwms: &[Motif], thresholds: &[ContribThreshold], rc: bool) -> Result<Vec<Filter>> {
    if thresholds.len() != 1 && thresholds.len() != cwms.len() {
        return arg_err(format!(
            "give one threshold or one per motif, not {} for {} motifs",
            thresholds.len(),
            cwms.len()
        ));
    }
    cwms.iter()
        .enumerate()
        .map(|(i, m)| {
            if m.motif_type != MotifType::Cwm {
                return arg_err(format!(
                    "contribution scores are matched with CWMs, but {} is a {}",
                    m.name, m.motif_type
                ));
            }
            if !m.alphabet.is_nucleotide() {
                return arg_err(format!("{} is not a DNA or RNA motif", m.name));
            }
            let t = thresholds[if thresholds.len() == 1 { 0 } else { i }];
            let width = m.width();
            let mut fwd = Vec::with_capacity(width * 4);
            for col in m.matrix.columns() {
                fwd.extend(col.iter().map(|&v| v as f32));
            }
            let l1: f64 = m.matrix.columns().flatten().map(|v| v.abs()).sum();
            let norm = fwd
                .iter()
                .map(|&v| f64::from(v) * f64::from(v))
                .sum::<f64>()
                .sqrt();
            if norm > 0.0 {
                fwd.iter_mut()
                    .for_each(|v| *v = (f64::from(*v) / norm) as f32);
            }
            let mut strands = Vec::with_capacity(2);
            if rc {
                // Columns in reverse order, each with its letters complemented.
                let rev: Vec<f32> = fwd
                    .as_chunks::<4>()
                    .0
                    .iter()
                    .rev()
                    .flat_map(|c| c.iter().rev().copied())
                    .collect();
                strands.push((Strand::Plus, fwd));
                strands.push((Strand::Minus, rev));
            } else {
                strands.push((Strand::Plus, fwd));
            }
            Ok(Filter {
                width,
                strands,
                min_similarity: t.min_similarity,
                min_importance: t.min_importance * l1,
            })
        })
        .collect()
}

/// The dot product of two slices of equal length, summed in eight lanes so
/// that it vectorizes.
fn dot(a: &[f32], b: &[f32]) -> f32 {
    let mut lanes = [0.0f32; 8];
    let ((ca, ra), (cb, rb)) = (a.as_chunks::<8>(), b.as_chunks::<8>());
    let tail: f32 = ra.iter().zip(rb).map(|(x, y)| x * y).sum();
    for (x, y) in ca.iter().zip(cb) {
        for k in 0..8 {
            lanes[k] += x[k] * y[k];
        }
    }
    lanes.iter().sum::<f32>() + tail
}

fn scan_one(
    sequence: usize,
    x: &[f32],
    width: usize,
    len: usize,
    filters: &[Filter],
) -> Vec<ContribHit> {
    // The scores position by position, and running sums of each position's
    // absolute and squared scores.
    let mut xt = vec![0.0f32; 4 * len];
    for a in 0..4 {
        for (p, &v) in x[a * width..a * width + len].iter().enumerate() {
            xt[p * 4 + a] = v;
        }
    }
    let mut abs = vec![0.0f64; len + 1];
    let mut sq = vec![0.0f64; len + 1];
    for p in 0..len {
        let col = &xt[p * 4..p * 4 + 4];
        abs[p + 1] = abs[p] + col.iter().map(|&v| f64::from(v.abs())).sum::<f64>();
        sq[p + 1] = sq[p]
            + col
                .iter()
                .map(|&v| f64::from(v) * f64::from(v))
                .sum::<f64>();
    }
    let mut hits = Vec::new();
    for (m, f) in filters.iter().enumerate() {
        let w = f.width;
        if w == 0 || w > len {
            continue;
        }
        for p in 0..=len - w {
            let importance = abs[p + w] - abs[p];
            if importance < f.min_importance {
                continue;
            }
            let ss = sq[p + w] - sq[p];
            if ss <= 0.0 {
                continue;
            }
            let norm = ss.sqrt();
            let window = &xt[p * 4..(p + w) * 4];
            for (strand, filt) in &f.strands {
                let score = f64::from(dot(window, filt));
                let similarity = score / norm;
                if similarity >= f.min_similarity {
                    hits.push(ContribHit {
                        motif: m,
                        sequence,
                        start: p,
                        end: p + w,
                        strand: *strand,
                        score,
                        similarity,
                        importance,
                    });
                }
            }
        }
    }
    hits
}

/// Find the windows of contribution scores that match each CWM.
///
/// `scores` holds the projected scores of `lengths.len()` sequences, those of
/// the letter present at each position and zero for the others, in the layout
/// (sequences, 4, `width`) with the letters in the order ACGT; sequence `i`
/// uses its first `lengths[i]` positions. Every window of every sequence is
/// compared with every CWM, and with its reverse complement with `opts.rc`,
/// and is a hit when its similarity and importance reach the motif's
/// threshold, given once for every motif or once per motif. A window that
/// matches both strands is a hit on both. Hits come in the order of the
/// motifs, then of the sequences, then of the starts, with `+` before `-`.
///
/// # Example
///
/// ```
/// use universalmotif::Alphabet;
/// use universalmotif::contributions::{ContribScanOptions, ContribThreshold, scan_contributions};
/// use universalmotif::create::{CreateOptions, from_typed_matrix};
/// use universalmotif::matrix::Matrix;
/// use universalmotif::motif::MotifType;
///
/// let opts = CreateOptions { name: "GAT".into(), motif_type: MotifType::Cwm, ..Default::default() };
/// let cols = [[0.0, 0.0, 0.9, 0.0], [0.8, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.7]];
/// let gat = from_typed_matrix(Matrix::from_columns(&cols)?, None, Some(Alphabet::Dna), MotifType::Cwm, &opts)?;
/// // Scores of CCGATCC, high over GAT, in the layout (sequences, 4, length).
/// let mut scores = vec![0.0f32; 4 * 7];
/// for (p, letter, v) in [(0, 1, 0.1), (1, 1, 0.1), (2, 2, 0.9), (3, 0, 0.8), (4, 3, 0.6), (5, 1, 0.1), (6, 1, 0.0)] {
///     scores[letter * 7 + p] = v;
/// }
/// // Windows that share two letters with the site look alike too, so keep to
/// // close matches.
/// let strict = ContribThreshold { min_similarity: 0.9, ..Default::default() };
/// let hits = scan_contributions(&[gat], &scores, &[7], 7, &[strict], &ContribScanOptions::default())?;
/// assert_eq!(hits.len(), 1);
/// assert_eq!((hits[0].start, hits[0].end), (2, 5));
/// assert!(hits[0].similarity > 0.99);
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn scan_contributions(
    cwms: &[Motif],
    scores: &[f32],
    lengths: &[usize],
    width: usize,
    thresholds: &[ContribThreshold],
    opts: &ContribScanOptions,
) -> Result<Vec<ContribHit>> {
    if scores.len() != lengths.len() * 4 * width {
        return arg_err(format!(
            "{} scores do not fill {} sequences of 4 letters and {} positions",
            scores.len(),
            lengths.len(),
            width
        ));
    }
    if let Some(&l) = lengths.iter().find(|&&l| l > width) {
        return arg_err(format!(
            "a sequence of length {l} is longer than the {width} positions of the scores"
        ));
    }
    let filters = filters(cwms, thresholds, opts.rc)?;
    let per_seq: Vec<Vec<ContribHit>> = with_threads(opts.nthreads, || {
        lengths
            .par_iter()
            .enumerate()
            .map(|(s, &len)| {
                scan_one(
                    s,
                    &scores[s * 4 * width..(s + 1) * 4 * width],
                    width,
                    len,
                    &filters,
                )
            })
            .collect()
    })?;
    let mut hits: Vec<ContribHit> = per_seq.into_iter().flatten().collect();
    hits.sort_by_key(|h| (h.motif, h.sequence, h.start, h.strand == Strand::Minus));
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Alphabet;
    use crate::create::{CreateOptions, from_typed_matrix};
    use crate::matrix::Matrix;

    fn cwm(cols: &[[f64; 4]], name: &str) -> Motif {
        let opts = CreateOptions {
            name: name.into(),
            motif_type: MotifType::Cwm,
            ..Default::default()
        };
        from_typed_matrix(
            Matrix::from_columns(cols).unwrap(),
            None,
            Some(Alphabet::Dna),
            MotifType::Cwm,
            &opts,
        )
        .unwrap()
    }

    /// Projected scores of a sequence with `values` at each position.
    fn project(seq: &str, values: &[f32], width: usize) -> Vec<f32> {
        let mut out = vec![0.0f32; 4 * width];
        for (p, (c, &v)) in seq.chars().zip(values).enumerate() {
            let a = "ACGT".find(c).unwrap();
            out[a * width + p] = v;
        }
        out
    }

    /// The best similarity and importance by brute force, with numpy's
    /// definitions written out.
    fn brute(
        x: &[f32],
        width: usize,
        len: usize,
        m: &[[f64; 4]],
        start: usize,
        rc: bool,
    ) -> (f64, f64) {
        let w = m.len();
        let col = |j: usize| -> [f64; 4] {
            if rc {
                let c = m[w - 1 - j];
                [c[3], c[2], c[1], c[0]]
            } else {
                m[j]
            }
        };
        let (mut num, mut nx, mut nm, mut imp) = (0.0, 0.0, 0.0, 0.0);
        for j in 0..w {
            let c = col(j);
            for a in 0..4 {
                let v = f64::from(x[a * width + start + j]);
                num += v * c[a];
                nx += v * v;
                nm += c[a] * c[a];
                imp += v.abs();
            }
        }
        assert!(start + w <= len);
        (num / (nx.sqrt() * nm.sqrt()), imp)
    }

    #[test]
    fn matches_brute_force() {
        let m = [
            [0.0, 0.0, 0.9, 0.0],
            [0.8, 0.0, -0.1, 0.0],
            [0.0, 0.0, 0.0, 0.7],
            [0.3, 0.1, 0.0, 0.0],
        ];
        let seq = "CCGATACCTATCCGG";
        let vals: Vec<f32> = seq
            .chars()
            .enumerate()
            .map(|(i, _)| ((i * 7) % 5) as f32 * 0.2 - 0.1)
            .collect();
        let width = 18;
        let x = project(seq, &vals, width);
        let all = ContribThreshold {
            min_similarity: -1.0,
            min_importance: 0.0,
        };
        let hits = scan_contributions(
            &[cwm(&m, "m")],
            &x,
            &[seq.len()],
            width,
            &[all],
            &ContribScanOptions::default(),
        )
        .unwrap();
        // Every window on both strands, except those of zero scores.
        assert!(hits.len() >= 2 * (seq.len() - 3) - 4);
        for h in &hits {
            let (sim, imp) = brute(&x, width, seq.len(), &m, h.start, h.strand == Strand::Minus);
            assert!((h.similarity - sim).abs() < 1e-5, "{h:?} {sim}");
            assert!((h.importance - imp).abs() < 1e-5);
            assert!(h.end <= seq.len());
        }
    }

    #[test]
    fn thresholds_and_strands() {
        let m = [
            [0.0, 0.0, 0.9, 0.0],
            [0.8, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.7],
        ];
        let width = 12;
        // GAT at 1 and its reverse complement ATC at 7, with weaker scores.
        let seq = "TGATGGCATCCC";
        let mut vals = vec![0.01f32; seq.len()];
        vals[1..4].copy_from_slice(&[0.9, 0.8, 0.7]);
        vals[7..10].copy_from_slice(&[0.3, 0.35, 0.4]);
        let x = project(seq, &vals, width);
        let motifs = [cwm(&m, "gat")];
        let opts = ContribScanOptions::default();
        // Windows overlapping a site by two letters still look alike, so a
        // high similarity keeps to whole sites.
        let strict = ContribThreshold {
            min_similarity: 0.9,
            ..Default::default()
        };
        let hits = scan_contributions(&motifs, &x, &[seq.len()], width, &[strict], &opts).unwrap();
        let got: Vec<(usize, Strand)> = hits.iter().map(|h| (h.start, h.strand)).collect();
        // The second site carries less than half the CWM's contribution.
        assert_eq!(got, [(1, Strand::Plus)]);
        let low = ContribThreshold {
            min_importance: 0.4,
            ..strict
        };
        let hits = scan_contributions(&motifs, &x, &[seq.len()], width, &[low], &opts).unwrap();
        let got: Vec<(usize, Strand)> = hits.iter().map(|h| (h.start, h.strand)).collect();
        assert_eq!(got, [(1, Strand::Plus), (7, Strand::Minus)]);
        let fwd_only = ContribScanOptions {
            rc: false,
            ..Default::default()
        };
        let hits = scan_contributions(&motifs, &x, &[seq.len()], width, &[low], &fwd_only).unwrap();
        assert_eq!(hits.len(), 1);
        // A shorter sequence leaves the second site out.
        let hits = scan_contributions(&motifs, &x, &[9], width, &[low], &opts).unwrap();
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn checks() {
        let ppm = crate::create::from_strings(&["GAT"], None, &CreateOptions::default()).unwrap();
        let x = vec![0.0f32; 4 * 5];
        let t = [ContribThreshold::default()];
        let opts = ContribScanOptions::default();
        assert!(scan_contributions(&[ppm], &x, &[5], 5, &t, &opts).is_err());
        let c = cwm(&[[0.0, 0.0, 0.9, 0.0]], "c");
        assert!(scan_contributions(std::slice::from_ref(&c), &x, &[5, 5], 5, &t, &opts).is_err());
        assert!(scan_contributions(std::slice::from_ref(&c), &x, &[6], 5, &t, &opts).is_err());
        assert!(scan_contributions(&[c.clone(), c], &x, &[5], 5, &[t[0]; 3], &opts).is_err());
    }
}
