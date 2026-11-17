//! Motif hits in the contribution scores a sequence model gives its inputs,
//! with [`scan_contributions`].
//!
//! After TF-MoDISco has gathered a model's contribution scores into patterns,
//! the next step is usually to find where each pattern occurs, by matching
//! its contribution weight matrix (CWM) against the scores of every input.
//! Each window is scored by two measures, as hit callers for TF-MoDISco
//! patterns usually score their hits: its similarity, the cosine between the
//! window's scores and the CWM, which says whether the contributions have the
//! shape of the pattern, and its importance, the sum of the window's absolute
//! scores, which says whether there is much contribution there at all. A
//! window is a hit when both reach the motif's thresholds.
//!
//! Every motif is scored on its own, in a single pass, so motifs that share a
//! site all have hits there, where a caller that fits all the motifs to each
//! input at once makes them compete for the same scores.
//!
//! When the model's hypothetical scores are at hand, the scores every letter
//! would have had at each position, [`scan_hypothetical`] matches them
//! instead. A window's hypothetical scores say which letters the model wants
//! there, including at the positions where a site has the wrong letter, so
//! they can be compared with the motif's preferences as a whole: with a CWM
//! of hypothetical scores, or with any probabilistic motif, such as one from
//! a database, through its log-odds. Importance still comes from the
//! observed scores.
//!
//! # References
//!
//! - Shrikumar A, Tian K, Avsec Ž, Shcherbina A, Banerjee A, Sharmin M, Nair S,
//!   Kundaje A (2018). Technical note on transcription factor motif discovery
//!   from importance scores (TF-MoDISco) version 0.5.6.5. arXiv:1811.00416.
//!   <https://doi.org/10.48550/arXiv.1811.00416>

use rayon::prelude::*;

use crate::convert::convert_type;
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
pub struct ContribScanOptions {
    /// Also match the reverse complement of each CWM.
    pub rc: bool,
    /// Match gapped CWMs with their gaps; otherwise their columns are
    /// matched as if they followed one another.
    pub use_gaps: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for ContribScanOptions {
    fn default() -> Self {
        ContribScanOptions {
            rc: true,
            use_gaps: true,
            nthreads: 1,
        }
    }
}

/// Combinations of gap sizes beyond which a gapped CWM is refused.
const MAX_GAP_COMBOS: usize = 1 << 16;

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
    /// The sum of the window's absolute scores, on the scale of
    /// [`ContribThreshold::min_importance`]: as a fraction of the sum of the
    /// CWM's absolute values, or of the importance its threshold is
    /// otherwise relative to.
    pub importance: f64,
    /// The size of each gap of a gapped CWM in this hit, in the CWM's own
    /// order, and empty otherwise.
    pub gaps: Vec<usize>,
}

/// A CWM ready for scanning: its columns scaled to unit norm, as `width * 4`
/// values position by position, for each strand searched, in pieces of
/// `widths` columns with a gap of one of `sizes` after each piece but the
/// last.
struct Filter {
    widths: Vec<usize>,
    sizes: Vec<(usize, usize)>,
    strands: Vec<(Strand, Vec<f32>)>,
    min_similarity: f64,
    min_importance: f64,
    /// The importance that `min_importance` was given as a fraction of.
    reference: f64,
}

impl Filter {
    /// Every combination of gap sizes, in the CWM's own order.
    fn combos(&self) -> Vec<Vec<usize>> {
        let mut out = vec![Vec::new()];
        for &(lo, hi) in &self.sizes {
            out = out
                .into_iter()
                .flat_map(|c| {
                    (lo..=hi).map(move |g| {
                        let mut c = c.clone();
                        c.push(g);
                        c
                    })
                })
                .collect();
        }
        out
    }
}

/// Piece widths, and the smallest and largest size of the gap after each
/// piece but the last.
type Pieces = (Vec<usize>, Vec<(usize, usize)>);

/// The pieces of a CWM between its gaps, and the smallest and largest size
/// of each gap; a CWM without gaps, or scanned without them, is one piece.
fn pieces(m: &Motif, use_gaps: bool) -> Result<Pieces> {
    let g = &m.gap;
    let w = m.width();
    if !use_gaps || !g.isgapped || g.gaploc.is_empty() {
        return Ok((vec![w], Vec::new()));
    }
    if g.gaploc.windows(2).any(|p| p[0] >= p[1]) || g.gaploc.iter().any(|&l| l == 0 || l >= w) {
        return arg_err(format!(
            "the gap locations of {} must be increasing and within the motif",
            m.name
        ));
    }
    let mut bounds = vec![0];
    bounds.extend(&g.gaploc);
    bounds.push(w);
    let sizes: Vec<(usize, usize)> = g
        .mingap
        .iter()
        .copied()
        .zip(g.maxgap.iter().copied())
        .collect();
    if sizes.len() != g.gaploc.len() || sizes.iter().any(|(lo, hi)| lo > hi) {
        return arg_err(format!(
            "the gap sizes of {} do not fit its gap locations",
            m.name
        ));
    }
    let combos: usize = sizes.iter().map(|(lo, hi)| hi - lo + 1).product();
    if combos > MAX_GAP_COMBOS {
        return arg_err(format!(
            "{} has {combos} combinations of gap sizes, more than {MAX_GAP_COMBOS} can be matched",
            m.name
        ));
    }
    Ok((bounds.windows(2).map(|b| b[1] - b[0]).collect(), sizes))
}

/// Pseudocount added to each probability of a motif matched against
/// hypothetical scores, so that letters it never saw give finite log-odds.
const HYPOTHETICAL_PSEUDOCOUNT: f64 = 0.01;

/// The values a motif is matched with, column by column, and the importance
/// its thresholds are relative to. Observed scores are matched with a CWM as
/// it is, relative to the sum of its absolute values. Hypothetical scores are
/// matched with values centred in each column, so that only the preferences
/// among letters count: a CWM's own values, relative to the sum of each
/// column's largest absolute value, the importance of a site with the
/// favoured letters, or another motif's log2 probabilities, which carry no
/// scale, relative to its width, so that the threshold is a mean per
/// position.
fn motif_values(m: &Motif, hypothetical: bool) -> Result<(Vec<f64>, f64)> {
    let is_cwm = m.motif_type == MotifType::Cwm;
    if !hypothetical {
        if !is_cwm {
            return arg_err(format!(
                "contribution scores are matched with CWMs, but {} is a {}",
                m.name, m.motif_type
            ));
        }
        let values: Vec<f64> = m.matrix.columns().flatten().copied().collect();
        let l1 = values.iter().map(|v| v.abs()).sum();
        return Ok((values, l1));
    }
    let (columns, reference): (Vec<Vec<f64>>, f64) = if is_cwm {
        let cols: Vec<Vec<f64>> = m.matrix.columns().map(<[f64]>::to_vec).collect();
        let reference = cols
            .iter()
            .map(|c| c.iter().fold(0.0f64, |a, &v| a.max(v.abs())))
            .sum();
        (cols, reference)
    } else {
        let ppm = convert_type(m, MotifType::Ppm, None, false)?;
        let cols = ppm
            .matrix
            .columns()
            .map(|c| {
                let total: f64 = c.iter().sum::<f64>() + HYPOTHETICAL_PSEUDOCOUNT * c.len() as f64;
                c.iter()
                    .map(|&v| ((v + HYPOTHETICAL_PSEUDOCOUNT) / total).log2())
                    .collect()
            })
            .collect();
        (cols, m.width() as f64)
    };
    let values = columns
        .iter()
        .flat_map(|c| {
            let mean = c.iter().sum::<f64>() / c.len() as f64;
            c.iter().map(move |&v| v - mean)
        })
        .collect();
    Ok((values, reference))
}

fn filters(
    motifs: &[Motif],
    thresholds: &[ContribThreshold],
    opts: &ContribScanOptions,
    hypothetical: bool,
) -> Result<Vec<Filter>> {
    if thresholds.len() != 1 && thresholds.len() != motifs.len() {
        return arg_err(format!(
            "give one threshold or one per motif, not {} for {} motifs",
            thresholds.len(),
            motifs.len()
        ));
    }
    motifs
        .iter()
        .enumerate()
        .map(|(i, m)| {
            if !m.alphabet.is_nucleotide() {
                return arg_err(format!("{} is not a DNA or RNA motif", m.name));
            }
            let t = thresholds[if thresholds.len() == 1 { 0 } else { i }];
            let (widths, sizes) = pieces(m, opts.use_gaps)?;
            let (values, reference) = motif_values(m, hypothetical)?;
            let mut fwd: Vec<f32> = values.iter().map(|&v| v as f32).collect();
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
            if opts.rc {
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
                widths,
                sizes,
                strands,
                min_similarity: t.min_similarity,
                min_importance: t.min_importance * reference,
                reference: if reference > 0.0 { reference } else { 1.0 },
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

/// The hits of one sequence. Windows are matched by the scores `x`, centred
/// in each position with `centre`, and their importance comes from the
/// scores `imp`, both in the layout (4, `width`).
fn scan_one(
    sequence: usize,
    x: &[f32],
    imp: &[f32],
    centre: bool,
    width: usize,
    len: usize,
    plans: &[Plan<'_>],
) -> Vec<ContribHit> {
    // The scores position by position, and running sums of each position's
    // absolute and squared scores.
    let mut xt = vec![0.0f32; 4 * len];
    for a in 0..4 {
        for (p, &v) in x[a * width..a * width + len].iter().enumerate() {
            xt[p * 4 + a] = v;
        }
    }
    if centre {
        for col in xt.as_chunks_mut::<4>().0 {
            let mean = col.iter().sum::<f32>() / 4.0;
            col.iter_mut().for_each(|v| *v -= mean);
        }
    }
    let mut abs = vec![0.0f64; len + 1];
    let mut sq = vec![0.0f64; len + 1];
    for p in 0..len {
        let col = &xt[p * 4..p * 4 + 4];
        let imp_sum: f64 = (0..4).map(|a| f64::from(imp[a * width + p].abs())).sum();
        abs[p + 1] = abs[p] + imp_sum;
        sq[p + 1] = sq[p]
            + col
                .iter()
                .map(|&v| f64::from(v) * f64::from(v))
                .sum::<f64>();
    }
    let mut hits = Vec::new();
    // The best combination of gap sizes at each start, for gapped motifs,
    // and the starts that hold one. Ungapped motifs have one combination and
    // keep every window that passes.
    let mut best: Vec<Option<ContribHit>> = Vec::new();
    let mut touched: Vec<usize> = Vec::new();
    for pl in plans {
        let single = pl.layouts.len() == 1;
        if !single && best.len() < len {
            best.resize(len, None);
        }
        for (gaps, offsets, span) in &pl.layouts {
            let span = *span;
            if span == 0 || span > len {
                continue;
            }
            for p in 0..=len - span {
                let sum = |acc: &[f64]| -> f64 {
                    offsets
                        .iter()
                        .zip(&pl.widths)
                        .map(|(&o, &w)| acc[p + o + w] - acc[p + o])
                        .sum()
                };
                let importance = sum(&abs);
                if importance < pl.min_importance {
                    continue;
                }
                let ss = sum(&sq);
                if ss <= 0.0 {
                    continue;
                }
                let score: f64 = offsets
                    .iter()
                    .zip(&pl.widths)
                    .zip(&pl.columns)
                    .map(|((&o, &w), &c)| {
                        f64::from(dot(
                            &xt[(p + o) * 4..(p + o + w) * 4],
                            &pl.filt[c * 4..(c + w) * 4],
                        ))
                    })
                    .sum();
                let similarity = score / ss.sqrt();
                if similarity < pl.min_similarity {
                    continue;
                }
                if !single && best[p].as_ref().is_some_and(|h| h.score >= score) {
                    continue;
                }
                let hit = ContribHit {
                    motif: pl.motif,
                    sequence,
                    start: p,
                    end: p + span,
                    strand: pl.strand,
                    score,
                    similarity,
                    importance: importance / pl.reference,
                    gaps: gaps.clone(),
                };
                if single {
                    hits.push(hit);
                } else {
                    if best[p].is_none() {
                        touched.push(p);
                    }
                    best[p] = Some(hit);
                }
            }
        }
        if !single {
            touched.sort_unstable();
            for p in touched.drain(..) {
                hits.extend(best[p].take());
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
/// A gapped CWM, with `opts.use_gaps`, is matched with every combination of
/// its gap sizes. The letters in a gap belong to no column of the CWM, so
/// they count towards neither the similarity nor the importance, and of the
/// combinations that reach the thresholds at a start, the one with the
/// highest score is kept, its gap sizes given in [`ContribHit::gaps`].
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
    let filters = filters(cwms, thresholds, opts, false)?;
    run(&filters, scores, scores, false, lengths, width, opts)
}

/// One strand of a motif ready for scanning: its pieces, and where they lie
/// for each combination of gap sizes, worked out once for every sequence.
struct Plan<'a> {
    motif: usize,
    strand: Strand,
    filt: &'a [f32],
    widths: Vec<usize>,
    /// The first column of each piece in `filt`.
    columns: Vec<usize>,
    /// For each combination of gap sizes, in the CWM's own order: the sizes,
    /// the offset of each piece in the window and the window's span.
    layouts: Vec<(Vec<usize>, Vec<usize>, usize)>,
    min_similarity: f64,
    min_importance: f64,
    reference: f64,
}

fn plans(filters: &[Filter]) -> Vec<Plan<'_>> {
    let mut out = Vec::new();
    for (m, f) in filters.iter().enumerate() {
        let combos = f.combos();
        for (strand, filt) in &f.strands {
            // On the minus strand the pieces and gaps come in reverse order.
            let minus = *strand == Strand::Minus;
            let widths: Vec<usize> = if minus {
                f.widths.iter().rev().copied().collect()
            } else {
                f.widths.clone()
            };
            let mut columns = vec![0usize; widths.len()];
            for k in 1..widths.len() {
                columns[k] = columns[k - 1] + widths[k - 1];
            }
            let layouts = combos
                .iter()
                .map(|gaps| {
                    let between: Vec<usize> = if minus {
                        gaps.iter().rev().copied().collect()
                    } else {
                        gaps.clone()
                    };
                    let mut offsets = vec![0usize; widths.len()];
                    for k in 1..widths.len() {
                        offsets[k] = offsets[k - 1] + widths[k - 1] + between[k - 1];
                    }
                    let span = offsets[widths.len() - 1] + widths[widths.len() - 1];
                    (gaps.clone(), offsets, span)
                })
                .collect();
            out.push(Plan {
                motif: m,
                strand: *strand,
                filt,
                widths,
                columns,
                layouts,
                min_similarity: f.min_similarity,
                min_importance: f.min_importance,
                reference: f.reference,
            });
        }
    }
    out
}

/// Scan each sequence on the thread pool and gather the hits in order.
fn run(
    filters: &[Filter],
    x: &[f32],
    imp: &[f32],
    centre: bool,
    lengths: &[usize],
    width: usize,
    opts: &ContribScanOptions,
) -> Result<Vec<ContribHit>> {
    let plans = plans(filters);
    let per_seq: Vec<Vec<ContribHit>> = with_threads(opts.nthreads, || {
        lengths
            .par_iter()
            .enumerate()
            .map(|(s, &len)| {
                let part = s * 4 * width..(s + 1) * 4 * width;
                scan_one(s, &x[part.clone()], &imp[part], centre, width, len, &plans)
            })
            .collect()
    })?;
    let mut hits: Vec<ContribHit> = per_seq.into_iter().flatten().collect();
    hits.sort_by_key(|h| (h.motif, h.sequence, h.start, h.strand == Strand::Minus));
    Ok(hits)
}

/// Find the windows of hypothetical contribution scores that match each
/// motif.
///
/// `hypothetical` holds the score each letter would have had at each
/// position, and `projected` those of the letters present, zero for the
/// others, both laid out as for [`scan_contributions`]. Each window is
/// matched by its hypothetical scores, centred in each position so that
/// only the preferences among letters count, against the motif's values
/// centred the same way, and its similarity is the cosine of the two. The
/// motifs can be CWMs, best built from hypothetical scores, or any other
/// type, whose probabilities are matched as log2 values after a small
/// pseudocount; a database motif can then be matched against a model's
/// scores directly. Importance is the sum of the window's absolute
/// projected scores, as in [`scan_contributions`].
///
/// A CWM's `min_importance` is a fraction of the importance of a site with
/// its favoured letters, the sum of each column's largest absolute value.
/// Other motifs say nothing about the size of contributions, so theirs is
/// the least mean absolute projected score per position of the window. The
/// similarities of probabilistic motifs run lower than those of CWMs, as
/// log-odds and contributions differ in shape where a motif is loose, so
/// they want a lower `min_similarity`; 0.5 suits them well.
///
/// Hits, gaps and strands are as in [`scan_contributions`].
///
/// # Example
///
/// ```
/// use universalmotif::contributions::{ContribScanOptions, ContribThreshold, scan_hypothetical};
/// use universalmotif::create::{CreateOptions, from_strings};
///
/// // GAT, as a PPM, and the hypothetical scores of CCGATCC, where the model
/// // wants G, A and T over the middle three positions.
/// let gat = from_strings(&["GAT"], None, &CreateOptions::default())?;
/// let (len, wanted) = (7, [(2, 2), (3, 0), (4, 3)]);
/// let mut hyp = vec![-0.05f32; 4 * len];
/// for (p, letter) in wanted {
///     hyp[letter * len + p] = 0.9;
/// }
/// // The observed scores are those of the letters present, C C G A T C C.
/// let seq = [1, 1, 2, 0, 3, 1, 1];
/// let mut projected = vec![0.0f32; 4 * len];
/// for (p, &letter) in seq.iter().enumerate() {
///     projected[letter * len + p] = hyp[letter * len + p];
/// }
/// let t = ContribThreshold { min_similarity: 0.5, min_importance: 0.3 };
/// let hits = scan_hypothetical(&[gat], &hyp, &projected, &[len], len, &[t], &ContribScanOptions::default())?;
/// assert_eq!((hits[0].start, hits[0].end), (2, 5));
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn scan_hypothetical(
    motifs: &[Motif],
    hypothetical: &[f32],
    projected: &[f32],
    lengths: &[usize],
    width: usize,
    thresholds: &[ContribThreshold],
    opts: &ContribScanOptions,
) -> Result<Vec<ContribHit>> {
    for (name, scores) in [("hypothetical", hypothetical), ("projected", projected)] {
        if scores.len() != lengths.len() * 4 * width {
            return arg_err(format!(
                "{} {name} scores do not fill {} sequences of 4 letters and {} positions",
                scores.len(),
                lengths.len(),
                width
            ));
        }
    }
    if let Some(&l) = lengths.iter().find(|&&l| l > width) {
        return arg_err(format!(
            "a sequence of length {l} is longer than the {width} positions of the scores"
        ));
    }
    let filters = filters(motifs, thresholds, opts, true)?;
    run(
        &filters,
        hypothetical,
        projected,
        true,
        lengths,
        width,
        opts,
    )
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
            // Importance is reported as a fraction of the CWM's absolute sum.
            let l1: f64 = m.iter().flatten().map(|v| v.abs()).sum();
            assert!((h.importance - imp / l1).abs() < 1e-5);
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
    fn gapped_cwms() {
        // GATA, a gap of 2 to 6 letters, then CACGTG.
        let core = "GATACACGTG";
        let cols: Vec<[f64; 4]> = core
            .chars()
            .map(|c| {
                let mut col = [0.0; 4];
                col["ACGT".find(c).unwrap()] = 1.0;
                col
            })
            .collect();
        let m = crate::ops::add_gap(&cwm(&cols, "gapped"), Some(&[4]), &[2], &[6]).unwrap();
        // The sites with a gap of 4, whose letters score highly but belong to
        // no column, and its reverse complement in a second sequence.
        let fwd = "TTGATATTTTCACGTGTT";
        let rev: String = fwd
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        let width = fwd.len();
        let mut vals = vec![0.01f32; width];
        vals[2..16].fill(0.5);
        let mut x = project(fwd, &vals, width);
        let rev_vals: Vec<f32> = vals.iter().rev().copied().collect();
        x.extend(project(&rev, &rev_vals, width));
        let strict = ContribThreshold {
            min_similarity: 0.99,
            min_importance: 0.0,
        };
        let opts = ContribScanOptions::default();
        let hits = scan_contributions(
            std::slice::from_ref(&m),
            &x,
            &[width, width],
            width,
            &[strict],
            &opts,
        )
        .unwrap();
        let got: Vec<(usize, usize, usize, Strand, Vec<usize>)> = hits
            .iter()
            .map(|h| (h.sequence, h.start, h.end, h.strand, h.gaps.clone()))
            .collect();
        assert_eq!(
            got,
            [
                (0, 2, 16, Strand::Plus, vec![4]),
                (1, 2, 16, Strand::Minus, vec![4])
            ]
        );
        // The gap's letters count towards neither measure.
        assert!((hits[0].similarity - 1.0).abs() < 1e-6);
        // Ten positions of 0.5, as a fraction of the ten columns of 1.
        assert!((hits[0].importance - 0.5).abs() < 1e-5);
        // Without the gaps the columns run together and match nothing.
        let flat = ContribScanOptions {
            use_gaps: false,
            ..Default::default()
        };
        let hits = scan_contributions(&[m], &x, &[width, width], width, &[strict], &flat).unwrap();
        assert!(hits.is_empty());
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
