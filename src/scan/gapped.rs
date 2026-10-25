//! Scanning with gapped motifs.
//!
//! A gapped motif is pieces of columns with a gap of variable size after
//! each piece but the last, so each combination of gap sizes amounts to a
//! separate motif with blank columns in the gaps. Rather than scanning each
//! of those, every piece is scored once per position, and the combinations
//! are tried by
//! adding the pieces' scores, piece by piece, dropping a partial combination
//! as soon as the pieces still to come can no longer reach the threshold.
//! The scores, and so the thresholds and P-values, are those of the motif
//! without its gaps, since blank columns score zero.

use super::{CodedScorer, RawHit, Scorer};
use crate::error::{Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{Motif, Strand};

/// Where a motif's gaps are and how large they can be, in the motif's own
/// orientation. Combinations of gap sizes are numbered with the first gap
/// varying fastest.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GapShape {
    /// The widths of the pieces between gaps.
    widths: Vec<usize>,
    /// The smallest and largest size of each gap.
    sizes: Vec<(usize, usize)>,
}

impl GapShape {
    /// The gaps of a gapped motif.
    pub(crate) fn from_motif(m: &Motif) -> Result<GapShape> {
        let g = &m.gap;
        let w = m.width();
        if g.gaploc.windows(2).any(|p| p[0] >= p[1]) {
            return arg_err(format!(
                "the gap locations of motif '{}' must be in increasing order",
                m.name
            ));
        }
        let mut bounds = vec![0];
        bounds.extend(&g.gaploc);
        bounds.push(w);
        let shape = GapShape {
            widths: bounds.windows(2).map(|b| b[1] - b[0]).collect(),
            sizes: g
                .mingap
                .iter()
                .copied()
                .zip(g.maxgap.iter().copied())
                .collect(),
        };
        if shape.combos() > u16::MAX as usize + 1 {
            return arg_err(format!(
                "motif '{}' has {} combinations of gap sizes, more than {} can be scanned",
                m.name,
                shape.combos(),
                u16::MAX as usize + 1
            ));
        }
        Ok(shape)
    }

    /// The widths of the pieces between gaps.
    #[cfg(test)]
    pub(crate) fn piece_widths(&self) -> &[usize] {
        &self.widths
    }

    /// The number of combinations of gap sizes.
    pub(crate) fn combos(&self) -> usize {
        self.sizes.iter().map(|(lo, hi)| hi - lo + 1).product()
    }

    /// The step in the combination number for each size of each gap.
    fn strides(&self) -> Vec<u32> {
        let mut s = 1;
        self.sizes
            .iter()
            .map(|(lo, hi)| {
                let this = s;
                s *= (hi - lo + 1) as u32;
                this
            })
            .collect()
    }

    /// The size of each gap in combination `combo`.
    pub(crate) fn gap_sizes(&self, combo: u16) -> Vec<usize> {
        let mut c = combo as usize;
        self.sizes
            .iter()
            .map(|(lo, hi)| {
                let n = hi - lo + 1;
                let g = lo + c % n;
                c /= n;
                g
            })
            .collect()
    }

    /// The full width of the motif with the gaps of combination `combo`.
    pub(crate) fn width(&self, combo: u16) -> usize {
        self.widths.iter().sum::<usize>() + self.gap_sizes(combo).iter().sum::<usize>()
    }

    /// The widest the motif can be.
    pub(crate) fn max_width(&self) -> usize {
        self.widths.iter().sum::<usize>() + self.sizes.iter().map(|s| s.1).sum::<usize>()
    }

    /// Replace the letters in the gaps of combination `combo` with dots;
    /// `letters` is the match in the motif's own orientation.
    pub(crate) fn mark_gaps(&self, combo: u16, letters: &mut [u8]) {
        let mut at = 0;
        for (w, g) in self.widths.iter().zip(self.gap_sizes(combo)) {
            at += w;
            letters[at..at + g].fill(b'.');
            at += g;
        }
    }
}

/// How one piece is scored: on packed letters for DNA and RNA, or on letter
/// codes for other alphabets.
pub(super) enum PieceScorer {
    Packed(Scorer),
    Coded(CodedScorer),
}

/// A gapped motif laid out for scanning one strand.
struct Layout {
    strand: Strand,
    pieces: Vec<PieceScorer>,
    widths: Vec<usize>,
    /// Smallest and largest size of each gap, in scanning order.
    sizes: Vec<(usize, usize)>,
    /// The combination number step of each gap, in scanning order.
    strides: Vec<u32>,
    /// The best total score of the pieces from each one on.
    best_rest: Vec<i32>,
    /// The narrowest the motif can be from each piece on.
    min_rest: Vec<usize>,
}

impl Layout {
    fn new(
        strand: Strand,
        pieces: Vec<Matrix>,
        sizes: Vec<(usize, usize)>,
        strides: Vec<u32>,
        packed: bool,
    ) -> Layout {
        let widths: Vec<usize> = pieces.iter().map(Matrix::ncol).collect();
        let best: Vec<i32> = pieces
            .iter()
            .map(|p| {
                p.columns()
                    .map(|c| c.iter().map(|&v| super::grid_value(v)).max().unwrap())
                    .sum()
            })
            .collect();
        let n = pieces.len();
        let mut best_rest = vec![0; n + 1];
        let mut min_rest = vec![0; n + 1];
        for j in (0..n).rev() {
            best_rest[j] = best_rest[j + 1] + best[j];
            min_rest[j] = min_rest[j + 1] + widths[j] + if j + 1 < n { sizes[j].0 } else { 0 };
        }
        let pieces = pieces
            .iter()
            .map(|p| {
                if packed {
                    PieceScorer::Packed(Scorer::new(p))
                } else {
                    PieceScorer::Coded(CodedScorer::new(p))
                }
            })
            .collect();
        Layout {
            strand,
            pieces,
            widths,
            sizes,
            strides,
            best_rest,
            min_rest,
        }
    }
}

/// A gapped motif ready for scanning, with a layout per strand.
pub(super) struct GappedScorer {
    pub(super) shape: GapShape,
    layouts: Vec<Layout>,
}

impl GappedScorer {
    /// Lay out `pwm`, the motif's scoring matrix without gaps, with the gaps of
    /// `shape`. With `rc`, the reverse strand of DNA and RNA uses the reverse
    /// complement of each piece, in reverse order.
    pub(super) fn new(pwm: &Matrix, shape: GapShape, nucleotide: bool, rc: bool) -> GappedScorer {
        let mut pieces = Vec::new();
        let mut at = 0;
        for &w in &shape.widths {
            pieces.push(pwm.slice_columns(at..at + w));
            at += w;
        }
        let strides = shape.strides();
        let mut layouts = vec![Layout::new(
            Strand::Plus,
            pieces.clone(),
            shape.sizes.clone(),
            strides.clone(),
            nucleotide,
        )];
        if nucleotide && rc {
            layouts.push(Layout::new(
                Strand::Minus,
                pieces.iter().rev().map(Matrix::reversed).collect(),
                shape.sizes.iter().rev().copied().collect(),
                strides.iter().rev().copied().collect(),
                true,
            ));
        }
        GappedScorer { shape, layouts }
    }

    /// Find the hits of motif `mi` starting in `starts` of a stretch of
    /// sequence `si` free of unknown letters and ending at `seg_end`.
    /// `score(piece, p)` scores a piece at position `p`, for pieces ending by
    /// `stop`, the end of the letters decoded for the chunk; `buf` holds the
    /// pieces' scores between calls.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn scan(
        &self,
        mi: usize,
        si: usize,
        starts: (usize, usize),
        stop: usize,
        seg_end: usize,
        min_score: i32,
        buf: &mut Vec<Vec<i32>>,
        score: impl Fn(&PieceScorer, usize) -> i32,
        out: &mut Vec<RawHit>,
    ) {
        let (c, end) = starts;
        for lay in &self.layouts {
            let n = lay.pieces.len();
            if buf.len() < n {
                buf.resize_with(n, Vec::new);
            }
            // Each piece's score wherever it fits.
            for (j, piece) in lay.pieces.iter().enumerate() {
                // A piece must end by `stop`, which is never past seg_end.
                let upto = (stop + 1).saturating_sub(lay.widths[j]).max(c);
                let b = &mut buf[j];
                b.clear();
                b.extend((c..upto).map(|p| score(piece, p)));
            }
            // The best that the next piece scores over every size of gap
            // `j`, from each position of piece `j`, which bounds the hits
            // starting there more tightly than its best possible score.
            let mut reach: Vec<Vec<i32>> = Vec::with_capacity(n.saturating_sub(1));
            for j in 0..n.saturating_sub(1) {
                let (lo, hi) = lay.sizes[j];
                let next = &buf[j + 1];
                let mut r = vec![i32::MIN; buf[j].len()];
                for g in lo..=hi {
                    let shift = lay.widths[j] + g;
                    if shift >= next.len() {
                        break;
                    }
                    let len = r.len().min(next.len() - shift);
                    for (x, &y) in r[..len].iter_mut().zip(&next[shift..shift + len]) {
                        *x = (*x).max(y);
                    }
                }
                reach.push(r);
            }
            let minus = lay.strand == Strand::Minus;
            let mut find = Finder {
                lay,
                buf,
                reach: &reach,
                c,
                seg_end,
                min_score,
                mi,
                si,
                minus,
                out: &mut *out,
            };
            for i in c..end {
                if i + lay.min_rest[0] > seg_end {
                    break;
                }
                let s = find.buf[0][i - c];
                if find.can_reach(0, i, s) {
                    find.extend(i, 0, i, s, 0);
                }
            }
        }
    }
}

/// The search for one strand's combinations from one window start.
struct Finder<'a> {
    lay: &'a Layout,
    buf: &'a [Vec<i32>],
    reach: &'a [Vec<i32>],
    c: usize,
    seg_end: usize,
    min_score: i32,
    mi: usize,
    si: usize,
    minus: bool,
    out: &'a mut Vec<RawHit>,
}

impl Finder<'_> {
    /// Whether pieces up to `j`, with piece `j` at `pos` and scoring
    /// `partial`, can still reach the threshold.
    fn can_reach(&self, j: usize, pos: usize, partial: i32) -> bool {
        let lay = self.lay;
        match self.reach.get(j) {
            Some(r) => {
                let next = r[pos - self.c];
                next != i32::MIN
                    && i64::from(partial) + i64::from(next) + i64::from(lay.best_rest[j + 2])
                        >= i64::from(self.min_score)
            }
            None => partial >= self.min_score,
        }
    }

    /// Piece `j` sits at `pos`, and the pieces so far score `partial`;
    /// try every size of the next gap.
    fn extend(&mut self, start: usize, j: usize, pos: usize, partial: i32, combo: u32) {
        let lay = self.lay;
        if j + 1 == lay.pieces.len() {
            if partial >= self.min_score {
                self.out.push(RawHit {
                    motif: self.mi as u32,
                    sequence: self.si as u32,
                    start: start as u32,
                    score: partial,
                    id: 0,
                    minus: self.minus,
                    combo: combo as u16,
                    pvalue: f64::NAN,
                });
            }
            return;
        }
        if !self.can_reach(j, pos, partial) {
            return;
        }
        let (lo, hi) = lay.sizes[j];
        for g in lo..=hi {
            let next = pos + lay.widths[j] + g;
            if next + lay.min_rest[j + 1] > self.seg_end {
                break;
            }
            let s = partial + self.buf[j + 1][next - self.c];
            self.extend(
                start,
                j + 1,
                next,
                s,
                combo + (g - lo) as u32 * lay.strides[j],
            );
        }
    }
}
