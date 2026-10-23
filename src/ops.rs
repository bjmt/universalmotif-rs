//! Operations that transform a motif: reverse complement, trimming,
//! rounding and renormalization.

use crate::alphabet::Alphabet;
use crate::convert::{convert_type, position_icscore};
use crate::error::{Error, Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{GapInfo, Motif, MotifType};

/// Reverse complement a DNA or RNA motif. Other alphabets are refused unless
/// `ignore_alphabet` is set, in which case the rows and columns are simply
/// reversed. The total IC and consensus are recomputed; the background is
/// kept as is.
pub fn reverse_complement(motif: &Motif, ignore_alphabet: bool) -> Result<Motif> {
    if !motif.alphabet.is_nucleotide() && !ignore_alphabet {
        return Err(Error::Unsupported(
            "only DNA/RNA motifs can be reverse complemented".into(),
        ));
    }
    let ncol = motif.width();
    let gap = GapInfo {
        isgapped: motif.gap.isgapped,
        gaploc: motif.gap.gaploc.iter().rev().map(|&g| ncol - g).collect(),
        mingap: motif.gap.mingap.iter().rev().copied().collect(),
        maxgap: motif.gap.maxgap.iter().rev().copied().collect(),
    };
    let p = motif.clone().into_parts();
    Motif::builder(p.matrix.reversed(), p.alphabet)
        .name(p.name)
        .altname(p.altname)
        .family(p.family)
        .organism(p.organism)
        .motif_type(Some(p.motif_type))
        .nsites(p.nsites)
        .pseudocount(p.pseudocount)
        .bkg(Some(p.bkg))
        .bkgsites(p.bkgsites)
        .strand(p.strand)
        .pval(p.pval)
        .qval(p.qval)
        .eval(p.eval)
        .extrainfo(p.extrainfo)
        .gap(gap)
        .build()
}

/// Which ends of a motif to trim.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TrimFrom {
    /// Trim both ends.
    #[default]
    Both,
    /// Trim only the start.
    Left,
    /// Trim only the end.
    Right,
}

impl TrimFrom {
    /// Parse `"both"`, `"left"` or `"right"`.
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "both" => Ok(TrimFrom::Both),
            "left" => Ok(TrimFrom::Left),
            "right" => Ok(TrimFrom::Right),
            _ => arg_err(format!(
                "incorrect trim_from '{s}' (must be one of both, left, right)"
            )),
        }
    }

    fn left(self) -> bool {
        matches!(self, TrimFrom::Both | TrimFrom::Left)
    }

    fn right(self) -> bool {
        matches!(self, TrimFrom::Both | TrimFrom::Right)
    }
}

/// Keep only the columns in `keep`, shifting gap positions to match and
/// dropping gaps that fall outside.
fn keep_columns(motif: &Motif, keep: std::ops::Range<usize>) -> (Matrix, GapInfo) {
    let width = keep.len();
    let mut gap = GapInfo {
        isgapped: false,
        ..GapInfo::default()
    };
    for (i, &g) in motif.gap.gaploc.iter().enumerate() {
        if g > keep.start && g - keep.start < width {
            gap.gaploc.push(g - keep.start);
            gap.mingap.push(motif.gap.mingap[i]);
            gap.maxgap.push(motif.gap.maxgap[i]);
        }
    }
    gap.isgapped = motif.gap.isgapped && !gap.gaploc.is_empty();
    (motif.matrix.slice_columns(keep), gap)
}

/// Rebuild the derived fields, the total IC and consensus, after the matrix
/// changed, leaving everything else as it was.
fn refresh(motif: &Motif, matrix: Matrix, gap: GapInfo, ic_type: MotifType) -> Result<Motif> {
    let mut p = motif.clone().into_parts();
    let nsites = Some(p.nsites.unwrap_or(100.0));
    p.icscore = Motif::compute_icscore(&matrix, &p.bkg, ic_type, p.pseudocount, nsites);
    p.consensus = Motif::compute_consensus(&matrix, &p.alphabet, ic_type, p.pseudocount);
    p.matrix = matrix;
    p.gap = gap;
    Motif::from_parts(p)
}

/// Trim low information content columns from the edges of a motif.
///
/// Columns with an IC below `min_ic` are removed from the chosen ends until a
/// column at or above it is reached. Returns `None` if no column is left.
pub fn trim_motif(motif: &Motif, min_ic: f64, trim_from: TrimFrom) -> Result<Option<Motif>> {
    let nsites = motif.nsites.unwrap_or(100.0);
    let scores: Vec<f64> = motif
        .matrix
        .columns()
        .map(|c| {
            position_icscore(
                c,
                &motif.bkg,
                motif.motif_type,
                motif.pseudocount,
                nsites,
                false,
            )
        })
        .collect();
    let n = scores.len();
    let mut left = 0;
    if trim_from.left() {
        while left < n && scores[left] < min_ic {
            left += 1;
        }
    }
    let mut right = n;
    if trim_from.right() {
        while right > left && scores[right - 1] < min_ic {
            right -= 1;
        }
    }
    if right <= left {
        return Ok(None);
    }
    let (matrix, gap) = keep_columns(motif, left..right);
    refresh(motif, matrix, gap, motif.motif_type).map(Some)
}

/// Trim the edges of a contribution weight matrix.
///
/// A column is kept if the sum of its absolute values reaches the cutoff, which
/// is `abs_threshold` if given and otherwise `trim_threshold` times the largest
/// column sum. Columns are only removed from the chosen ends, up to the first
/// kept column from each, so everything in between is retained. Returns `None`
/// if no column is left.
pub fn trim_cwm(
    motif: &Motif,
    trim_threshold: f64,
    abs_threshold: Option<f64>,
    trim_from: TrimFrom,
) -> Result<Option<Motif>> {
    if motif.motif_type != MotifType::Cwm {
        return Err(Error::Unsupported(format!(
            "trim_cwm requires CWM motifs; got type {}",
            motif.motif_type
        )));
    }
    if trim_threshold.is_nan() || trim_threshold < 0.0 {
        return arg_err("trim_threshold must be a single non-negative number");
    }
    if let Some(a) = abs_threshold
        && (a.is_nan() || a < 0.0)
    {
        return arg_err("abs_threshold must be a single non-negative number");
    }

    let col_abs: Vec<f64> = motif
        .matrix
        .columns()
        .map(|c| c.iter().fold(0.0, |a, &v| a + v.abs()))
        .collect();
    let cutoff = match abs_threshold {
        Some(a) => a,
        None => {
            let peak = col_abs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            if peak == 0.0 {
                return Ok(None);
            }
            trim_threshold * peak
        }
    };
    let surviving: Vec<usize> = (0..col_abs.len())
        .filter(|&i| col_abs[i] >= cutoff)
        .collect();
    let (Some(&first), Some(&last)) = (surviving.first(), surviving.last()) else {
        return Ok(None);
    };
    let left = if trim_from.left() { first } else { 0 };
    let right = if trim_from.right() {
        last + 1
    } else {
        col_abs.len()
    };
    let (matrix, gap) = keep_columns(motif, left..right);
    refresh(motif, matrix, gap, MotifType::Cwm).map(Some)
}

/// Set probabilities below `pct_tolerance` to zero and renormalize each column.
/// The motif is converted to a PPM without a pseudocount for the rounding and
/// then back to its original type.
pub fn round_motif(motif: &Motif, pct_tolerance: f64) -> Result<Motif> {
    if !(0.0..=1.0).contains(&pct_tolerance) {
        return arg_err("pct_tolerance must be a numeric value between 0 and 1");
    }
    let ppm = convert_type(motif, MotifType::Ppm, Some(0.0), false)?;
    let mut p = ppm.into_parts();
    p.matrix = p.matrix.map_columns(|c| {
        let kept: Vec<f64> = c
            .iter()
            .map(|&v| if v < pct_tolerance { 0.0 } else { v })
            .collect();
        let s = crate::matrix::sum(&kept);
        kept.iter().map(|&v| v / s).collect()
    });
    let rounded = Motif::from_parts(p)?;
    let mut out = convert_type(&rounded, motif.motif_type, None, false)?.into_parts();
    out.recompute_derived();
    Motif::from_parts(out)
}

/// Re-apply the pseudocount by converting to counts and back. A pseudocount of
/// zero is replaced by 1. PCM motifs are returned unchanged.
pub fn normalize(motif: &Motif) -> Result<Motif> {
    let pseudo = if motif.pseudocount == 0.0 {
        1.0
    } else {
        motif.pseudocount
    };
    let pcm = convert_type(motif, MotifType::Pcm, None, false)?;
    let mut out = convert_type(&pcm, motif.motif_type, Some(pseudo), false)?.into_parts();
    out.recompute_derived();
    Motif::from_parts(out)
}

/// Add gaps to a motif: gap `i` follows column `gaploc[i]`, counting from 1,
/// and spans `mingap[i]` to `maxgap[i]` positions. Without `gaploc`, the gap
/// goes in the middle of the motif. Shorter lists are recycled to the length of
/// the longest.
pub fn add_gap(
    motif: &Motif,
    gaploc: Option<&[usize]>,
    mingap: &[usize],
    maxgap: &[usize],
) -> Result<Motif> {
    let middle = [motif.width() / 2];
    let gaploc = gaploc.unwrap_or(&middle);
    if gaploc.is_empty() || mingap.is_empty() || maxgap.is_empty() {
        return arg_err("gaploc, mingap and maxgap must not be empty");
    }
    let n = gaploc.len().max(mingap.len()).max(maxgap.len());
    let recycle = |v: &[usize]| (0..n).map(|i| v[i % v.len()]).collect::<Vec<usize>>();
    let gap = GapInfo {
        isgapped: true,
        gaploc: recycle(gaploc),
        mingap: recycle(mingap),
        maxgap: recycle(maxgap),
    };
    if gap.mingap.iter().zip(&gap.maxgap).any(|(lo, hi)| lo > hi) {
        return arg_err("mingap cannot be greater than the corresponding maxgap");
    }
    if gap.maxgap.contains(&0) {
        return arg_err("maxgap must be greater than 0");
    }
    let mut p = motif.clone().into_parts();
    p.gap = gap;
    Motif::from_parts(p)
}

/// Turn a motif's gaps off. Unless `delete`, the gap positions and sizes are
/// kept so that [`add_gap`] can restore them.
pub fn ungap(motif: &Motif, delete: bool) -> Motif {
    let mut p = motif.clone().into_parts();
    if delete {
        p.gap = GapInfo::default();
    } else {
        p.gap.isgapped = false;
    }
    Motif::from_parts(p).expect("removing gaps keeps a motif valid")
}

/// Switch a motif between the DNA and RNA alphabets. The matrix and
/// background stay as they are.
pub fn switch_alph(motif: &Motif) -> Result<Motif> {
    let (to, from_letter, to_letter) = match motif.alphabet {
        Alphabet::Dna => (Alphabet::Rna, 'T', 'U'),
        Alphabet::Rna => (Alphabet::Dna, 'U', 'T'),
        _ => {
            return arg_err(format!(
                "motif '{}' has alphabet '{}'; will not convert",
                motif.name, motif.alphabet
            ));
        }
    };
    let mut p = motif.clone().into_parts();
    p.alphabet = to;
    p.consensus = p
        .consensus
        .map(|c| c.replace(from_letter, &to_letter.to_string()));
    Motif::from_parts(p)
}
