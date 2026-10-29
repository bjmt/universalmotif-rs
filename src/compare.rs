//! Motif comparison with [`compare_motifs`].
//!
//! Two motifs are compared by sliding one along the other and summing, over the
//! overlapping columns, a score for each column pair; for DNA and RNA the
//! reverse complement is tried as well. The best alignment is kept. Its
//! significance comes from the distribution of that sum under a null model of
//! target columns, computed exactly by convolution, Bonferroni-corrected for
//! the number of alignments tried.
//!
//! Columns are scored by Pearson correlation by default, centred on 1/n for an
//! n-letter alphabet so that it is the correlation of probability columns, or
//! by Euclidean distance, which is negated so that the best alignment always
//! has the highest sum. With [`Overhang::Background`], the columns of either
//! motif hanging over the ends of the other are compared with the background
//! too.
//!
//! Motifs of any alphabet can be compared. Reverse complements are only tried,
//! and the parametric null is only available, for DNA and RNA.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::convert::{consensus_letter, consensus_letter_aa, convert_type};
use crate::error::{Error, Result, arg_err};
use crate::matrix::which_max;
use crate::metrics::ColumnMetric;
use crate::motif::{Motif, MotifType, Strand};
use crate::pvalue::with_threads;

/// How close to the edge between two bins, as a fraction of a bin, a score
/// must be for [`Scorer::careful`] to decide its bin.
const EDGE: f64 = 1e-9;
const ENUM_GRID_K: usize = 5;
const ENUM_DIRICHLET_N: f64 = 4.0;
const TIE_TOLERANCE: f64 = 1e-10;

/// A motif column of probabilities. Fixed-size arrays of 4 letters for DNA and
/// RNA and of 20 for amino acids are used so that code is compiled for those
/// sizes; other alphabets use vectors.
trait Col: Clone + Send + Sync {
    fn from_slice(v: &[f64]) -> Self;
    fn values(&self) -> &[f64];
}

impl<const N: usize> Col for [f64; N] {
    fn from_slice(v: &[f64]) -> Self {
        let mut a = [0.0; N];
        a.copy_from_slice(v);
        a
    }
    fn values(&self) -> &[f64] {
        self
    }
}

impl Col for Vec<f64> {
    fn from_slice(v: &[f64]) -> Self {
        v.to_vec()
    }
    fn values(&self) -> &[f64] {
        self
    }
}

/// A way of scoring column pairs: columns are prepared once, so that scoring
/// a pair does as little as possible, and scores lie in `LO..=HI`, higher
/// being more alike.
trait Scorer {
    const LO: f64;
    const HI: f64;
    /// Bins that column scores are rounded to, evenly over `LO..=HI`, for
    /// the null distributions.
    const BINS: usize;
    /// Whether alignments of overlaps alone are chosen by their mean score
    /// rather than their sum, as scores that are never positive would
    /// otherwise favour the shortest overlaps.
    const BY_MEAN: bool;
    fn prepare<C: Col>(c: &C) -> C;
    fn score<C: Col>(a: &C, b: &C) -> f64;
    /// Columns prepared for [`Scorer::careful`].
    type Careful<C: Col>: Send + Sync;
    fn careful_prepare<C: Col>(c: &C) -> Self::Careful<C>;
    /// The score of two columns by the plainest formula. Column medians and
    /// the alignments kept are scored this way, and so are column scores on
    /// the edge between two bins, so that rounding errors cannot send a
    /// score that is on the edge in theory, such as a correlation of -0.5,
    /// either way, and results do not depend on the quicker scores.
    fn careful<C: Col>(a: &Self::Careful<C>, b: &Self::Careful<C>) -> f64;
}

/// A probability column centred on the uniform value, with its sum of
/// squares, for correlating it carefully.
struct Centred<C> {
    values: C,
    ss: f64,
}

/// Pearson correlation. Columns are centred on the uniform value, 1/n (0.25
/// for DNA and RNA), and scaled to unit length, which makes a correlation a
/// dot product. Uniform columns become zeros, so that they correlate 0 with
/// anything.
struct PccScore;

impl Scorer for PccScore {
    const LO: f64 = -1.0;
    const HI: f64 = 1.0;
    // Correlations between real columns pile up near 1 and -1. Among all
    // pairs of JASPAR vertebrate motifs, rounding to 51 bins lost 6% of
    // the comparisons significant at q 0.05 and moved 1% of P-values from
    // 1e-6 to 0.05 by more than an order of magnitude; 201 bins lose a
    // third as many and move them by at most about half that, at about
    // three times the cost of the P-values.
    const BINS: usize = 201;
    const BY_MEAN: bool = false;

    fn prepare<C: Col>(c: &C) -> C {
        let v = c.values();
        let centre = 1.0 / v.len() as f64;
        let ss = v.iter().fold(0.0, |a, &x| a + (x - centre) * (x - centre));
        let norm = ss.sqrt();
        let unit: Vec<f64> = v
            .iter()
            .map(|&x| if ss < 1e-12 { 0.0 } else { (x - centre) / norm })
            .collect();
        C::from_slice(&unit)
    }

    #[inline(always)]
    fn score<C: Col>(a: &C, b: &C) -> f64 {
        let (a, b) = (a.values(), b.values());
        let mut s = 0.0;
        for i in 0..a.len() {
            s += a[i] * b[i];
        }
        s
    }

    type Careful<C: Col> = Centred<C>;

    fn careful_prepare<C: Col>(c: &C) -> Centred<C> {
        let v = c.values();
        let centre = 1.0 / v.len() as f64;
        let mut ss = 0.0;
        let centred: Vec<f64> = v
            .iter()
            .map(|&x| {
                let a = x - centre;
                ss += a * a;
                a
            })
            .collect();
        Centred {
            values: C::from_slice(&centred),
            ss,
        }
    }

    fn careful<C: Col>(a: &Centred<C>, b: &Centred<C>) -> f64 {
        if a.ss < 1e-12 || b.ss < 1e-12 {
            return 0.0;
        }
        Self::score(&a.values, &b.values) / (a.ss * b.ss).sqrt()
    }
}

/// Euclidean distance between probability columns, negated so that more
/// alike columns score higher. No two columns are further apart than the
/// square root of 2.
struct EuclScore;

impl Scorer for EuclScore {
    const LO: f64 = -std::f64::consts::SQRT_2;
    const HI: f64 = 0.0;
    // Distances between real columns spread over the whole range. Among
    // all pairs of JASPAR vertebrate motifs with P-values from 1e-6 to
    // 0.05, rounding to 51 bins moved 1% of the P-values by more than 0.13
    // orders of magnitude and 201 bins by more than 0.03, at about three
    // times the cost.
    const BINS: usize = 201;
    const BY_MEAN: bool = true;

    fn prepare<C: Col>(c: &C) -> C {
        c.clone()
    }

    #[inline(always)]
    fn score<C: Col>(a: &C, b: &C) -> f64 {
        let (a, b) = (a.values(), b.values());
        let mut s = 0.0;
        for i in 0..a.len() {
            let d = a[i] - b[i];
            s += d * d;
        }
        -s.sqrt()
    }

    type Careful<C: Col> = C;

    fn careful_prepare<C: Col>(c: &C) -> C {
        c.clone()
    }

    fn careful<C: Col>(a: &C, b: &C) -> f64 {
        Self::score(a, b)
    }
}

/// The bin of a sum of `len` column scores.
fn sum_score_to_bin<S: Scorer>(s: f64, len: usize) -> usize {
    let per_col = ((s / len as f64 - S::LO) / (S::HI - S::LO)).clamp(0.0, 1.0);
    let bin = (per_col * (S::BINS - 1) as f64 * len as f64 + 0.5).floor() as i64;
    bin.clamp(0, (len * (S::BINS - 1)) as i64) as usize
}

/// A column prepared both for quick scoring and for careful scoring.
struct Both<C: Col, S: Scorer> {
    quick: C,
    careful: S::Careful<C>,
}

impl<C: Col, S: Scorer> Both<C, S> {
    fn new(c: &C) -> Self {
        Both {
            quick: S::prepare(c),
            careful: S::careful_prepare(c),
        }
    }

    /// The bin of the score of the two columns. Away from the edges between
    /// bins, quicker arithmetic than [`sum_score_to_bin`] finds the same bin.
    #[inline(always)]
    fn bin(&self, t: &Both<C, S>) -> usize {
        let s = S::score(&self.quick, &t.quick);
        let x = (s - S::LO) * ((S::BINS - 1) as f64 / (S::HI - S::LO)) + 0.5;
        let f = x.floor();
        if x - f < EDGE || f + 1.0 - x < EDGE {
            return sum_score_to_bin::<S>(S::careful(&self.careful, &t.careful), 1);
        }
        (f as i64).clamp(0, S::BINS as i64 - 1) as usize
    }
}

/// Reverse complement of DNA or RNA columns. In A, C, G, T order the complement
/// reverses each column.
fn rc_columns<C: Col>(m: &[C]) -> Vec<C> {
    m.iter()
        .rev()
        .map(|c| {
            let mut v = c.values().to_vec();
            v.reverse();
            C::from_slice(&v)
        })
        .collect()
}

/// The 56 columns of the K = 5 simplex grid with Dirichlet-multinomial
/// weights, for the parametric null.
fn enum_columns(bkg: &[f64; 4]) -> Vec<([f64; 4], f64)> {
    let lgamma = libm::lgamma;
    let alpha: Vec<f64> = bkg
        .iter()
        .map(|b| (ENUM_DIRICHLET_N * b).max(0.1))
        .collect();
    let alpha_sum: f64 = alpha.iter().sum();
    let k = ENUM_GRID_K as f64;
    let mut cols = Vec::with_capacity(56);
    for a in 0..=ENUM_GRID_K {
        for b in 0..=ENUM_GRID_K - a {
            for c in 0..=ENUM_GRID_K - a - b {
                let d = ENUM_GRID_K - a - b - c;
                let counts = [a, b, c, d].map(|x| x as f64);
                let mut lw = lgamma(k + 1.0);
                for i in 0..4 {
                    lw -= lgamma(counts[i] + 1.0);
                    lw += lgamma(counts[i] + alpha[i]) - lgamma(alpha[i]);
                }
                lw -= lgamma(k + alpha_sum) - lgamma(alpha_sum);
                cols.push((counts.map(|x| x / k), lw));
            }
        }
    }
    let max = cols.iter().map(|c| c.1).fold(f64::NEG_INFINITY, f64::max);
    let mut z = 0.0;
    for c in cols.iter_mut() {
        c.1 = (c.1 - max).exp();
        z += c.1;
    }
    for c in cols.iter_mut() {
        c.1 /= z;
    }
    cols
}

/// The null model column correlations are compared against.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NullModel {
    /// Columns of all the motifs being compared.
    #[default]
    Empirical,
    /// A Dirichlet-multinomial grid of columns around the background.
    Parametric,
}

/// What the columns of either motif hanging over the ends of the other are
/// compared with.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Overhang {
    /// Nothing: only the overlap is scored.
    #[default]
    Ignore,
    /// The background, so that every column of both motifs is compared with
    /// something in every alignment. A query column is compared with the
    /// target column facing it or, past the target's ends, with the
    /// background, and the same goes for each target column, which makes as
    /// many column scores as the two motifs have columns between them,
    /// whatever the alignment. Each overlapping pair is thus counted from
    /// both sides. The alignment with the highest total is kept, and the
    /// score reported is its mean, so that partial overlaps are penalized
    /// by however much their overhanging columns differ from the
    /// background. With Pearson correlation and a uniform background, the
    /// overhanging columns score 0, which leaves the choice of alignment as
    /// it is without them and only lowers the reported score of partial
    /// overlaps.
    Background,
}

/// Options for comparisons.
#[derive(Debug, Clone, PartialEq)]
pub struct CompareOptions {
    /// Fewest overlapping columns an alignment may have, lowered to the width
    /// of the narrower motif when needed.
    pub min_overlap: usize,
    /// Also try the reverse complement of each query.
    pub rc: bool,
    /// The null model for P-values; the empirical one by default.
    pub null: NullModel,
    /// Background for the parametric null and for
    /// [`Overhang::Background`], in A, C, G, T order. Other alphabets use a
    /// uniform background.
    pub bkg: [f64; 4],
    /// How column pairs are scored: [`ColumnMetric::Pcc`], their Pearson
    /// correlation, by default, or [`ColumnMetric::Eucl`], their Euclidean
    /// distance. Distances are negated while aligning, so that the best
    /// alignment always has the highest sum of column scores, but reported
    /// as they are, lower meaning more alike. Since summed distances favour
    /// short overlaps, Euclidean alignments are chosen by their mean column
    /// distance unless complete scores or background overhangs are used,
    /// both of which make alignments cover the same number of columns.
    pub method: ColumnMetric,
    /// What overhanging columns are compared with; nothing by default.
    /// Cannot be combined with `complete`.
    pub overhang: Overhang,
    /// Score whole queries rather than overlaps, as TOMTOM does by default.
    /// Each query column's scores are measured against the median of its
    /// scores with the null columns, so that a typical column
    /// scores zero, and query columns hanging over the ends of the target
    /// count as typical columns. The best alignment is then the one with the
    /// highest sum, and its P-value is that of the best of all the alignments
    /// tried, assuming they are independent, rather than a Bonferroni
    /// correction. This favours alignments covering more of the query and
    /// makes comparisons depend on which motif is the query. Off by default,
    /// when only the overlap counts.
    pub complete: bool,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
}

impl Default for CompareOptions {
    fn default() -> Self {
        CompareOptions {
            min_overlap: 5,
            rc: true,
            null: NullModel::Empirical,
            bkg: [0.25; 4],
            method: ColumnMetric::Pcc,
            overhang: Overhang::Ignore,
            complete: false,
            nthreads: 1,
        }
    }
}

/// The best alignment of a target motif to a query motif.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Alignment {
    /// Index of the query motif.
    pub query: usize,
    /// Index of the target motif.
    pub target: usize,
    /// Sum of column scores over the overlap, with Euclidean distances
    /// negated.
    pub score: f64,
    /// The score the alignment was chosen by: `score`, its mean over the
    /// overlap for Euclidean distances, with [`CompareOptions::complete`] the
    /// sum over the overlap of each query column's score less its median, or
    /// with [`Overhang::Background`] the sum over every column of both
    /// motifs.
    pub complete_score: f64,
    /// Number of overlapping columns, or 0 if no alignment was possible.
    pub overlap: usize,
    /// Column of the oriented query facing target column 0; can be negative.
    pub offset: i64,
    /// `Minus` if the query's reverse complement aligned best.
    pub strand: Strand,
    /// First overlapping column of the oriented query.
    pub query_start: usize,
    /// Number of alignments tried.
    pub n_tested: usize,
}

/// The alphabet shared by all the motifs.
fn shared_alphabet(motifs: &[Motif]) -> Result<Alphabet> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs to compare");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if first.alphabet.len() < 2 {
        return arg_err("motifs need at least two letters to be compared");
    }
    Ok(first.alphabet.clone())
}

/// Check the options that the alignments and P-values both depend on.
fn check_options(opts: &CompareOptions, alphabet: &Alphabet) -> Result<()> {
    if !matches!(opts.method, ColumnMetric::Pcc | ColumnMetric::Eucl) {
        return Err(Error::Unsupported(format!(
            "motifs can only be compared by PCC or EUCL, not {}",
            opts.method
        )));
    }
    if opts.overhang == Overhang::Background {
        if opts.complete {
            return arg_err("complete scores and background overhangs cannot be combined");
        }
        let b = &opts.bkg;
        if alphabet.is_nucleotide()
            && (b.iter().any(|v| !(v.is_finite() && *v >= 0.0))
                || (b.iter().sum::<f64>() - 1.0).abs() > 1e-6)
        {
            return arg_err("bkg must be four non-negative values summing to 1");
        }
    }
    Ok(())
}

/// Run `$f::<C>` with the column type that suits the alphabet.
macro_rules! by_alphabet {
    ($alphabet:expr, $f:ident($($arg:expr),*)) => {
        match $alphabet.len() {
            4 => $f::<[f64; 4]>($($arg),*),
            20 => $f::<[f64; 20]>($($arg),*),
            _ => $f::<Vec<f64>>($($arg),*),
        }
    };
}

/// Run `$f::<C, S>` with the column type that suits the alphabet and the
/// scorer of `opts.method`, which [`check_options`] has made PCC or EUCL.
macro_rules! by_method {
    ($alphabet:expr, $opts:expr, $f:ident($($arg:expr),*)) => {
        match ($alphabet.len(), $opts.method) {
            (4, ColumnMetric::Eucl) => $f::<[f64; 4], EuclScore>($($arg),*),
            (20, ColumnMetric::Eucl) => $f::<[f64; 20], EuclScore>($($arg),*),
            (_, ColumnMetric::Eucl) => $f::<Vec<f64>, EuclScore>($($arg),*),
            (4, _) => $f::<[f64; 4], PccScore>($($arg),*),
            (20, _) => $f::<[f64; 20], PccScore>($($arg),*),
            _ => $f::<Vec<f64>, PccScore>($($arg),*),
        }
    };
}

/// PPM columns of the motifs.
fn columns<C: Col>(motifs: &[Motif]) -> Result<Vec<Vec<C>>> {
    motifs
        .iter()
        .map(|m| {
            let ppm = convert_type(m, MotifType::Ppm, None, false)?;
            Ok(ppm.matrix.columns().map(C::from_slice).collect())
        })
        .collect()
}

/// The background column that overhanging columns are compared with.
fn bkg_column<C: Col>(opts: &CompareOptions, n: usize) -> C {
    if n == 4 {
        C::from_slice(&opts.bkg)
    } else {
        C::from_slice(&vec![1.0 / n as f64; n])
    }
}

/// The median score of each column of each motif with the null columns,
/// which complete scores subtract so that a typical column scores zero. The
/// null is the columns of all the motifs, in both orientations when reverse
/// complements are tried, or the parametric grid with its weights.
fn column_medians<C: Col, S: Scorer>(
    cols: &[Vec<C>],
    opts: &CompareOptions,
    rc_ok: bool,
) -> Vec<Vec<f64>> {
    let null: Vec<(S::Careful<C>, f64)> = match opts.null {
        NullModel::Parametric => enum_columns(&opts.bkg)
            .into_iter()
            .map(|(c, w)| (S::careful_prepare(&C::from_slice(&c)), w))
            .collect(),
        NullModel::Empirical => {
            let mut v: Vec<(S::Careful<C>, f64)> = cols
                .iter()
                .flatten()
                .map(|c| (S::careful_prepare(c), 1.0))
                .collect();
            if rc_ok {
                v.extend(
                    cols.iter()
                        .flat_map(|m| rc_columns(m))
                        .map(|c| (S::careful_prepare(&c), 1.0)),
                );
            }
            v
        }
    };
    cols.par_iter()
        .map(|m| {
            m.iter()
                .map(|q| {
                    let q = S::careful_prepare(q);
                    let mut v: Vec<(f64, f64)> =
                        null.iter().map(|(t, w)| (S::careful(&q, t), *w)).collect();
                    v.sort_by(|a, b| a.0.total_cmp(&b.0));
                    let half = v.iter().map(|x| x.1).sum::<f64>() / 2.0;
                    let mut cum = 0.0;
                    v.iter()
                        .find(|x| {
                            cum += x.1;
                            cum >= half
                        })
                        .map_or(0.0, |x| x.0)
                })
                .collect()
        })
        .collect()
}

/// Running sums of some per-column values, so that the sum over any run of
/// columns is a difference of two entries.
fn running_sums(v: impl Iterator<Item = f64>) -> Vec<f64> {
    let mut out = vec![0.0];
    for x in v {
        out.push(out.last().unwrap() + x);
    }
    out
}

/// The sum of the entries of `sums`, running sums of some column values,
/// outside the `len` columns from `start`.
fn sum_outside(sums: &[f64], start: usize, len: usize) -> f64 {
    sums[start] + (sums[sums.len() - 1] - sums[start + len])
}

/// How alignments are chosen beyond the sum of column scores over the
/// overlap, with each motif's running sums in its own orientation and as its
/// reverse complement has them.
enum Extras {
    /// Nothing more: the overlap's sum, or its mean for Euclidean distances.
    None,
    /// Running sums of column medians, for complete scores.
    Medians(Vec<(Vec<f64>, Vec<f64>)>),
    /// Running sums of column scores with the background, for background
    /// overhangs.
    Background(Vec<(Vec<f64>, Vec<f64>)>),
}

impl Extras {
    fn new<C: Col, S: Scorer>(
        cols: &[Vec<C>],
        n_letters: usize,
        opts: &CompareOptions,
        rc_ok: bool,
    ) -> Extras {
        let both = |m: Vec<f64>, minus: Vec<f64>| {
            (running_sums(m.into_iter()), running_sums(minus.into_iter()))
        };
        if opts.complete {
            Extras::Medians(
                column_medians::<C, S>(cols, opts, rc_ok)
                    .into_iter()
                    .map(|m| {
                        let rev = m.iter().rev().copied().collect();
                        both(m, rev)
                    })
                    .collect(),
            )
        } else if opts.overhang == Overhang::Background {
            let b = S::prepare(&bkg_column::<C>(opts, n_letters));
            let scores =
                |m: &[C]| -> Vec<f64> { m.iter().map(|c| S::score(&S::prepare(c), &b)).collect() };
            Extras::Background(
                cols.iter()
                    .map(|m| {
                        let minus = if rc_ok {
                            scores(&rc_columns(m))
                        } else {
                            Vec::new()
                        };
                        both(scores(m), minus)
                    })
                    .collect(),
            )
        } else {
            Extras::None
        }
    }

    fn rule<S: Scorer>(&self, q: usize, t: usize) -> Rule<'_> {
        match self {
            Extras::None if S::BY_MEAN => Rule::Mean,
            Extras::None => Rule::Sum,
            Extras::Medians(m) => Rule::Complete(&m[q].0, &m[q].1),
            Extras::Background(b) => Rule::Background {
                query: (&b[q].0, &b[q].1),
                target: &b[t].0,
            },
        }
    }
}

/// What an alignment is chosen by: the sum of column scores over the
/// overlap, its mean, the sum less the query's column medians over the
/// overlap, or the sum over every column of both motifs, each overlapping
/// pair counted from both sides and overhanging columns compared with the
/// background, given the running sums of these background scores.
#[derive(Clone, Copy)]
enum Rule<'a> {
    Sum,
    Mean,
    Complete(&'a [f64], &'a [f64]),
    Background {
        query: (&'a [f64], &'a [f64]),
        target: &'a [f64],
    },
}

impl Rule<'_> {
    /// The value of an alignment with sum `overlap_sum` over the `len`
    /// columns from `q_start` of the query in orientation `strand` and from
    /// `t_start` of the target.
    #[inline(always)]
    fn value(
        &self,
        overlap_sum: f64,
        strand: Strand,
        q_start: usize,
        t_start: usize,
        len: usize,
    ) -> f64 {
        let (qs, ts, l) = (q_start, t_start, len);
        match *self {
            Rule::Sum => overlap_sum,
            Rule::Mean => overlap_sum / len as f64,
            Rule::Complete(plus, minus) => {
                let b = if strand == Strand::Plus { plus } else { minus };
                overlap_sum - (b[qs + l] - b[qs])
            }
            Rule::Background {
                query: (plus, minus),
                target,
            } => {
                let b = if strand == Strand::Plus { plus } else { minus };
                2.0 * overlap_sum + sum_outside(b, qs, l) + sum_outside(target, ts, l)
            }
        }
    }
}

/// The best alignment of a target to a query.
struct Best {
    /// The score it was chosen by.
    chosen: f64,
    overlap: usize,
    offset: i64,
    strand: Strand,
    /// First overlapping column of the oriented query.
    query_start: usize,
    /// Number of alignments tried.
    n_tested: usize,
    /// Sum of column scores over the overlap.
    sum: f64,
}

/// The best alignment of `t` to the query, chosen by `rule`.
fn best_alignment<C: Col, S: Scorer>(
    q: &[C],
    q_rc: Option<&[C]>,
    t: &[C],
    min_overlap: usize,
    rule: Rule<'_>,
) -> Best {
    let (wq, wt) = (q.len() as i64, t.len() as i64);
    let mo = (min_overlap as i64).min(wq).min(wt);
    let mut best = Best {
        chosen: f64::NEG_INFINITY,
        overlap: 0,
        offset: 0,
        strand: Strand::Plus,
        query_start: 0,
        n_tested: 0,
        sum: f64::NEG_INFINITY,
    };
    let orientations: Vec<(&[C], Strand)> = match q_rc {
        Some(r) => vec![(q, Strand::Plus), (r, Strand::Minus)],
        None => vec![(q, Strand::Plus)],
    };
    for (qq, strand) in orientations {
        for d in mo - wt..=wq - mo {
            let t_start = if d >= 0 { 0 } else { -d };
            let t_end = (wq - d).min(wt);
            let len = t_end - t_start;
            if len < mo || len <= 0 {
                continue;
            }
            let q_start = d + t_start;
            let overlap_sum = (0..len).fold(0.0, |acc, j| {
                acc + S::score(&qq[(q_start + j) as usize], &t[(t_start + j) as usize])
            });
            let (qs, ts, l) = (q_start as usize, t_start as usize, len as usize);
            let s = rule.value(overlap_sum, strand, qs, ts, l);
            best.n_tested += 1;
            // Ties (such as both strands of a palindrome) go to the longer
            // overlap, so that comparing the two motifs either way round gives
            // the same mean score, and then to the first alignment, rather than
            // to whichever rounds higher.
            let tied_longer = (s - best.chosen).abs() <= TIE_TOLERANCE && l > best.overlap;
            if s > best.chosen + TIE_TOLERANCE || tied_longer || best.n_tested == 1 {
                best = Best {
                    chosen: s,
                    overlap: l,
                    offset: d,
                    strand,
                    query_start: qs,
                    n_tested: best.n_tested,
                    sum: overlap_sum,
                };
            }
        }
    }
    best
}

/// Best alignments for `(query, target)` pairs of motif indices. Reverse
/// complements are only tried for DNA and RNA motifs.
#[doc(hidden)]
pub fn align(
    motifs: &[Motif],
    pairs: &[(usize, usize)],
    opts: &CompareOptions,
) -> Result<Vec<Alignment>> {
    let alphabet = shared_alphabet(motifs)?;
    check_options(opts, &alphabet)?;
    by_method!(
        alphabet,
        opts,
        align_with(motifs, pairs, opts, alphabet.is_nucleotide())
    )
}

fn align_with<C: Col, S: Scorer>(
    motifs: &[Motif],
    pairs: &[(usize, usize)],
    opts: &CompareOptions,
    nucleotide: bool,
) -> Result<Vec<Alignment>> {
    if opts.min_overlap < 1 {
        return arg_err("min_overlap must be at least 1");
    }
    let cols: Vec<Vec<C>> = columns(motifs)?;
    let rc_ok = opts.rc && nucleotide;
    let prepare = |m: &[C]| -> (Vec<C>, Vec<S::Careful<C>>) {
        (
            m.iter().map(S::prepare).collect(),
            m.iter().map(S::careful_prepare).collect(),
        )
    };
    let fwd: Vec<_> = cols.iter().map(|m| prepare(m)).collect();
    let rc: Vec<_> = if rc_ok {
        cols.iter().map(|m| prepare(&rc_columns(m))).collect()
    } else {
        Vec::new()
    };
    if pairs
        .iter()
        .any(|&(q, t)| q >= cols.len() || t >= cols.len())
    {
        return arg_err("motif index out of range");
    }
    with_threads(opts.nthreads, || {
        let extras = Extras::new::<C, S>(&cols, motifs[0].alphabet.len(), opts, rc_ok);
        pairs
            .par_iter()
            .map(|&(q, t)| {
                let rule = extras.rule::<S>(q, t);
                let b = best_alignment::<C, S>(
                    &fwd[q].0,
                    rc_ok.then(|| rc[q].0.as_slice()),
                    &fwd[t].0,
                    opts.min_overlap,
                    rule,
                );
                // The alignment kept is scored again by careful scores, so
                // that the scores reported and binned for P-values do not
                // depend on how quickly they were found.
                let qq = if b.strand == Strand::Minus {
                    &rc[q].1
                } else {
                    &fwd[q].1
                };
                let ts = (-b.offset).max(0) as usize;
                let sum = (0..b.overlap).fold(0.0, |acc, j| {
                    acc + S::careful(&qq[b.query_start + j], &fwd[t].1[ts + j])
                });
                let chosen = if b.overlap > 0 {
                    rule.value(sum, b.strand, b.query_start, ts, b.overlap)
                } else {
                    b.chosen
                };
                Alignment {
                    query: q,
                    target: t,
                    score: if b.overlap > 0 { sum } else { b.sum },
                    complete_score: chosen,
                    overlap: b.overlap,
                    offset: b.offset,
                    strand: b.strand,
                    query_start: b.query_start,
                    n_tested: b.n_tested,
                }
            })
            .collect()
    })
}

fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    // Zero products add nothing, so only the span of nonzero bins of `b` is
    // visited.
    let Some(b_lo) = b.iter().position(|&v| v != 0.0) else {
        return out;
    };
    let b_hi = b.iter().rposition(|&v| v != 0.0).unwrap();
    let b = &b[b_lo..=b_hi];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0.0 {
            continue;
        }
        let dst = &mut out[i + b_lo..i + b_lo + b.len()];
        for (d, &bj) in dst.iter_mut().zip(b) {
            *d += ai * bj;
        }
    }
    out
}

/// The probability of a score in `bin` or above. A score beyond what the null
/// can reach gets the probability of the null's highest score rather than 0,
/// since a null built from a coarse grid of columns may not reach scores that
/// real columns do.
fn tail_prob(pmf: &[f64], bin: usize) -> f64 {
    let Some(top) = pmf.iter().rposition(|&p| p > 0.0) else {
        return 0.0;
    };
    pmf[bin.min(top)..]
        .iter()
        .fold(0.0, |a, &b| a + b)
        .clamp(0.0, 1.0)
}

/// Where a motif sits in a shared column frame, as laid out by
/// [`anchor_layout`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct Placement {
    /// Empty columns to the left of the motif.
    pub left: usize,
    /// Whether the motif is shown as its reverse complement.
    pub rc: bool,
}

/// Motifs laid out in a common column frame by [`anchor_layout`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Layout {
    /// Where each motif sits, in input order.
    pub placements: Vec<Placement>,
    /// The width of the frame, in columns.
    pub width: usize,
}

/// Lay motifs out in a common column frame by aligning each one to
/// `motifs[anchor]`. Motifs that do not overlap the anchor by at least
/// `opts.min_overlap` columns, lowered to the width of the narrower of the
/// two as in comparisons, are placed at the anchor's left edge.
pub fn anchor_layout(motifs: &[Motif], anchor: usize, opts: &CompareOptions) -> Result<Layout> {
    if anchor >= motifs.len() {
        return arg_err("anchor index out of range");
    }
    let pairs: Vec<(usize, usize)> = (0..motifs.len()).map(|t| (anchor, t)).collect();
    let al = align(motifs, &pairs, opts)?;
    let widths: Vec<i64> = motifs.iter().map(|m| m.width() as i64).collect();
    let wa = widths[anchor];
    let mut frame = vec![0i64; motifs.len()];
    let mut rc = vec![false; motifs.len()];
    for (i, a) in al.iter().enumerate() {
        let needed = opts
            .min_overlap
            .min(widths[anchor] as usize)
            .min(widths[i] as usize);
        if i == anchor || a.overlap < needed {
            continue;
        }
        let l = a.overlap as i64;
        let q_start = a.query_start as i64;
        let t_start = (-a.offset).max(0);
        if a.strand == Strand::Minus {
            frame[i] = (wa - q_start - l) - (widths[i] - t_start - l);
            rc[i] = true;
        } else {
            frame[i] = q_start - t_start;
        }
    }
    let lo = *frame.iter().min().unwrap();
    let hi = frame.iter().zip(&widths).map(|(f, w)| f + w).max().unwrap();
    let placements = frame
        .iter()
        .zip(rc)
        .map(|(&f, rc)| Placement {
            left: (f - lo) as usize,
            rc,
        })
        .collect();
    Ok(Layout {
        placements,
        width: (hi - lo) as usize,
    })
}

/// Bonferroni-corrected P-values of alignments. The parametric null is only
/// available for DNA and RNA motifs.
pub(crate) fn alignment_pvalues(
    motifs: &[Motif],
    alignments: &[Alignment],
    opts: &CompareOptions,
) -> Result<Vec<f64>> {
    let alphabet = shared_alphabet(motifs)?;
    check_options(opts, &alphabet)?;
    if opts.null == NullModel::Parametric && !alphabet.is_nucleotide() {
        return Err(Error::Unsupported(format!(
            "the parametric null is only available for DNA and RNA motifs; use the empirical null for {alphabet} motifs"
        )));
    }
    by_method!(
        alphabet,
        opts,
        alignment_pvalues_with(motifs, alignments, opts)
    )
}

fn alignment_pvalues_with<C: Col, S: Scorer>(
    motifs: &[Motif],
    alignments: &[Alignment],
    opts: &CompareOptions,
) -> Result<Vec<f64>> {
    let cols: Vec<Vec<C>> = columns(motifs)?;
    let flat: Vec<Both<C, S>> = cols.iter().flatten().map(Both::new).collect();
    let grid: Option<Vec<(Both<C, S>, f64)>> = (opts.null == NullModel::Parametric).then(|| {
        enum_columns(&opts.bkg)
            .into_iter()
            .map(|(c, w)| (Both::new(&C::from_slice(&c)), w))
            .collect()
    });
    let col_pmf = |q: &C| -> Vec<f64> {
        let mut pmf = vec![0.0; S::BINS];
        let q = Both::<C, S>::new(q);
        match &grid {
            None => {
                if flat.is_empty() {
                    return pmf;
                }
                for t in &flat {
                    pmf[q.bin(t)] += 1.0;
                }
                let inv = 1.0 / flat.len() as f64;
                pmf.iter_mut().for_each(|p| *p *= inv);
            }
            Some(g) => {
                for (c, w) in g {
                    pmf[q.bin(c)] += w;
                }
            }
        }
        pmf
    };
    // Null distributions of an oriented query by start and length, up to
    // the longest length needed from each start (every length if `need` is
    // None).
    let cache = |q: &[C], need: Option<&[usize]>| -> Vec<Vec<Vec<f64>>> {
        let pmfs: Vec<Vec<f64>> = q.iter().map(col_pmf).collect();
        (0..q.len())
            .map(|start| {
                let longest = need.map_or(q.len() - start, |n| n[start]);
                let mut row: Vec<Vec<f64>> = vec![Vec::new(), pmfs[start].clone()];
                for len in 2..=longest {
                    let next = convolve(&row[len - 1], &pmfs[start + len - 1]);
                    row.push(next);
                }
                row
            })
            .collect()
    };

    let mut by_query: Vec<Vec<usize>> = vec![Vec::new(); cols.len()];
    for (k, a) in alignments.iter().enumerate() {
        if a.query >= cols.len() || a.target >= cols.len() {
            return arg_err("motif index out of range");
        }
        by_query[a.query].push(k);
    }
    let rc_ok = opts.rc && motifs.first().is_some_and(|m| m.alphabet.is_nucleotide());
    let n_letters = motifs[0].alphabet.len();
    // With background overhangs, the null distributions of the sums of
    // target columns' scores with the background, by number of columns.
    let overhang_pmfs: Vec<Vec<f64>> = if opts.overhang == Overhang::Background {
        let widest = cols.iter().map(Vec::len).max().unwrap_or(0);
        let one = col_pmf(&bkg_column::<C>(opts, n_letters));
        let mut v = vec![vec![1.0]];
        for k in 1..=widest {
            let next = convolve(&v[k - 1], &one);
            v.push(next);
        }
        v
    } else {
        Vec::new()
    };
    let per_query: Vec<Vec<(usize, f64)>> = with_threads(opts.nthreads, || {
        let extras = Extras::new::<C, S>(&cols, n_letters, opts, rc_ok);
        by_query
            .par_iter()
            .enumerate()
            .map(|(q, ks)| {
                if ks.is_empty() {
                    return Vec::new();
                }
                // Outside complete scoring, only the null of each
                // alignment's own start and overlap is looked up.
                let width = cols[q].len();
                let (mut need_fwd, mut need_rev) = (vec![0; width], vec![0; width]);
                for &k in ks {
                    let a = &alignments[k];
                    if a.overlap > 0 && a.query_start < width && a.n_tested > 0 {
                        let need = if a.strand == Strand::Minus {
                            &mut need_rev
                        } else {
                            &mut need_fwd
                        };
                        let longest = a.overlap.min(width - a.query_start);
                        need[a.query_start] = need[a.query_start].max(longest);
                    }
                }
                let (need_fwd, need_rev) = if opts.complete {
                    (None, None)
                } else {
                    (Some(need_fwd.as_slice()), Some(need_rev.as_slice()))
                };
                let fwd = cache(&cols[q], need_fwd);
                let need_rc = (opts.complete && rc_ok)
                    || ks.iter().any(|&k| alignments[k].strand == Strand::Minus);
                let rev = need_rc.then(|| cache(&rc_columns(&cols[q]), need_rev));
                ks.iter()
                    .map(|&k| {
                        let a = &alignments[k];
                        let mut p = 1.0;
                        if a.overlap > 0 && a.query_start < cols[q].len() && a.n_tested > 0 {
                            let own = if a.strand == Strand::Plus {
                                Some(&fwd)
                            } else {
                                rev.as_ref()
                            };
                            let pmf = own
                                .and_then(|c| c.get(a.query_start))
                                .and_then(|r| r.get(a.overlap));
                            match &extras {
                                Extras::Medians(sums) => {
                                    let (plus, minus) = &sums[q];
                                    let mut orientations = vec![(&fwd, plus.as_slice())];
                                    if let Some(r) = rev.as_ref().filter(|_| rc_ok) {
                                        orientations.push((r, minus.as_slice()));
                                    }
                                    // An observed score beyond what the null
                                    // reaches keeps the probability of the
                                    // null's highest.
                                    let floor = pmf.map_or(0.0, |pmf| {
                                        tail_prob(pmf, sum_score_to_bin::<S>(a.score, a.overlap))
                                    });
                                    p = best_of_all::<S>(
                                        a.complete_score,
                                        cols[q].len(),
                                        cols[a.target].len(),
                                        opts.min_overlap,
                                        &orientations,
                                    )
                                    .max(floor);
                                }
                                Extras::Background(sums) => {
                                    if let Some(pmf) = pmf {
                                        // The query's overhanging columns
                                        // score the same whatever the target,
                                        // so only the target's columns, facing
                                        // the query or the background, vary.
                                        let n_out = cols[a.target].len() - a.overlap;
                                        let t_start = (-a.offset).max(0) as usize;
                                        let t_out =
                                            sum_outside(&sums[a.target].0, t_start, a.overlap);
                                        let tail = overhang_tail::<S>(
                                            pmf,
                                            &overhang_pmfs[n_out],
                                            2.0 * a.score + t_out,
                                            2 * a.overlap + n_out,
                                        );
                                        p = (tail * a.n_tested as f64).clamp(0.0, 1.0);
                                    }
                                }
                                Extras::None => {
                                    if let Some(pmf) = pmf {
                                        let bin = sum_score_to_bin::<S>(a.score, a.overlap);
                                        p = (tail_prob(pmf, bin) * a.n_tested as f64)
                                            .clamp(0.0, 1.0);
                                    }
                                }
                            }
                        }
                        (k, p)
                    })
                    .collect()
            })
            .collect()
    })?;
    let mut out = vec![1.0; alignments.len()];
    for (k, p) in per_query.into_iter().flatten() {
        out[k] = p;
    }
    Ok(out)
}

/// The probability that twice a sum drawn from `overlap`, the null of the
/// overlap's column scores, plus a sum drawn from `outside`, the null of
/// the overhanging target columns' scores with the background, reaches
/// `score`, a sum of `units` column scores. As in [`tail_prob`], a score
/// beyond what the null can reach gets the probability of the null's
/// highest.
fn overhang_tail<S: Scorer>(overlap: &[f64], outside: &[f64], score: f64, units: usize) -> f64 {
    let (Some(top_x), Some(top_y)) = (
        overlap.iter().rposition(|&p| p > 0.0),
        outside.iter().rposition(|&p| p > 0.0),
    ) else {
        return 0.0;
    };
    let bin = sum_score_to_bin::<S>(score, units).min(2 * top_x + top_y);
    let mut at_least = vec![0.0; overlap.len() + 1];
    for i in (0..overlap.len()).rev() {
        at_least[i] = at_least[i + 1] + overlap[i];
    }
    let mut p = 0.0;
    for (j, &py) in outside.iter().enumerate() {
        if py == 0.0 {
            continue;
        }
        // Twice the overlap's bin must reach what the outside leaves.
        let i = bin.saturating_sub(j).div_ceil(2);
        if let Some(&x) = at_least.get(i) {
            p += py * x;
        }
    }
    p.clamp(0.0, 1.0)
}

/// Null distributions of an oriented query's column sums, by first column and
/// number of columns.
type NullCache = Vec<Vec<Vec<f64>>>;

/// The probability that the best of all the alignments of a target of width
/// `wt` to a query of width `wq` reaches a complete score of `score`,
/// treating the alignments as independent. `orientations` holds, for each
/// orientation of the query tried, its null distributions by first column
/// and overlap and the running sums of its column medians.
fn best_of_all<S: Scorer>(
    score: f64,
    wq: usize,
    wt: usize,
    min_overlap: usize,
    orientations: &[(&NullCache, &[f64])],
) -> f64 {
    let (wq, wt) = (wq as i64, wt as i64);
    let mo = (min_overlap as i64).min(wq).min(wt);
    // The log of the probability that no alignment reaches the score.
    let mut log_none = 0.0;
    for (cache, medians) in orientations {
        for d in mo - wt..=wq - mo {
            let t_start = if d >= 0 { 0 } else { -d };
            let len = (wq - d).min(wt) - t_start;
            if len < mo || len <= 0 {
                continue;
            }
            let (start, len) = ((d + t_start) as usize, len as usize);
            // The complete score is the overlap's sum less its medians, so
            // reaching it takes an overlap sum of the score plus the medians.
            let raw = score + medians[start + len] - medians[start];
            // A sum of column scores cannot exceed the best each column can
            // score.
            if raw > len as f64 * S::HI + 1e-9 {
                continue;
            }
            let pmf = &cache[start][len];
            let bin = sum_score_to_bin::<S>(raw, len);
            let tail = pmf
                .get(bin..)
                .map_or(0.0, |t| t.iter().sum::<f64>().min(1.0));
            if tail >= 1.0 {
                return 1.0;
            }
            log_none += (-tail).ln_1p();
        }
    }
    (-log_none.exp_m1()).clamp(0.0, 1.0)
}

/// One row of a comparison: a target motif's best alignment to a query.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Comparison {
    /// Index of the query motif.
    pub query: usize,
    /// Index of the target motif.
    pub target: usize,
    /// Column of the oriented query facing the target's first column, which can
    /// be negative.
    pub offset: i64,
    /// `Minus` if the query's reverse complement aligned best.
    pub strand: Strand,
    /// Number of overlapping columns.
    pub overlap: usize,
    /// Mean column correlation over the overlap, or for
    /// [`ColumnMetric::Eucl`] the mean column distance, lower meaning more
    /// alike. With [`Overhang::Background`], the mean is over every column
    /// of both motifs instead, as [`Overhang::Background`] describes.
    pub score: f64,
    /// The P-value of the alignment's score.
    pub pvalue: f64,
    /// Benjamini-Hochberg adjusted P-value, among the comparisons of the
    /// same query.
    pub qvalue: f64,
    /// Consensus of the aligned query columns, in the query's own orientation.
    pub query_consensus: String,
    /// Consensus of the aligned target columns.
    pub target_consensus: String,
}

/// Consensus letters of some columns, by the same rules as a motif's own
/// consensus for DNA, RNA and proteins, and the most probable letter for other
/// alphabets.
/// The consensus letter of each column of each motif.
fn column_letters<C: Col>(cols: &[Vec<C>], alphabet: &Alphabet) -> Vec<Vec<char>> {
    let letters = alphabet.letters();
    cols.iter()
        .map(|m| {
            m.iter()
                .map(|c| {
                    let v = c.values();
                    match alphabet {
                        Alphabet::Dna => consensus_letter(v, false, MotifType::Ppm, 0.0),
                        Alphabet::Rna => consensus_letter(v, true, MotifType::Ppm, 0.0),
                        Alphabet::Aa => consensus_letter_aa(v, MotifType::Ppm, 0.0),
                        _ => letters[which_max(v)],
                    }
                })
                .collect()
        })
        .collect()
}

/// The consensus of `len` columns from `start`, clipped to the motif.
fn consensus_range(letters: &[char], start: i64, len: usize) -> String {
    let start = (start.max(0) as usize).min(letters.len());
    let len = len.min(letters.len() - start);
    letters[start..start + len].iter().collect()
}

/// Compare each query motif with every motif, itself included.
///
/// Returns one [`Comparison`] per query and target pair, ordered by query and
/// then target, without any filtering.
///
/// # Example
///
/// ```
/// use universalmotif::compare::{CompareOptions, compare_motifs};
/// use universalmotif::create::{CreateOptions, from_strings};
///
/// let opts = CreateOptions { nsites: Some(20.0), pseudocount: 1.0, ..Default::default() };
/// let motifs = [
///     from_strings(&["TGACTCA", "TGAGTCA"], None, &opts)?,
///     from_strings(&["TGACGTCA"], None, &opts)?,
///     from_strings(&["CACGTG"], None, &opts)?,
/// ];
/// // The first motif against all three, best matches first.
/// for c in compare_motifs(&motifs, &[0], &CompareOptions::default())? {
///     println!("{} vs {}: {:.3} (P = {:.2e})", c.query, c.target, c.score, c.pvalue);
/// }
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn compare_motifs(
    motifs: &[Motif],
    queries: &[usize],
    opts: &CompareOptions,
) -> Result<Vec<Comparison>> {
    let n = motifs.len();
    let pairs: Vec<(usize, usize)> = queries
        .iter()
        .flat_map(|&q| (0..n).map(move |t| (q, t)))
        .collect();
    let al = align(motifs, &pairs, opts)?;
    let pvals = alignment_pvalues(motifs, &al, opts)?;
    let alphabet = shared_alphabet(motifs)?;
    by_alphabet!(
        alphabet,
        finish_comparisons(motifs, queries, &al, &pvals, &alphabet, opts)
    )
}

/// The score reported for an alignment: the mean column score, over the
/// overlap or with background overhangs over every column of both motifs,
/// with Euclidean distances as they are, or NaN if no alignment was made.
fn mean_score(a: &Alignment, motifs: &[Motif], opts: &CompareOptions) -> f64 {
    if a.overlap == 0 {
        return f64::NAN;
    }
    let mean = match opts.overhang {
        Overhang::Ignore => a.score / a.overlap as f64,
        Overhang::Background => {
            a.complete_score / (motifs[a.query].width() + motifs[a.target].width()) as f64
        }
    };
    // Subtracting from 0 rather than negating keeps a distance of 0 from
    // being -0.
    if opts.method == ColumnMetric::Eucl {
        0.0 - mean
    } else {
        mean
    }
}

fn finish_comparisons<C: Col>(
    motifs: &[Motif],
    queries: &[usize],
    al: &[Alignment],
    pvals: &[f64],
    alphabet: &Alphabet,
    opts: &CompareOptions,
) -> Result<Vec<Comparison>> {
    let cols: Vec<Vec<C>> = columns(motifs)?;
    let letters = column_letters(&cols, alphabet);

    // The alignments of each query, gathered in one pass.
    let mut by_query: Vec<Vec<usize>> = vec![Vec::new(); motifs.len()];
    for (k, a) in al.iter().enumerate() {
        by_query[a.query].push(k);
    }
    let mut qvals = vec![f64::NAN; al.len()];
    for &q in queries {
        let rows = std::mem::take(&mut by_query[q]);
        let adj = crate::stats::p_adjust_bh(&rows.iter().map(|&k| pvals[k]).collect::<Vec<_>>());
        for (k, v) in rows.into_iter().zip(adj) {
            qvals[k] = v;
        }
    }

    with_threads(opts.nthreads, || {
        al.par_iter()
            .enumerate()
            .map(|(k, a)| {
                let wq = cols[a.query].len() as i64;
                let len = a.overlap as i64;
                let q_start = if a.strand == Strand::Minus {
                    wq - 1 - (a.query_start as i64 + len - 1)
                } else {
                    a.query_start as i64
                };
                let t_start = (-a.offset).max(0);
                Comparison {
                    query: a.query,
                    target: a.target,
                    offset: a.offset,
                    strand: a.strand,
                    overlap: a.overlap,
                    score: mean_score(a, motifs, opts),
                    pvalue: pvals[k],
                    qvalue: qvals[k],
                    query_consensus: consensus_range(&letters[a.query], q_start, a.overlap),
                    target_consensus: consensus_range(&letters[a.target], t_start, a.overlap),
                }
            })
            .collect()
    })
}

/// All-against-all comparison results, indexed `[query][target]`.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ComparisonMatrices {
    /// Mean column scores, as [`Comparison::score`] has them.
    pub scores: Vec<Vec<f64>>,
    /// P-values of the best alignments.
    pub pvalues: Vec<Vec<f64>>,
    /// Benjamini-Hochberg q-values, adjusted within each query.
    pub qvalues: Vec<Vec<f64>>,
}

/// All-against-all comparison. With `need_pvalues`, P-values and q-values
/// are filled in; otherwise they are NaN.
pub fn compare_matrix(
    motifs: &[Motif],
    opts: &CompareOptions,
    need_pvalues: bool,
) -> Result<ComparisonMatrices> {
    let n = motifs.len();
    if n < 2 {
        return arg_err("need at least 2 motifs for an all-against-all comparison");
    }
    let queries: Vec<usize> = (0..n).collect();
    let (scores, pvals, qvals) = if need_pvalues {
        let rows = compare_motifs(motifs, &queries, opts)?;
        let get = |f: &dyn Fn(&Comparison) -> f64| -> Vec<Vec<f64>> {
            rows.chunks(n).map(|r| r.iter().map(f).collect()).collect()
        };
        (get(&|c| c.score), get(&|c| c.pvalue), get(&|c| c.qvalue))
    } else {
        let pairs: Vec<(usize, usize)> = queries
            .iter()
            .flat_map(|&q| (0..n).map(move |t| (q, t)))
            .collect();
        let al = align(motifs, &pairs, opts)?;
        let s: Vec<Vec<f64>> = al
            .chunks(n)
            .map(|r| r.iter().map(|a| mean_score(a, motifs, opts)).collect())
            .collect();
        let nan = vec![vec![f64::NAN; n]; n];
        (s, nan.clone(), nan)
    };
    Ok(ComparisonMatrices {
        scores,
        pvalues: pvals,
        qvalues: qvals,
    })
}
