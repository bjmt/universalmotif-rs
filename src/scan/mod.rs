//! Scanning sequences for motif hits (`scan_sequences_lite()`).
//!
//! Each motif is turned into a PWM, a score threshold is found from the
//! P-value cutoff (see [`crate::pvalue`]), and every window of every sequence
//! scoring at least the threshold is reported, on both strands by default
//! for DNA and RNA and on the forward strand for other alphabets.
//! Scores use the same integer grid as the P-value computations, so a hit's
//! P-value is always within the cutoff.
//!
//! Sequences are packed at two bits per letter (see [`crate::packed`]), and
//! each group of four motif positions is scored with one lookup in a table of
//! all 256 four-letter combinations. Other alphabets are coded at a byte
//! per letter and scored two positions per lookup. Work is split into chunks
//! of the sequences, so a single long sequence still uses every thread, and
//! motifs are taken in batches so that the P-value tables held at once stay
//! small.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::cdf::{ScoreDistribution, ScoreTail, grid_value, score_grid_ceiling};
use crate::convert::convert_type;
use crate::dedup::{Interval, dedup};
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::ops::normalize;
use crate::packed::{CodedSeq, EncodedSequences, PackedSeq, Seqs};
use crate::pvalue::{scoring_matrix, with_threads};

mod gapped;

use gapped::{GapShape, GappedScorer, PieceScorer};

/// Which value decides between overlapping hits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OverlapPriority {
    /// Keep the hit with the lowest P-value.
    #[default]
    Pvalue,
    /// Keep the hit with the highest score.
    Score,
}

/// Removal of overlapping hits (see [`crate::dedup`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct NoOverlaps {
    pub by: OverlapPriority,
    /// Compare hits on different strands too.
    pub ignore_strand: bool,
    /// Compare hits of different motifs too.
    pub ignore_motif: bool,
}

/// Options for [`scan_sequences`].
#[derive(Debug, Clone, PartialEq)]
pub struct ScanOptions {
    /// P-value cutoff, in (0, 1).
    pub pvalue: f64,
    /// Also scan the reverse strand.
    pub rc: bool,
    /// Remove overlapping hits.
    pub no_overlaps: Option<NoOverlaps>,
    /// Threads to use; 0 means all available.
    pub nthreads: usize,
    /// Scan gapped motifs with their gaps (otherwise their gaps are
    /// ignored).
    pub use_gaps: bool,
    /// Bonferroni-correct gapped motifs' P-values for the number of
    /// combinations of gap sizes tried: the score threshold comes from the
    /// cutoff divided by that number, and P-values are multiplied by it (up
    /// to 1). Off by default, as in the R package.
    pub correct_gaps: bool,
}

impl Default for ScanOptions {
    fn default() -> Self {
        ScanOptions {
            pvalue: 1e-4,
            rc: true,
            no_overlaps: None,
            nthreads: 1,
            use_gaps: true,
            correct_gaps: false,
        }
    }
}

/// A motif hit.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    /// Index of the motif.
    pub motif: usize,
    /// Index of the sequence.
    pub sequence: usize,
    /// First position of the hit on the forward strand (1-based).
    pub start: usize,
    /// Last position of the hit on the forward strand (1-based, inclusive).
    pub end: usize,
    pub strand: Strand,
    pub score: f64,
    /// Score as a percentage of the motif's maximum score.
    pub score_pct: f64,
    /// The matched sequence, reverse complemented for `-` strand hits.
    pub matched: String,
    pub pvalue: f64,
}

/// Hits as columns, one entry per hit, ordered by motif, sequence and start
/// (with `+` strand hits first). Indices and positions are `i64` so the
/// columns can be handed to Arrow as they are. The matched letters are not
/// stored; [`HitTable::matched`] reads them from the sequences (and
/// [`HitTable::matched_encoded`] from encoded ones).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct HitTable {
    /// Index of the motif.
    pub motif: Vec<i64>,
    /// Index of the sequence.
    pub sequence: Vec<i64>,
    /// First position on the forward strand (1-based).
    pub start: Vec<i64>,
    /// Last position on the forward strand (1-based, inclusive).
    pub end: Vec<i64>,
    /// Whether the hit is on the `-` strand.
    pub minus: Vec<bool>,
    pub score: Vec<f64>,
    /// Score as a percentage of the motif's maximum score.
    pub score_pct: Vec<f64>,
    pub pvalue: Vec<f64>,
    rna: bool,
    /// Whether matched letters are reported in upper case (all but custom
    /// alphabets).
    upper: bool,
    /// Each hit's combination of gap sizes (0 for motifs without gaps).
    combo: Vec<u16>,
    /// The gaps of each motif scanned with them.
    gaps: Vec<Option<GapShape>>,
}

impl HitTable {
    pub fn len(&self) -> usize {
        self.motif.len()
    }

    pub fn is_empty(&self) -> bool {
        self.motif.is_empty()
    }

    /// The size of each gap in hit `i`, if its motif was scanned with gaps.
    pub fn gap_sizes(&self, i: usize) -> Option<Vec<usize>> {
        self.gaps[self.motif[i] as usize]
            .as_ref()
            .map(|g| g.gap_sizes(self.combo[i]))
    }

    /// Show the gaps of hit `i` in its matched letters as dots.
    fn mark_gaps(&self, i: usize, out: &mut [u8]) {
        if let Some(g) = &self.gaps[self.motif[i] as usize] {
            g.mark_gaps(self.combo[i], out);
        }
    }

    /// Write the letters matched by hit `i` (upper case, reverse
    /// complemented for `-` strand hits, and with dots in the gaps of
    /// gapped motifs) into `out`, replacing its contents.
    pub fn matched_into<S: AsRef<str>>(&self, i: usize, sequences: &[S], out: &mut Vec<u8>) {
        let seq = sequences[self.sequence[i] as usize].as_ref().as_bytes();
        let bytes = &seq[self.start[i] as usize - 1..self.end[i] as usize];
        out.clear();
        if self.minus[i] {
            out.extend(
                bytes
                    .iter()
                    .rev()
                    .map(|c| complement(c.to_ascii_uppercase(), self.rna)),
            );
        } else if self.upper {
            out.extend(bytes.iter().map(u8::to_ascii_uppercase));
        } else {
            out.extend_from_slice(bytes);
        }
        self.mark_gaps(i, out);
    }

    /// [`HitTable::matched_into`] for hits found by [`scan_encoded`],
    /// reading the letters from the encoded sequences.
    pub fn matched_encoded_into(&self, i: usize, seqs: &EncodedSequences, out: &mut Vec<u8>) {
        out.clear();
        let (a, b) = (self.start[i] as usize - 1, self.end[i] as usize);
        seqs.decode_into(self.sequence[i] as usize, a, b, out);
        if self.minus[i] {
            out.reverse();
            out.iter_mut().for_each(|c| *c = complement(*c, self.rna));
        }
        self.mark_gaps(i, out);
    }

    /// The letters matched by hit `i` (see
    /// [`HitTable::matched_encoded_into`]).
    pub fn matched_encoded(&self, i: usize, seqs: &EncodedSequences) -> String {
        let mut out = Vec::new();
        self.matched_encoded_into(i, seqs, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    /// The letters matched by hit `i` (see [`HitTable::matched_into`]).
    pub fn matched<S: AsRef<str>>(&self, i: usize, sequences: &[S]) -> String {
        let mut out = Vec::new();
        self.matched_into(i, sequences, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    /// One [`Hit`] per row.
    pub fn to_hits<S: AsRef<str>>(&self, sequences: &[S]) -> Vec<Hit> {
        (0..self.len())
            .map(|i| Hit {
                motif: self.motif[i] as usize,
                sequence: self.sequence[i] as usize,
                start: self.start[i] as usize,
                end: self.end[i] as usize,
                strand: if self.minus[i] {
                    Strand::Minus
                } else {
                    Strand::Plus
                },
                score: self.score[i],
                score_pct: self.score_pct[i],
                matched: self.matched(i, sequences),
                pvalue: self.pvalue[i],
            })
            .collect()
    }
}

fn complement(c: u8, rna: bool) -> u8 {
    match c {
        b'A' => {
            if rna {
                b'U'
            } else {
                b'T'
            }
        }
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

/// Reverse complement a DNA or RNA sequence, IUPAC codes included. The
/// result is upper case.
pub fn reverse_complement_seq(seq: &str, rna: bool) -> String {
    rc_bytes(seq.as_bytes(), rna)
}

fn rc_bytes(seq: &[u8], rna: bool) -> String {
    let out: Vec<u8> = seq
        .iter()
        .rev()
        .map(|c| complement(c.to_ascii_uppercase(), rna))
        .collect();
    String::from_utf8_lossy(&out).into_owned()
}

/// Integer scores for each group of four positions, indexed by the group's
/// letters packed into one byte (see [`PackedSeq::quad`]). A last group of
/// fewer than four positions ignores the extra letters.
struct Scorer {
    groups: Vec<[i32; 256]>,
}

impl Scorer {
    fn new(m: &crate::matrix::Matrix) -> Scorer {
        let cols: Vec<Vec<i32>> = m
            .columns()
            .map(|c| c.iter().map(|&v| grid_value(v)).collect())
            .collect();
        let groups = cols
            .chunks(4)
            .map(|g| {
                let mut t = [0i32; 256];
                for (code, slot) in t.iter_mut().enumerate() {
                    *slot = g
                        .iter()
                        .enumerate()
                        .map(|(j, col)| col[(code >> (6 - 2 * j)) & 3])
                        .sum();
                }
                t
            })
            .collect();
        Scorer { groups }
    }
}

/// Integer scores for groups of `span` positions (two, or one for alphabets
/// of more than 64 letters) for alphabets other than DNA and RNA, indexed by
/// the group's letter codes as a base-n number. A last group of fewer
/// positions ignores the extra letters.
struct CodedScorer {
    groups: Vec<Vec<i32>>,
}

impl CodedScorer {
    fn new(m: &crate::matrix::Matrix) -> CodedScorer {
        let n = m.nrow();
        let span = if n <= 64 { 2 } else { 1 };
        let cols: Vec<Vec<i32>> = m
            .columns()
            .map(|c| c.iter().map(|&v| grid_value(v)).collect())
            .collect();
        let groups = cols
            .chunks(span)
            .map(|g| {
                (0..n.pow(span as u32))
                    .map(|code| {
                        g.iter()
                            .enumerate()
                            .map(|(j, col)| col[(code / n.pow((span - 1 - j) as u32)) % n])
                            .sum()
                    })
                    .collect()
            })
            .collect();
        CodedScorer { groups }
    }
}

/// How each motif is scored: DNA and RNA on packed sequences, one scorer
/// per strand; other alphabets letter by letter, forward strand only;
/// gapped motifs piece by piece (see [`gapped`]).
enum Scorers {
    Packed(Vec<(Strand, Scorer)>),
    Coded(CodedScorer),
    Gapped(GappedScorer),
}

/// Window starts per unit of work; the decoded letters stay in cache.
const CHUNK: usize = 1 << 16;

/// Grid points of P-value tails held at once while scanning (32 MB).
const TAIL_BUDGET: usize = 1 << 22;

/// Motifs prepared in parallel at a time while filling a batch.
const PREPARE_GROUP: usize = 32;

/// Letters of sequence read at a time by [`scan_fasta`] (a 16 MB batch for
/// DNA and RNA).
const STREAM_BATCH: usize = 1 << 26;

/// Buffers [`scan_fasta`] keeps for reuse.
const SPARE_BUFFERS: usize = 4;

/// What scanning needs to know about one motif.
struct Prepared {
    width: usize,
    /// The widest window the motif covers (its width unless it is gapped).
    span: usize,
    /// What P-values are multiplied by (see [`ScanOptions::correct_gaps`]).
    pvalue_scale: f64,
    min_score: i32,
    max_score: f64,
    scorers: Scorers,
    /// P-values of the scores at or above the threshold, unless they were
    /// left out to save memory.
    tail: Option<ScoreTail>,
}

/// What the hit table needs to know about a motif.
struct Shape {
    width: usize,
    max_score: f64,
    gaps: Option<GapShape>,
}

impl Shape {
    fn of(p: &Prepared) -> Shape {
        Shape {
            width: p.width,
            max_score: p.max_score,
            gaps: match &p.scorers {
                Scorers::Gapped(g) => Some(g.shape.clone()),
                _ => None,
            },
        }
    }

    /// The window a hit covers, with its combination of gap sizes.
    fn span(&self, combo: u16) -> usize {
        self.gaps.as_ref().map_or(self.width, |g| g.width(combo))
    }
}

/// A hit before it is turned into a [`Hit`].
struct RawHit {
    motif: u32,
    sequence: u32,
    start: u32,
    score: i32,
    /// Position in the order the hits were found, for [`scan_fasta`].
    id: u32,
    minus: bool,
    /// The combination of gap sizes, for gapped motifs.
    combo: u16,
    /// Filled in once the motif's tail is at hand.
    pvalue: f64,
}

/// Scan the windows starting in `starts` of one stretch of a sequence free
/// of unknown letters (ending at `seg_end`) with every motif. The packed
/// letters are decoded once into one byte per position, holding the four
/// letters that start there.
#[allow(clippy::too_many_arguments)]
// Kept out of line: inlined into the work loop, it compiles to a slower
// scanning loop.
#[inline(never)]
fn scan_chunk(
    prepared: &[Option<Prepared>],
    seq: &PackedSeq,
    si: usize,
    starts: (usize, usize),
    seg_end: usize,
    quads: &mut Vec<u8>,
    pieces: &mut Vec<Vec<i32>>,
    out: &mut Vec<RawHit>,
) {
    let (c, end) = starts;
    let max_w = prepared.iter().flatten().map(|p| p.span).max().unwrap_or(0);
    quads.clear();
    quads.extend((c..(end + max_w).min(seg_end + 4)).map(|i| seq.quad(i)));
    for (mi, p) in prepared.iter().enumerate() {
        let Some(p) = p else { continue };
        if seg_end < c + p.width {
            continue;
        }
        let last = end.min(seg_end + 1 - p.width);
        let Scorers::Packed(scorers) = &p.scorers else {
            continue;
        };
        for (strand, scorer) in scorers {
            for i in c..last {
                let q = &quads[i - c..];
                let mut s = 0i32;
                for (g, t) in scorer.groups.iter().enumerate() {
                    s += t[q[4 * g] as usize];
                }
                if s >= p.min_score {
                    out.push(RawHit {
                        motif: mi as u32,
                        sequence: si as u32,
                        start: i as u32,
                        minus: *strand == Strand::Minus,
                        score: s,
                        id: 0,
                        combo: 0,
                        pvalue: f64::NAN,
                    });
                }
            }
        }
    }
    if prepared
        .iter()
        .flatten()
        .any(|p| matches!(p.scorers, Scorers::Gapped(_)))
    {
        let stop = (end + max_w).min(seg_end);
        scan_chunk_gapped(prepared, quads, si, starts, stop, seg_end, pieces, out);
    }
}

/// [`scan_chunk`] for alphabets other than DNA and RNA: the letter codes of
/// each group of `span` positions are combined once per position.
#[allow(clippy::too_many_arguments)]
// Kept out of line: inlined into the work loop, it compiles to a slower
// scanning loop.
#[inline(never)]
fn scan_chunk_coded(
    prepared: &[Option<Prepared>],
    seq: &CodedSeq,
    si: usize,
    starts: (usize, usize),
    seg_end: usize,
    n: usize,
    span: usize,
    keys: &mut Vec<u16>,
    pieces: &mut Vec<Vec<i32>>,
    out: &mut Vec<RawHit>,
) {
    let (c, end) = starts;
    let max_w = prepared.iter().flatten().map(|p| p.span).max().unwrap_or(0);
    keys.clear();
    // The last letter any window here can use is at seg_end - 1.
    let stop = (end + max_w).min(seg_end);
    if span == 2 {
        keys.extend((c..stop).map(|i| (seq.get(i) as usize * n + seq.get(i + 1) as usize) as u16));
    } else {
        keys.extend((c..stop).map(|i| seq.get(i) as u16));
    }
    for (mi, p) in prepared.iter().enumerate() {
        let Some(p) = p else { continue };
        let Scorers::Coded(scorer) = &p.scorers else {
            continue;
        };
        if seg_end < c + p.width {
            continue;
        }
        let last = end.min(seg_end + 1 - p.width);
        for i in c..last {
            let k = &keys[i - c..];
            let mut s = 0i32;
            for (g, t) in scorer.groups.iter().enumerate() {
                s += t[k[span * g] as usize];
            }
            if s >= p.min_score {
                out.push(RawHit {
                    motif: mi as u32,
                    sequence: si as u32,
                    start: i as u32,
                    minus: false,
                    score: s,
                    id: 0,
                    combo: 0,
                    pvalue: f64::NAN,
                });
            }
        }
    }
    if prepared
        .iter()
        .flatten()
        .any(|p| matches!(p.scorers, Scorers::Gapped(_)))
    {
        scan_chunk_coded_gapped(prepared, keys, si, starts, stop, seg_end, span, pieces, out);
    }
}

/// The gapped motifs' part of [`scan_chunk`], kept apart so the loop over
/// other motifs compiles as it would without them.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn scan_chunk_gapped(
    prepared: &[Option<Prepared>],
    quads: &[u8],
    si: usize,
    starts: (usize, usize),
    stop: usize,
    seg_end: usize,
    pieces: &mut Vec<Vec<i32>>,
    out: &mut Vec<RawHit>,
) {
    let c = starts.0;
    let score = |piece: &PieceScorer, at: usize| match piece {
        PieceScorer::Packed(s) => {
            let q = &quads[at - c..];
            s.groups
                .iter()
                .enumerate()
                .map(|(g, t)| t[q[4 * g] as usize])
                .sum()
        }
        PieceScorer::Coded(_) => unreachable!("DNA and RNA pieces are packed"),
    };
    for (mi, p) in prepared.iter().enumerate() {
        if let Some(Prepared {
            scorers: Scorers::Gapped(g),
            min_score,
            ..
        }) = p
        {
            g.scan(
                mi, si, starts, stop, seg_end, *min_score, pieces, score, out,
            );
        }
    }
}

/// The gapped motifs' part of [`scan_chunk_coded`] (see
/// [`scan_chunk_gapped`]).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn scan_chunk_coded_gapped(
    prepared: &[Option<Prepared>],
    keys: &[u16],
    si: usize,
    starts: (usize, usize),
    stop: usize,
    seg_end: usize,
    span: usize,
    pieces: &mut Vec<Vec<i32>>,
    out: &mut Vec<RawHit>,
) {
    let c = starts.0;
    let score = |piece: &PieceScorer, at: usize| match piece {
        PieceScorer::Coded(s) => {
            let k = &keys[at - c..];
            s.groups
                .iter()
                .enumerate()
                .map(|(g, t)| t[k[span * g] as usize])
                .sum()
        }
        PieceScorer::Packed(_) => unreachable!("other alphabets are coded"),
    };
    for (mi, p) in prepared.iter().enumerate() {
        if let Some(Prepared {
            scorers: Scorers::Gapped(g),
            min_score,
            ..
        }) = p
        {
            g.scan(
                mi, si, starts, stop, seg_end, *min_score, pieces, score, out,
            );
        }
    }
}

/// The alphabet of the motifs, after checking the motifs and options.
fn check_scan(motifs: &[Motif], opts: &ScanOptions) -> Result<Alphabet> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs to scan with");
    };
    if !(opts.pvalue > 0.0 && opts.pvalue < 1.0) {
        return arg_err("the P-value cutoff must be between 0 and 1");
    }
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    Ok(first.alphabet.clone())
}

/// A motif's scoring tables and threshold, and its P-value tail if
/// `keep_tail`; `None` if no window can reach the threshold.
fn prepare_motif(m: &Motif, opts: &ScanOptions, keep_tail: bool) -> Result<Option<Prepared>> {
    let mut pwm = convert_type(m, MotifType::Pwm, None, false)?;
    if pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
        pwm = convert_type(&normalize(m)?, MotifType::Pwm, None, false)?;
    }
    let sm = scoring_matrix(&pwm, None)?;
    let gaps = if opts.use_gaps && m.gap.isgapped && !m.gap.gaploc.is_empty() {
        Some(GapShape::from_motif(m)?)
    } else {
        None
    };
    let pvalue_scale = match &gaps {
        Some(g) if opts.correct_gaps => g.combos() as f64,
        _ => 1.0,
    };
    let dist = ScoreDistribution::new(&sm.pwm, &sm.bkg);
    let threshold = dist.score(opts.pvalue / pvalue_scale);
    let min = score_grid_ceiling(threshold);
    if min >= i32::MAX as f64 {
        return Ok(None);
    }
    let nucleotide = m.alphabet.is_nucleotide();
    let scorers = if let Some(shape) = gaps {
        Scorers::Gapped(GappedScorer::new(&sm.pwm, shape, nucleotide, opts.rc))
    } else if nucleotide {
        let mut s = vec![(Strand::Plus, Scorer::new(&sm.pwm))];
        if opts.rc {
            s.push((Strand::Minus, Scorer::new(&sm.pwm.reversed())));
        }
        Scorers::Packed(s)
    } else {
        Scorers::Coded(CodedScorer::new(&sm.pwm))
    };
    let span = match &scorers {
        Scorers::Gapped(g) => g.shape.max_width(),
        _ => sm.pwm.ncol(),
    };
    Ok(Some(Prepared {
        width: sm.pwm.ncol(),
        span,
        pvalue_scale,
        min_score: min as i32,
        max_score: sm
            .pwm
            .columns()
            .map(|c| c.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .sum(),
        scorers,
        tail: keep_tail.then(|| dist.tail(threshold)),
    }))
}

/// Units of work: chunks of window starts within the stretches of each
/// sequence free of unknown letters, as (sequence, starts, stretch end).
fn work_items(seqs: &EncodedSequences) -> Vec<(usize, (usize, usize), usize)> {
    let segments: Vec<Vec<(usize, usize)>> = match seqs.seqs() {
        Seqs::Packed(v) => v.iter().map(|s| s.known_segments().collect()).collect(),
        Seqs::Coded(v) => v.iter().map(|s| s.known_segments().collect()).collect(),
    };
    let mut work = Vec::new();
    for (si, segs) in segments.iter().enumerate() {
        for &(a, b) in segs {
            let mut c = a;
            while c < b {
                let end = (c + CHUNK).min(b);
                work.push((si, (c, end), b));
                c = end;
            }
        }
    }
    work
}

/// Every hit of the prepared motifs in `seqs`, in no particular order and
/// without P-values. Runs on the current thread pool.
fn scan_raw(
    prepared: &[Option<Prepared>],
    seqs: &EncodedSequences,
    work: &[(usize, (usize, usize), usize)],
) -> Vec<RawHit> {
    let n_letters = seqs.alphabet().len();
    let span = if n_letters <= 64 { 2 } else { 1 };
    work.par_iter()
        .fold(
            || (Vec::new(), Vec::new(), Vec::new(), Vec::new()),
            |(mut quads, mut keys, mut pieces, mut out), &(si, starts, seg_end)| {
                match seqs.seqs() {
                    Seqs::Packed(v) => scan_chunk(
                        prepared,
                        &v[si],
                        si,
                        starts,
                        seg_end,
                        &mut quads,
                        &mut pieces,
                        &mut out,
                    ),
                    Seqs::Coded(v) => scan_chunk_coded(
                        prepared,
                        &v[si],
                        si,
                        starts,
                        seg_end,
                        n_letters,
                        span,
                        &mut keys,
                        &mut pieces,
                        &mut out,
                    ),
                }
                (quads, keys, pieces, out)
            },
        )
        .map(|(_, _, _, out)| out)
        .reduce(Vec::new, |mut a, mut b| {
            a.append(&mut b);
            a
        })
}

/// Order the hits, drop overlapping ones if asked and build the table,
/// returning with it each row's [`RawHit::id`]. `shapes` holds each motif's
/// width and best possible score.
fn finish(
    mut raw: Vec<RawHit>,
    shapes: &[Option<Shape>],
    opts: &ScanOptions,
    alphabet: &Alphabet,
) -> (HitTable, Vec<u32>) {
    raw.sort_unstable_by_key(|r| (r.motif, r.sequence, r.start, r.minus, r.combo));
    if let Some(no) = opts.no_overlaps {
        let intervals: Vec<Interval> = raw
            .iter()
            .map(|r| {
                let w = shapes[r.motif as usize].as_ref().unwrap().span(r.combo) as i64;
                Interval {
                    sequence: r.sequence as usize,
                    motif: r.motif as usize,
                    strand: r.minus as i32,
                    start: r.start as i64 + 1,
                    end: r.start as i64 + w,
                    priority: match no.by {
                        OverlapPriority::Pvalue => -r.pvalue,
                        OverlapPriority::Score => r.score as f64,
                    },
                }
            })
            .collect();
        let keep = dedup(&intervals, no.ignore_strand, no.ignore_motif);
        let mut k = keep.iter();
        raw.retain(|_| *k.next().unwrap());
    }

    let n = raw.len();
    let mut t = HitTable {
        motif: Vec::with_capacity(n),
        sequence: Vec::with_capacity(n),
        start: Vec::with_capacity(n),
        end: Vec::with_capacity(n),
        minus: Vec::with_capacity(n),
        score: Vec::with_capacity(n),
        score_pct: Vec::with_capacity(n),
        pvalue: Vec::with_capacity(n),
        rna: *alphabet == Alphabet::Rna,
        upper: !matches!(alphabet, Alphabet::Custom(_)),
        combo: Vec::with_capacity(n),
        gaps: shapes
            .iter()
            .map(|s| s.as_ref().and_then(|s| s.gaps.clone()))
            .collect(),
    };
    for r in &raw {
        let shape = shapes[r.motif as usize].as_ref().unwrap();
        let (width, max_score) = (shape.span(r.combo), shape.max_score);
        let score = r.score as f64 / 1000.0;
        t.motif.push(r.motif as i64);
        t.sequence.push(r.sequence as i64);
        t.start.push(r.start as i64 + 1);
        t.end.push((r.start as usize + width) as i64);
        t.combo.push(r.combo);
        t.minus.push(r.minus);
        t.score.push(score);
        t.score_pct.push(100.0 * score / max_score);
        t.pvalue.push(r.pvalue);
    }
    (t, raw.iter().map(|r| r.id).collect())
}

/// Scan sequences for hits of the motifs, returning them as columns.
///
/// Letters of DNA, RNA and amino acid sequences are matched
/// case-insensitively, and letters of custom alphabets exactly; other
/// letters (such as `N` or `X`) never contribute to a hit. Only DNA and RNA
/// have a reverse strand. Sequences shorter than a motif simply have no hits
/// for it. Hits are ordered by motif, sequence and start, with `+` strand
/// hits first.
///
/// The sequences are encoded for scanning first (see [`crate::packed`]);
/// to scan the same sequences more than once, or to scan a FASTA file
/// without holding it as text, see [`scan_encoded`] and [`scan_fasta`].
pub fn scan_table<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    opts: &ScanOptions,
) -> Result<HitTable> {
    let alphabet = check_scan(motifs, opts)?;
    let seqs = EncodedSequences::new(sequences, &alphabet, opts.nthreads)?;
    scan_encoded(motifs, &seqs, opts)
}

/// [`scan_table`] for sequences already encoded, which must be for the
/// motifs' alphabet. [`HitTable::matched_encoded`] gives the matched
/// letters.
///
/// # Example
///
/// Pack a genome once and scan it repeatedly, at a quarter of its size as
/// text:
///
/// ```no_run
/// use universalmotif::Alphabet;
/// use universalmotif::io::{fasta, meme};
/// use universalmotif::scan::{ScanOptions, scan_encoded};
///
/// let motifs = meme::read(&std::fs::read_to_string("motifs.meme")?, &Default::default())?;
/// let (names, genome) = fasta::read_encoded("genome.fa.gz", &Alphabet::Dna, 0)?;
/// let opts = ScanOptions { pvalue: 1e-5, nthreads: 8, ..Default::default() };
/// let hits = scan_encoded(&motifs, &genome, &opts)?;
/// for i in 0..hits.len() {
///     let chrom = &names[hits.sequence[i] as usize];
///     println!("{chrom}\t{}\t{}\t{}", hits.start[i], hits.end[i], hits.matched_encoded(i, &genome));
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn scan_encoded(
    motifs: &[Motif],
    seqs: &EncodedSequences,
    opts: &ScanOptions,
) -> Result<HitTable> {
    let alphabet = check_scan(motifs, opts)?;
    if *seqs.alphabet() != alphabet {
        return arg_err(format!(
            "the sequences were encoded as {} but the motifs are {alphabet}",
            seqs.alphabet()
        ));
    }
    let mut raw: Vec<RawHit> = Vec::new();
    // Each motif's width and best possible score, or None if it cannot hit.
    let mut shapes: Vec<Option<Shape>> = Vec::with_capacity(motifs.len());
    with_threads(opts.nthreads, || -> Result<()> {
        let work = work_items(seqs);
        // Motifs are scanned in batches whose P-value tails fit in
        // TAIL_BUDGET, since a wide score range (as amino acid PWMs have)
        // makes for long tails. The hits of each batch get their P-values
        // before its tails are dropped.
        let mut base = 0;
        while base < motifs.len() {
            let mut batch: Vec<Option<Prepared>> = Vec::new();
            let mut size = 0;
            while base + batch.len() < motifs.len() && size < TAIL_BUDGET {
                let lo = base + batch.len();
                let hi = (lo + PREPARE_GROUP).min(motifs.len());
                let group = motifs[lo..hi]
                    .par_iter()
                    .map(|m| prepare_motif(m, opts, true))
                    .collect::<Result<Vec<Option<Prepared>>>>()?;
                size += group
                    .iter()
                    .flatten()
                    .map(|p| p.tail.as_ref().map_or(0, ScoreTail::len))
                    .sum::<usize>();
                batch.extend(group);
            }
            let mut hits = scan_raw(&batch, seqs, &work);
            hits.par_iter_mut().for_each(|r| {
                let p = batch[r.motif as usize].as_ref().unwrap();
                let raw = p.tail.as_ref().unwrap().pvalue(r.score as f64 / 1000.0);
                r.pvalue = (raw * p.pvalue_scale).min(1.0);
                r.motif += base as u32;
            });
            raw.append(&mut hits);
            shapes.extend(batch.iter().map(|p| p.as_ref().map(Shape::of)));
            base += batch.len();
        }
        Ok(())
    })??;
    Ok(finish(raw, &shapes, opts, &alphabet).0)
}

/// Hits from [`scan_fasta`], with the letters each one matched (since the
/// sequences themselves are not kept).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FastaHits {
    /// The name of every sequence in the file, in order.
    pub names: Vec<String>,
    /// The length of every sequence in the file.
    pub lengths: Vec<usize>,
    pub table: HitTable,
    letters: Vec<u8>,
    /// Where each row's letters start in `letters`.
    offsets: Vec<usize>,
}

impl FastaHits {
    /// The number of hits.
    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }

    /// The letters matched by hit `i` (as [`HitTable::matched_into`]).
    pub fn matched(&self, i: usize) -> &str {
        let t = &self.table;
        let w = (t.end[i] - t.start[i] + 1) as usize;
        std::str::from_utf8(&self.letters[self.offsets[i]..self.offsets[i] + w]).unwrap()
    }
}

/// [`scan_table`] for the sequences of a FASTA file (which may be gzipped,
/// bgzip included), read a batch at a time on a thread of its own while the
/// last batch is scanned, so the file is never held in memory whole.
/// Sequence indices follow the file's records.
///
/// # Example
///
/// ```no_run
/// use universalmotif::io::meme;
/// use universalmotif::scan::{ScanOptions, scan_fasta};
///
/// let motifs = meme::read(&std::fs::read_to_string("motifs.meme")?, &Default::default())?;
/// let opts = ScanOptions { pvalue: 1e-5, nthreads: 8, ..Default::default() };
/// let hits = scan_fasta(&motifs, "genome.fa.gz", &opts)?;
/// for i in 0..hits.len() {
///     let t = &hits.table;
///     println!("{}\t{}\t{}", hits.names[t.sequence[i] as usize], t.start[i], hits.matched(i));
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn scan_fasta(
    motifs: &[Motif],
    path: impl AsRef<std::path::Path>,
    opts: &ScanOptions,
) -> Result<FastaHits> {
    scan_fasta_with(motifs, path.as_ref(), opts, STREAM_BATCH, TAIL_BUDGET)
}

/// [`scan_fasta`] with the batch size and tail budget given.
fn scan_fasta_with(
    motifs: &[Motif],
    path: &std::path::Path,
    opts: &ScanOptions,
    batch_size: usize,
    tail_budget: usize,
) -> Result<FastaHits> {
    let alphabet = check_scan(motifs, opts)?;
    let reader = crate::io::fasta::open(path)?;
    let mut out = FastaHits::default();
    let mut raw: Vec<RawHit> = Vec::new();
    // Each found hit's letters start at letter_starts[id].
    let mut letters: Vec<u8> = Vec::new();
    let mut letter_starts: Vec<usize> = Vec::new();

    // Unlike scan_encoded, every motif is needed for every batch, so tails
    // beyond the budget are dropped and rebuilt at the end instead.
    let prepared: Vec<Option<Prepared>> = with_threads(opts.nthreads, || {
        let mut all: Vec<Option<Prepared>> = Vec::with_capacity(motifs.len());
        let mut size = 0;
        for group in motifs.chunks(PREPARE_GROUP) {
            let group = group
                .par_iter()
                .map(|m| prepare_motif(m, opts, true))
                .collect::<Result<Vec<Option<Prepared>>>>()?;
            for mut p in group {
                if let Some(p) = p.as_mut() {
                    size += p.tail.as_ref().map_or(0, ScoreTail::len);
                    if size > tail_budget {
                        p.tail = None;
                    }
                }
                all.push(p);
            }
        }
        Ok::<_, Error>(all)
    })??;

    let shapes: Vec<Option<Shape>> = prepared.iter().map(|p| p.as_ref().map(Shape::of)).collect();

    let (tx, rx) = std::sync::mpsc::sync_channel::<(Vec<String>, EncodedSequences)>(1);
    // Scanned batches hand their buffers back to the reader. Allocating
    // fresh ones instead leaves freed chromosome-sized blocks to the
    // allocator, which may keep them (macOS does), so memory would grow
    // with the file rather than with the batch.
    let (spare_tx, spare_rx) = std::sync::mpsc::channel::<Vec<u8>>();
    let (scanned, read) = std::thread::scope(|s| {
        // The reader owns the sender, so the batches end when it does.
        let alph = alphabet.clone();
        let reader = s.spawn(move || {
            let mut pool: Vec<Vec<u8>> = Vec::new();
            let spare = || {
                pool.extend(spare_rx.try_iter());
                // Keep the largest few; any more would only hold memory.
                pool.sort_unstable_by_key(Vec::capacity);
                let excess = pool.len().saturating_sub(SPARE_BUFFERS);
                pool.drain(..excess);
                pool.pop()
            };
            crate::io::fasta::read_encoded_batches_reusing(
                reader,
                &alph,
                0,
                batch_size,
                spare,
                |n, b| {
                    tx.send((n, b))
                        .map_err(|_| Error::Argument("scanning stopped".into()))
                },
            )
        });
        let scanned = with_threads(opts.nthreads, || {
            for (names, batch) in rx {
                let base = out.names.len();
                let work = work_items(&batch);
                let mut hits = scan_raw(&prepared, &batch, &work);
                for h in &mut hits {
                    let shape = shapes[h.motif as usize].as_ref().unwrap();
                    let start = letters.len();
                    let (a, b) = (h.start as usize, h.start as usize + shape.span(h.combo));
                    batch.decode_into(h.sequence as usize, a, b, &mut letters);
                    if h.minus {
                        let rna = alphabet == Alphabet::Rna;
                        letters[start..].reverse();
                        letters[start..]
                            .iter_mut()
                            .for_each(|c| *c = complement(*c, rna));
                    }
                    if let Some(g) = &shape.gaps {
                        g.mark_gaps(h.combo, &mut letters[start..]);
                    }
                    h.id = letter_starts.len() as u32;
                    letter_starts.push(start);
                    h.sequence += base as u32;
                }
                raw.append(&mut hits);
                out.lengths
                    .extend((0..batch.len()).map(|i| batch.seq_len(i)));
                out.names.extend(names);
                for buf in batch.into_buffers() {
                    if buf.capacity() >= 1 << 20 {
                        // The reader is gone if this fails; nothing to do.
                        let _ = spare_tx.send(buf);
                    }
                }
            }
        });
        (scanned, reader.join().expect("the FASTA reader panicked"))
    });
    scanned?;
    read?;

    // P-values, motif by motif, rebuilding the tails that were dropped.
    raw.sort_unstable_by_key(|r| (r.motif, r.sequence, r.start, r.minus));
    let groups: Vec<&mut [RawHit]> = raw.chunk_by_mut(|a, b| a.motif == b.motif).collect();
    with_threads(opts.nthreads, || {
        groups.into_par_iter().try_for_each(|g| -> Result<()> {
            let m = g[0].motif as usize;
            let scale = prepared[m].as_ref().unwrap().pvalue_scale;
            let rebuilt;
            let tail = match &prepared[m].as_ref().unwrap().tail {
                Some(t) => t,
                None => {
                    rebuilt = prepare_motif(&motifs[m], opts, true)?
                        .unwrap()
                        .tail
                        .unwrap();
                    &rebuilt
                }
            };
            for r in g {
                r.pvalue = (tail.pvalue(r.score as f64 / 1000.0) * scale).min(1.0);
            }
            Ok(())
        })
    })??;

    let (table, ids) = finish(raw, &shapes, opts, &alphabet);
    let mut kept = Vec::new();
    let mut offsets = Vec::with_capacity(ids.len());
    for (i, &id) in ids.iter().enumerate() {
        let w = (table.end[i] - table.start[i] + 1) as usize;
        let a = letter_starts[id as usize];
        offsets.push(kept.len());
        kept.extend_from_slice(&letters[a..a + w]);
    }
    out.table = table;
    out.letters = kept;
    out.offsets = offsets;
    Ok(out)
}

/// Scan sequences for hits of the motifs, one [`Hit`] per hit
/// (see [`scan_table`], which this wraps, for the details). For many hits,
/// [`scan_table`] takes far less memory.
pub fn scan_sequences<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    opts: &ScanOptions,
) -> Result<Vec<Hit>> {
    Ok(scan_table(motifs, sequences, opts)?.to_hits(sequences))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::create::{CreateOptions, from_strings};
    use crate::shuffle::create_sequences;

    /// Every window scored the plain way, for comparison with the chunked,
    /// packed scanner.
    fn naive(motif: &Motif, seq: &str, pvalue: f64) -> Vec<(usize, bool, i32)> {
        let sm = scoring_matrix(
            &convert_type(motif, MotifType::Pwm, None, false).unwrap(),
            None,
        )
        .unwrap();
        let min = score_grid_ceiling(ScoreDistribution::new(&sm.pwm, &sm.bkg).score(pvalue)) as i32;
        let index = |b: u8| b"ACGT".iter().position(|&l| l == b.to_ascii_uppercase());
        let w = sm.pwm.ncol();
        let mut out = Vec::new();
        for (minus, m) in [(false, sm.pwm.clone()), (true, sm.pwm.reversed())] {
            for i in 0..=seq.len() - w {
                let s: Option<i32> = (0..w)
                    .map(|j| index(seq.as_bytes()[i + j]).map(|l| grid_value(m.get(l, j))))
                    .sum();
                if let Some(s) = s.filter(|&s| s >= min) {
                    out.push((i + 1, minus, s));
                }
            }
        }
        out.sort();
        out
    }

    /// Forward-strand windows of any alphabet scored the plain way.
    fn naive_any(motif: &Motif, seq: &str, pvalue: f64, exact_case: bool) -> Vec<(usize, i32)> {
        let sm = scoring_matrix(
            &convert_type(motif, MotifType::Pwm, None, false).unwrap(),
            None,
        )
        .unwrap();
        let min = score_grid_ceiling(ScoreDistribution::new(&sm.pwm, &sm.bkg).score(pvalue)) as i32;
        let letters = motif.alphabet.letters();
        let index = |c: char| {
            letters.iter().position(|&l| {
                if exact_case {
                    l == c
                } else {
                    l == c.to_ascii_uppercase()
                }
            })
        };
        let chars: Vec<char> = seq.chars().collect();
        let w = sm.pwm.ncol();
        let mut out = Vec::new();
        for i in 0..=chars.len() - w {
            let s: Option<i32> = (0..w)
                .map(|j| index(chars[i + j]).map(|l| grid_value(sm.pwm.get(l, j))))
                .sum();
            if let Some(s) = s.filter(|&s| s >= min) {
                out.push((i + 1, s));
            }
        }
        out
    }

    #[test]
    fn other_alphabets_match_plain_scoring() {
        use crate::alphabet::Alphabet;
        let opts = CreateOptions {
            nsites: Some(20.0),
            pseudocount: 1.0,
            ..Default::default()
        };
        let aa_letters: Vec<char> = Alphabet::Aa.letters().to_vec();
        let aa = [
            from_strings(
                &["MKVLAWHE", "MKVLGWHE", "MRVLAWHD"],
                Some(Alphabet::Aa),
                &opts,
            )
            .unwrap(),
            from_strings(&["CPYQ"], Some(Alphabet::Aa), &opts).unwrap(),
            from_strings(&["GSTNDEKRHQAA"], Some(Alphabet::Aa), &opts).unwrap(),
        ];
        let mut seq = create_sequences(&aa_letters, 1, 80_000, None, 3, 1)
            .unwrap()
            .remove(0)
            .into_bytes();
        seq[100..103].copy_from_slice(b"XXX");
        seq[CHUNK - 2..CHUNK + 1].copy_from_slice(b"*X*");
        seq[5_000..6_000].make_ascii_lowercase();
        seq[7_000..7_008].copy_from_slice(b"MKVLAWHE");
        let seq = String::from_utf8(seq).unwrap();
        let sopts = ScanOptions {
            pvalue: 0.001,
            nthreads: 2,
            ..Default::default()
        };
        let hits = scan_table(&aa, &[&seq], &sopts).unwrap();
        assert!(hits.minus.iter().all(|m| !m), "proteins have one strand");
        for (mi, m) in aa.iter().enumerate() {
            let got: Vec<(usize, i32)> = (0..hits.len())
                .filter(|&i| hits.motif[i] == mi as i64)
                .map(|i| {
                    (
                        hits.start[i] as usize,
                        (hits.score[i] * 1000.0).round() as i32,
                    )
                })
                .collect();
            assert_eq!(got, naive_any(m, &seq, 0.001, false), "protein motif {mi}");
        }
        let planted = (0..hits.len()).find(|&i| hits.start[i] == 7_001).unwrap();
        assert_eq!(hits.matched(planted, &[&seq]), "MKVLAWHE");

        // A custom alphabet matches letters exactly, including case.
        let custom = Alphabet::Custom(vec!['a', 'b', 'c', 'D', 'e']);
        let m = from_strings(&["abcDe", "abcDa", "abDDe"], Some(custom.clone()), &opts).unwrap();
        let letters = custom.letters().to_vec();
        let mut cseq = create_sequences(&letters, 1, 30_000, None, 4, 1)
            .unwrap()
            .remove(0)
            .into_bytes();
        cseq[50..55].copy_from_slice(b"ABCDE");
        cseq[60..65].copy_from_slice(b"abcDe");
        let cseq = String::from_utf8(cseq).unwrap();
        let hits = scan_table(std::slice::from_ref(&m), &[&cseq], &sopts).unwrap();
        let got: Vec<(usize, i32)> = (0..hits.len())
            .map(|i| {
                (
                    hits.start[i] as usize,
                    (hits.score[i] * 1000.0).round() as i32,
                )
            })
            .collect();
        assert_eq!(got, naive_any(&m, &cseq, 0.001, true));
        let planted = (0..hits.len()).find(|&i| hits.start[i] == 61).unwrap();
        assert_eq!(hits.matched(planted, &[&cseq]), "abcDe");
    }

    /// Every window and combination of gap sizes of a gapped motif scored
    /// the plain way: `(start, minus, gap sizes, score, matched letters)`.
    fn naive_gapped(
        motif: &Motif,
        seq: &str,
        pvalue: f64,
        rc: bool,
    ) -> Vec<(usize, bool, Vec<usize>, i32, String)> {
        let sm = scoring_matrix(
            &convert_type(motif, MotifType::Pwm, None, false).unwrap(),
            None,
        )
        .unwrap();
        let min = score_grid_ceiling(ScoreDistribution::new(&sm.pwm, &sm.bkg).score(pvalue)) as i32;
        let letters: Vec<u8> = motif.alphabet.letters().iter().map(|&c| c as u8).collect();
        let shape = GapShape::from_motif(motif).unwrap();
        let rna = motif.alphabet == crate::alphabet::Alphabet::Rna;
        let bytes = seq.as_bytes();
        let mut out = Vec::new();
        for combo in 0..shape.combos() as u16 {
            let gaps = shape.gap_sizes(combo);
            let w = shape.width(combo);
            // Motif columns at each window offset, None in the gaps.
            let mut cols: Vec<Option<usize>> = Vec::new();
            let mut col = 0;
            for (piece, &width) in shape.piece_widths().iter().enumerate() {
                for _ in 0..width {
                    cols.push(Some(col));
                    col += 1;
                }
                if piece < gaps.len() {
                    cols.extend(std::iter::repeat_n(None, gaps[piece]));
                }
            }
            if bytes.len() < w {
                continue;
            }
            for i in 0..=bytes.len() - w {
                let window = &bytes[i..i + w];
                let index = |b: u8| letters.iter().position(|&l| l == b.to_ascii_uppercase());
                if window.iter().any(|&b| index(b).is_none()) {
                    continue;
                }
                for minus in [false, true] {
                    if minus && !rc {
                        continue;
                    }
                    let oriented: Vec<u8> = if minus {
                        window
                            .iter()
                            .rev()
                            .map(|&b| complement(b.to_ascii_uppercase(), rna))
                            .collect()
                    } else {
                        window.iter().map(u8::to_ascii_uppercase).collect()
                    };
                    let s: i32 = cols
                        .iter()
                        .zip(&oriented)
                        .filter_map(|(c, &b)| {
                            c.map(|c| grid_value(sm.pwm.get(index(b).unwrap(), c)))
                        })
                        .sum();
                    if s >= min {
                        let mut m = oriented.clone();
                        shape.mark_gaps(combo, &mut m);
                        out.push((i + 1, minus, gaps.clone(), s, String::from_utf8(m).unwrap()));
                    }
                }
            }
        }
        out.sort();
        out
    }

    #[test]
    fn gapped_scan_matches_plain_scoring() {
        use crate::alphabet::Alphabet;
        use crate::ops::add_gap;
        let copts = CreateOptions {
            nsites: Some(20.0),
            pseudocount: 1.0,
            ..Default::default()
        };
        let dna = [
            add_gap(
                &from_strings(&["TTGACTCA"], None, &copts).unwrap(),
                Some(&[3]),
                &[0],
                &[4],
            )
            .unwrap(),
            add_gap(
                &from_strings(&["CACGTGAAATT"], None, &copts).unwrap(),
                Some(&[2, 6]),
                &[1, 0],
                &[2, 5],
            )
            .unwrap(),
            add_gap(
                &from_strings(&["GGGCGGGA"], None, &copts).unwrap(),
                Some(&[1]),
                &[2],
                &[2],
            )
            .unwrap(),
        ];
        let mut seq = create_sequences(&['A', 'C', 'G', 'T'], 1, CHUNK + 4_000, None, 11, 1)
            .unwrap()
            .remove(0)
            .into_bytes();
        // Sites with each gap size, unknown letters inside a gap and across a
        // chunk boundary, and lower case.
        for (at, site) in [
            (100, &b"TTGACTCA"[..]),
            (200, b"TTGAAAACTCA"),
            (300, b"TTGNNACTCA"),
            (400, b"CAGCGTGAAAGGTT"),
            (CHUNK - 5, b"TTGACCCTCA"),
            (CHUNK + 100, b"tgagtgtgcaa"),
        ] {
            seq[at..at + site.len()].copy_from_slice(site);
        }
        seq[1_000..1_050].fill(b'N');
        let seq = String::from_utf8(seq).unwrap();
        for rc in [true, false] {
            let opts = ScanOptions {
                pvalue: 0.005,
                rc,
                nthreads: 2,
                ..Default::default()
            };
            let t = scan_table(&dna, &[&seq], &opts).unwrap();
            for (mi, m) in dna.iter().enumerate() {
                let mut got: Vec<_> = (0..t.len())
                    .filter(|&i| t.motif[i] == mi as i64)
                    .map(|i| {
                        let w = (t.end[i] - t.start[i] + 1) as usize;
                        assert_eq!(w, GapShape::from_motif(m).unwrap().width(t.combo[i]));
                        (
                            t.start[i] as usize,
                            t.minus[i],
                            t.gap_sizes(i).unwrap(),
                            (t.score[i] * 1000.0).round() as i32,
                            t.matched(i, &[&seq]),
                        )
                    })
                    .collect();
                got.sort();
                let want = naive_gapped(m, &seq, 0.005, rc);
                assert!(want.len() > 20, "motif {mi}: only {} hits", want.len());
                assert_eq!(got, want, "motif {mi}, rc {rc}");
            }
        }
        // Corrected for the number of gap sizes: the hits whose P-value is
        // within the cutoff divided by that number, their P-values scaled.
        let base = ScanOptions {
            pvalue: 0.005,
            nthreads: 2,
            ..Default::default()
        };
        let raw = scan_table(&dna, &[&seq], &base).unwrap();
        let corrected = ScanOptions {
            correct_gaps: true,
            ..base.clone()
        };
        let adj = scan_table(&dna, &[&seq], &corrected).unwrap();
        assert!(!adj.is_empty() && adj.len() < raw.len());
        let key = |t: &HitTable, i: usize| (t.motif[i], t.start[i], t.minus[i], t.combo[i]);
        let mut expected = Vec::new();
        for i in 0..raw.len() {
            let n = GapShape::from_motif(&dna[raw.motif[i] as usize])
                .unwrap()
                .combos() as f64;
            if raw.pvalue[i] <= 0.005 / n {
                expected.push((key(&raw, i), (raw.pvalue[i] * n).min(1.0)));
            }
        }
        let got: Vec<_> = (0..adj.len())
            .map(|i| (key(&adj, i), adj.pvalue[i]))
            .collect();
        assert_eq!(got, expected);
        assert!(adj.pvalue.iter().all(|&p| p <= 0.005));

        // Without use_gaps the motifs are scanned whole.
        let opts = ScanOptions {
            pvalue: 0.005,
            use_gaps: false,
            ..Default::default()
        };
        let ungapped: Vec<Motif> = dna.iter().map(|m| crate::ops::ungap(m, true)).collect();
        assert_eq!(
            scan_table(&dna, &[&seq], &opts).unwrap().start,
            scan_table(&ungapped, &[&seq], &opts).unwrap().start
        );

        // Other alphabets take the coded path, forward strand only.
        let aa = add_gap(
            &from_strings(&["MKVLAW"], Some(Alphabet::Aa), &copts).unwrap(),
            Some(&[2, 4]),
            &[0, 1],
            &[3, 2],
        )
        .unwrap();
        let aa_letters: Vec<char> = Alphabet::Aa.letters().to_vec();
        let mut aseq = create_sequences(&aa_letters, 1, 20_000, None, 12, 1)
            .unwrap()
            .remove(0);
        aseq.replace_range(500..509, "MKAAVLXAW");
        aseq.replace_range(700..708, "mkvvlqaw");
        let opts = ScanOptions {
            pvalue: 0.001,
            nthreads: 2,
            ..Default::default()
        };
        let t = scan_table(std::slice::from_ref(&aa), &[&aseq], &opts).unwrap();
        let mut got: Vec<_> = (0..t.len())
            .map(|i| {
                (
                    t.start[i] as usize,
                    t.minus[i],
                    t.gap_sizes(i).unwrap(),
                    (t.score[i] * 1000.0).round() as i32,
                    t.matched(i, &[&aseq]),
                )
            })
            .collect();
        got.sort();
        let want = naive_gapped(&aa, &aseq, 0.001, false);
        assert!(want.len() > 5);
        assert_eq!(got, want);
    }

    #[test]
    fn fasta_and_encoded_scans_match_scan_table() {
        use crate::alphabet::Alphabet;
        let copts = CreateOptions {
            nsites: Some(20.0),
            pseudocount: 1.0,
            ..Default::default()
        };
        let dna = [
            from_strings(&["TTGACGTCAAGG"], None, &copts).unwrap(),
            from_strings(&["CACGTG"], None, &copts).unwrap(),
            from_strings(&["WGATAR"], None, &copts).unwrap(),
            crate::ops::add_gap(
                &from_strings(&["TTGACTCA"], None, &copts).unwrap(),
                Some(&[3]),
                &[0],
                &[3],
            )
            .unwrap(),
        ];
        let aa = [
            from_strings(&["MKVLAW"], Some(Alphabet::Aa), &copts).unwrap(),
            from_strings(&["CPYQ"], Some(Alphabet::Aa), &copts).unwrap(),
        ];
        let mut dna_seqs = create_sequences(&['A', 'C', 'G', 'T'], 6, 3_000, None, 5, 1).unwrap();
        dna_seqs[1].replace_range(100..400, &"N".repeat(300));
        dna_seqs[2] = dna_seqs[2].to_ascii_lowercase();
        dna_seqs[4].truncate(3);
        let aa_letters: Vec<char> = Alphabet::Aa.letters().to_vec();
        let mut aa_seqs = create_sequences(&aa_letters, 4, 2_000, None, 6, 1).unwrap();
        aa_seqs[0].replace_range(50..56, "MKVLAW");
        aa_seqs[1] = aa_seqs[1].to_ascii_lowercase();
        aa_seqs[3].replace_range(10..20, "XXXXX*XXXX");

        let dir = std::env::temp_dir().join(format!("um-scan-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for (motifs, seqs, label) in [(&dna[..], &dna_seqs, "dna"), (&aa[..], &aa_seqs, "aa")] {
            let path = dir.join(format!("{label}.fa"));
            let records: Vec<(String, &String)> = seqs
                .iter()
                .enumerate()
                .map(|(i, s)| (format!("s{i} extra"), s))
                .collect();
            std::fs::write(&path, crate::io::fasta::write(&records, 61)).unwrap();
            for no_overlaps in [None, Some(NoOverlaps::default())] {
                let opts = ScanOptions {
                    pvalue: 0.01,
                    nthreads: 2,
                    no_overlaps,
                    ..Default::default()
                };
                let want = scan_table(motifs, seqs, &opts).unwrap();
                let matched: Vec<String> = (0..want.len()).map(|i| want.matched(i, seqs)).collect();
                assert!(want.len() > 20, "{label}: {} hits", want.len());

                let enc = EncodedSequences::new(seqs, &motifs[0].alphabet, 1).unwrap();
                let got = scan_encoded(motifs, &enc, &opts).unwrap();
                assert_eq!(got, want, "{label}: encoded");
                let got_matched: Vec<String> = (0..got.len())
                    .map(|i| got.matched_encoded(i, &enc))
                    .collect();
                assert_eq!(got_matched, matched, "{label}: encoded letters");

                // Small batches, and no room for any P-value tail.
                let text = std::fs::read(&path).unwrap();
                let gz = dir.join(format!("{label}.fa.gz"));
                let split = text.len() / 2;
                std::fs::write(&gz, crate::io::fasta::tests::gzip_members(&text, split)).unwrap();
                for (file, batch, budget) in [
                    (&path, STREAM_BATCH, TAIL_BUDGET),
                    (&path, 2_500, 0),
                    (&path, 1, 0),
                    (&gz, 2_500, TAIL_BUDGET),
                ] {
                    let f = scan_fasta_with(motifs, file, &opts, batch, budget).unwrap();
                    assert_eq!(f.table, want, "{label}: fasta, batches of {batch}");
                    let fm: Vec<&str> = (0..f.len()).map(|i| f.matched(i)).collect();
                    assert_eq!(fm, matched, "{label}: fasta letters");
                    assert_eq!(f.names[0], "s0 extra");
                    let lengths: Vec<usize> = seqs.iter().map(String::len).collect();
                    assert_eq!(f.lengths, lengths);
                }
            }
        }
        let rna = from_strings(&["ACGU"], Some(Alphabet::Rna), &copts).unwrap();
        let enc = EncodedSequences::new(&dna_seqs, &Alphabet::Dna, 1).unwrap();
        assert!(scan_encoded(&[rna], &enc, &ScanOptions::default()).is_err());
        assert!(scan_fasta(&dna, dir.join("missing.fa"), &ScanOptions::default()).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn chunked_scan_matches_plain_scoring() {
        let opts = CreateOptions {
            nsites: Some(20.0),
            pseudocount: 1.0,
            ..Default::default()
        };
        let motifs = [
            from_strings(&["TTGACGTCAAGG"], None, &opts).unwrap(),
            from_strings(&["CACGTG"], None, &opts).unwrap(),
            from_strings(&["GGGGCGGGGCCGGTTAACCA"], None, &opts).unwrap(),
        ];
        let mut seq = create_sequences(&['A', 'C', 'G', 'T'], 1, 150_000, None, 7, 1)
            .unwrap()
            .remove(0)
            .into_bytes();
        // Unknown letters, including at a chunk boundary, and lower case.
        for (a, b) in [
            (10, 13),
            (CHUNK - 3, CHUNK + 2),
            (90_000, 90_400),
            (149_998, 150_000),
        ] {
            seq[a..b].fill(b'N');
        }
        seq[20_000..21_000].make_ascii_lowercase();
        let seq = String::from_utf8(seq).unwrap();
        let opts = ScanOptions {
            pvalue: 0.01,
            nthreads: 2,
            ..Default::default()
        };
        let hits = scan_sequences(&motifs, &[&seq], &opts).unwrap();
        for (mi, m) in motifs.iter().enumerate() {
            let mut got: Vec<(usize, bool, i32)> = hits
                .iter()
                .filter(|h| h.motif == mi)
                .map(|h| {
                    (
                        h.start,
                        h.strand == Strand::Minus,
                        (h.score * 1000.0).round() as i32,
                    )
                })
                .collect();
            got.sort();
            assert_eq!(got, naive(m, &seq, 0.01), "motif {mi}");
        }
    }
}
