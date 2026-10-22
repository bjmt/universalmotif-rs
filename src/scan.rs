//! Scanning sequences for motif hits (`scan_sequences_lite()`).
//!
//! Each motif is turned into a PWM, a score threshold is found from the
//! P-value cutoff (see [`crate::pvalue`]), and every window of every sequence
//! scoring at least the threshold is reported, on both strands by default.
//! Scores use the same integer grid as the P-value computations, so a hit's
//! P-value is always within the cutoff.
//!
//! Sequences are packed at two bits per letter (see [`crate::packed`]), and
//! each group of four motif positions is scored with one lookup in a table of
//! all 256 four-letter combinations. Work is split into chunks of the
//! sequences, so a single long sequence still uses every thread.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::cdf::{ScoreDistribution, ScoreTail, grid_value, score_grid_ceiling};
use crate::convert::convert_type;
use crate::dedup::{Interval, dedup};
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::ops::normalize;
use crate::packed::PackedSeq;
use crate::pvalue::{scoring_matrix, with_threads};

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
}

impl Default for ScanOptions {
    fn default() -> Self {
        ScanOptions {
            pvalue: 1e-4,
            rc: true,
            no_overlaps: None,
            nthreads: 1,
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
/// stored; [`HitTable::matched`] reads them from the sequences.
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
}

impl HitTable {
    pub fn len(&self) -> usize {
        self.motif.len()
    }

    pub fn is_empty(&self) -> bool {
        self.motif.is_empty()
    }

    /// Write the letters matched by hit `i` (upper case, reverse
    /// complemented for `-` strand hits) into `out`, replacing its contents.
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
        } else {
            out.extend(bytes.iter().map(u8::to_ascii_uppercase));
        }
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

/// Window starts per unit of work; the decoded letters stay in cache.
const CHUNK: usize = 1 << 16;

/// What scanning needs to know about one motif.
struct Prepared {
    width: usize,
    min_score: i32,
    max_score: f64,
    /// One scorer per strand scanned.
    scorers: Vec<(Strand, Scorer)>,
    /// P-values of the scores at or above the threshold.
    tail: ScoreTail,
}

/// A hit before it is turned into a [`Hit`].
struct RawHit {
    motif: u32,
    sequence: u32,
    start: u32,
    minus: bool,
    score: i32,
}

/// Scan the windows starting in `starts` of one stretch of a sequence free
/// of unknown letters (ending at `seg_end`) with every motif. The packed
/// letters are decoded once into one byte per position, holding the four
/// letters that start there.
fn scan_chunk(
    prepared: &[Option<Prepared>],
    seq: &PackedSeq,
    si: usize,
    starts: (usize, usize),
    seg_end: usize,
    quads: &mut Vec<u8>,
    out: &mut Vec<RawHit>,
) {
    let (c, end) = starts;
    let max_w = prepared
        .iter()
        .flatten()
        .map(|p| p.width)
        .max()
        .unwrap_or(0);
    quads.clear();
    quads.extend((c..(end + max_w).min(seg_end + 4)).map(|i| seq.quad(i)));
    for (mi, p) in prepared.iter().enumerate() {
        let Some(p) = p else { continue };
        if seg_end < c + p.width {
            continue;
        }
        let last = end.min(seg_end + 1 - p.width);
        for (strand, scorer) in &p.scorers {
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
                    });
                }
            }
        }
    }
}

/// Scan DNA or RNA sequences for hits of the motifs, returning them as
/// columns.
///
/// Sequence letters are matched case-insensitively; other letters (such as
/// `N`) never contribute to a hit. Sequences shorter than a motif simply
/// have no hits for it. Hits are ordered by motif, sequence and start, with
/// `+` strand hits first.
pub fn scan_table<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    opts: &ScanOptions,
) -> Result<HitTable> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs to scan with");
    };
    if !(opts.pvalue > 0.0 && opts.pvalue < 1.0) {
        return arg_err("the P-value cutoff must be between 0 and 1");
    }
    let alphabet = first.alphabet.clone();
    if motifs.iter().any(|m| m.alphabet != alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if !alphabet.is_nucleotide() {
        return Err(Error::Unsupported(format!(
            "only DNA and RNA motifs can be scanned; got {alphabet}"
        )));
    }
    let rna = alphabet == Alphabet::Rna;

    // Each motif's score distribution is only needed for its threshold and
    // for the P-values above it, so just that tail is kept.
    let prepare = |m: &Motif| -> Result<Option<Prepared>> {
        let mut pwm = convert_type(m, MotifType::Pwm, None, false)?;
        if pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
            pwm = convert_type(&normalize(m)?, MotifType::Pwm, None, false)?;
        }
        let sm = scoring_matrix(&pwm, None)?;
        let dist = ScoreDistribution::new(&sm.pwm, &sm.bkg);
        let threshold = dist.score(opts.pvalue);
        let min = score_grid_ceiling(threshold);
        if min >= i32::MAX as f64 {
            return Ok(None);
        }
        let mut scorers = vec![(Strand::Plus, Scorer::new(&sm.pwm))];
        if opts.rc {
            scorers.push((Strand::Minus, Scorer::new(&sm.pwm.reversed())));
        }
        Ok(Some(Prepared {
            width: sm.pwm.ncol(),
            min_score: min as i32,
            max_score: sm
                .pwm
                .columns()
                .map(|c| c.iter().copied().fold(f64::NEG_INFINITY, f64::max))
                .sum(),
            scorers,
            tail: dist.tail(threshold),
        }))
    };
    let letters = if rna { *b"ACGU" } else { *b"ACGT" };
    let (prepared, packed) = with_threads(opts.nthreads, || {
        let prepared = motifs
            .par_iter()
            .map(prepare)
            .collect::<Result<Vec<Option<Prepared>>>>();
        let packed: Vec<PackedSeq> = sequences
            .par_iter()
            .map(|s| PackedSeq::new(s.as_ref().as_bytes(), letters))
            .collect();
        (prepared, packed)
    })?;
    let prepared = prepared?;

    // Work comes in chunks of window starts, so long sequences are shared
    // between threads too.
    let mut work: Vec<(usize, (usize, usize), usize)> = Vec::new();
    for (si, seq) in packed.iter().enumerate() {
        for (a, b) in seq.known_segments() {
            let mut c = a;
            while c < b {
                let end = (c + CHUNK).min(b);
                work.push((si, (c, end), b));
                c = end;
            }
        }
    }
    let mut raw: Vec<RawHit> = with_threads(opts.nthreads, || {
        work.par_iter()
            .fold(
                || (Vec::new(), Vec::new()),
                |(mut quads, mut out), &(si, starts, seg_end)| {
                    scan_chunk(
                        &prepared,
                        &packed[si],
                        si,
                        starts,
                        seg_end,
                        &mut quads,
                        &mut out,
                    );
                    (quads, out)
                },
            )
            .map(|(_, out)| out)
            .reduce(Vec::new, |mut a, mut b| {
                a.append(&mut b);
                a
            })
    })?;

    raw.sort_unstable_by_key(|r| (r.motif, r.sequence, r.start, r.minus));
    let pvalue_of = |r: &RawHit| {
        let p = prepared[r.motif as usize].as_ref().unwrap();
        p.tail.pvalue(r.score as f64 / 1000.0)
    };
    if let Some(no) = opts.no_overlaps {
        let intervals: Vec<Interval> = raw
            .iter()
            .map(|r| {
                let w = prepared[r.motif as usize].as_ref().unwrap().width as i64;
                Interval {
                    sequence: r.sequence as usize,
                    motif: r.motif as usize,
                    strand: r.minus as i32,
                    start: r.start as i64 + 1,
                    end: r.start as i64 + w,
                    priority: match no.by {
                        OverlapPriority::Pvalue => -pvalue_of(r),
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
        rna,
    };
    for r in &raw {
        let p = prepared[r.motif as usize].as_ref().unwrap();
        let score = r.score as f64 / 1000.0;
        t.motif.push(r.motif as i64);
        t.sequence.push(r.sequence as i64);
        t.start.push(r.start as i64 + 1);
        t.end.push((r.start as usize + p.width) as i64);
        t.minus.push(r.minus);
        t.score.push(score);
        t.score_pct.push(100.0 * score / p.max_score);
        t.pvalue.push(p.tail.pvalue(score));
    }
    Ok(t)
}

/// Scan DNA or RNA sequences for hits of the motifs, one [`Hit`] per hit
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
