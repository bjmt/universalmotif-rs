//! Scanning sequences for motif hits (`scan_sequences_lite()`).
//!
//! Each motif is turned into a PWM, a score threshold is found from the
//! P-value cutoff (see [`crate::pvalue`]), and every window of every sequence
//! scoring at least the threshold is reported, on both strands by default.
//! Scores use the same integer grid as the P-value computations, so a hit's
//! P-value is always within the cutoff.

use rayon::prelude::*;

use crate::alphabet::Alphabet;
use crate::cdf::{ScoreDistribution, grid_value, score_grid_ceiling};
use crate::convert::convert_type;
use crate::dedup::{Interval, dedup};
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType, Strand};
use crate::ops::normalize;
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

/// Integer score matrix (position-major) and its letter lookup table.
struct Scorer {
    cols: Vec<Vec<i32>>,
}

const NOT_A_LETTER: i32 = -999_999;

fn scan_one(scorer: &Scorer, seq: &[i32], min_score: i32, mut emit: impl FnMut(usize, i32)) {
    let w = scorer.cols.len();
    if seq.len() < w {
        return;
    }
    for i in 0..=seq.len() - w {
        let mut s = 0i32;
        for (j, col) in scorer.cols.iter().enumerate() {
            let l = seq[i + j];
            s += if l < 0 { NOT_A_LETTER } else { col[l as usize] };
        }
        if s >= min_score {
            emit(i, s);
        }
    }
}

/// Scan DNA or RNA sequences for hits of the motifs.
///
/// Sequence letters are matched case-insensitively; other letters (such as
/// `N`) never contribute to a hit. Sequences shorter than a motif simply
/// have no hits for it. Hits are ordered by motif, sequence and start, with
/// `+` strand hits first.
pub fn scan_sequences<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    opts: &ScanOptions,
) -> Result<Vec<Hit>> {
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

    // PWMs (reapplying the pseudocount where needed), score distributions,
    // thresholds and maximum scores.
    let mut dists = Vec::with_capacity(motifs.len());
    let mut forward = Vec::with_capacity(motifs.len());
    let mut thresholds = Vec::with_capacity(motifs.len());
    let mut max_scores = Vec::with_capacity(motifs.len());
    for m in motifs {
        let mut pwm = convert_type(m, MotifType::Pwm, None, false)?;
        if pwm.matrix.as_slice().iter().any(|v| v.is_infinite()) {
            pwm = convert_type(&normalize(m)?, MotifType::Pwm, None, false)?;
        }
        let sm = scoring_matrix(&pwm, None)?;
        let dist = ScoreDistribution::new(&sm.pwm, &sm.bkg);
        let min = score_grid_ceiling(dist.score(opts.pvalue));
        thresholds.push(if min >= i32::MAX as f64 {
            None
        } else {
            Some(min as i32)
        });
        max_scores.push(
            sm.pwm
                .columns()
                .map(|c| c.iter().copied().fold(f64::NEG_INFINITY, f64::max))
                .sum::<f64>(),
        );
        forward.push(sm.pwm);
        dists.push(dist);
    }

    let mut lookup = [-1i32; 256];
    for (i, &l) in alphabet.letters().iter().enumerate() {
        lookup[l as usize] = i as i32;
        lookup[l.to_ascii_lowercase() as usize] = i as i32;
    }
    let encoded: Vec<Vec<i32>> = sequences
        .iter()
        .map(|s| s.as_ref().bytes().map(|b| lookup[b as usize]).collect())
        .collect();

    let to_scorer = |m: &crate::matrix::Matrix| Scorer {
        cols: m
            .columns()
            .map(|c| c.iter().map(|&v| grid_value(v)).collect())
            .collect(),
    };

    let mut hits: Vec<Hit> = with_threads(opts.nthreads, || {
        (0..motifs.len())
            .into_par_iter()
            .flat_map_iter(|mi| {
                let mut out = Vec::new();
                let Some(min) = thresholds[mi] else {
                    return out;
                };
                let mut strands = vec![(Strand::Plus, to_scorer(&forward[mi]))];
                if opts.rc {
                    strands.push((Strand::Minus, to_scorer(&forward[mi].reversed())));
                }
                let w = forward[mi].ncol();
                for (si, seq) in encoded.iter().enumerate() {
                    for (strand, scorer) in &strands {
                        scan_one(scorer, seq, min, |i, s| {
                            let score = s as f64 / 1000.0;
                            let raw = &sequences[si].as_ref().as_bytes()[i..i + w];
                            let matched = if *strand == Strand::Minus {
                                rc_bytes(raw, rna)
                            } else {
                                String::from_utf8_lossy(raw).to_ascii_uppercase()
                            };
                            out.push(Hit {
                                motif: mi,
                                sequence: si,
                                start: i + 1,
                                end: i + w,
                                strand: *strand,
                                score,
                                score_pct: 100.0 * score / max_scores[mi],
                                matched,
                                pvalue: dists[mi].pvalue(score),
                            });
                        });
                    }
                }
                out
            })
            .collect()
    })?;

    hits.sort_by_key(|h| (h.motif, h.sequence, h.start, h.strand == Strand::Minus));

    if let Some(no) = opts.no_overlaps {
        let intervals: Vec<Interval> = hits
            .iter()
            .map(|h| Interval {
                sequence: h.sequence,
                motif: h.motif,
                strand: if h.strand == Strand::Minus { 1 } else { 0 },
                start: h.start as i64,
                end: h.end as i64,
                priority: match no.by {
                    OverlapPriority::Pvalue => -h.pvalue,
                    OverlapPriority::Score => h.score,
                },
            })
            .collect();
        let keep = dedup(&intervals, no.ignore_strand, no.ignore_motif);
        let mut k = keep.iter();
        hits.retain(|_| *k.next().unwrap());
    }
    Ok(hits)
}
