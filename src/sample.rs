//! Sampling sites from motifs with [`sample_sites`] and implanting them into
//! sequences with [`implant_motifs`].

use rand::RngExt;
use rand::rngs::StdRng;
use rand_distr::{Distribution, Poisson};

use crate::convert::convert_type;
use crate::error::{Error, Result, arg_err};
use crate::motif::{Motif, MotifType};
use crate::scan::reverse_complement_seq;
use crate::shuffle::{rng_for, weighted_index};

fn sample_site(ppm: &Motif, rng: &mut StdRng) -> String {
    let letters = ppm.alphabet.letters();
    ppm.matrix
        .columns()
        .map(|c| letters[weighted_index(rng, c).unwrap_or(0)])
        .collect()
}

/// Draw `n` sites from a motif, one letter per position with the motif's
/// probabilities.
pub fn sample_sites(motif: &Motif, n: usize, seed: u64) -> Result<Vec<String>> {
    let ppm = convert_type(motif, MotifType::Ppm, None, false)?;
    let mut rng = rng_for(seed, 0);
    Ok((0..n).map(|_| sample_site(&ppm, &mut rng)).collect())
}

/// How many sites [`implant_motifs`] places in each sequence.
#[derive(Debug, Clone, PartialEq)]
pub enum ImplantMode {
    /// The same number in every sequence.
    Count(usize),
    /// A Poisson number with mean `rate * sequence length`.
    Rate(f64),
    /// Given 1-based start positions for each sequence.
    Positions(Vec<Vec<usize>>),
}

/// Options for [`implant_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct ImplantOptions {
    pub mode: ImplantMode,
    /// Implant on either strand at random, rather than only on `+`.
    pub both_strands: bool,
    /// Random positions are the mean of this many uniform draws, so values
    /// above 1 favour the middle of sequences.
    pub centre_bias: usize,
    /// Fewest letters between implants in a sequence.
    pub min_spacing: usize,
    /// Attempts at finding a free position before giving up on an implant.
    pub max_retries: usize,
    pub seed: u64,
}

impl Default for ImplantOptions {
    fn default() -> Self {
        ImplantOptions {
            mode: ImplantMode::Count(1),
            both_strands: true,
            centre_bias: 1,
            min_spacing: 0,
            max_retries: 100,
            seed: 1,
        }
    }
}

/// One implanted site.
#[derive(Debug, Clone, PartialEq)]
pub struct Implant {
    pub sequence: usize,
    pub motif: usize,
    /// 1-based start.
    pub start: usize,
    pub width: usize,
    pub minus_strand: bool,
    /// The planted letters as they read on the forward strand.
    pub planted: String,
}

/// Result of [`implant_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct Implanted {
    pub sequences: Vec<String>,
    pub implants: Vec<Implant>,
    /// Problems worth reporting, such as positions that did not fit.
    pub warnings: Vec<String>,
}

fn collides(s: usize, e: usize, placed: &[(usize, usize)], gap: usize) -> bool {
    placed
        .iter()
        .any(|&(ps, pe)| s <= pe + gap && e + gap >= ps)
}

/// Plant sites sampled from DNA or RNA motifs into sequences. Each implant uses
/// a motif chosen at random when there are several, and overwrites the sequence
/// at its position.
pub fn implant_motifs<S: AsRef<str>>(
    motifs: &[Motif],
    sequences: &[S],
    opts: &ImplantOptions,
) -> Result<Implanted> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs given");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if !first.alphabet.is_nucleotide() {
        return Err(Error::Unsupported(
            "only DNA and RNA motifs can be implanted".into(),
        ));
    }
    if opts.centre_bias < 1 || opts.max_retries < 1 {
        return arg_err("centre_bias and max_retries must be at least 1");
    }
    if let ImplantMode::Positions(p) = &opts.mode
        && p.len() != sequences.len()
    {
        return arg_err("give one list of positions per sequence");
    }
    if let ImplantMode::Rate(r) = opts.mode
        && (r.is_nan() || r < 0.0)
    {
        return arg_err("rate must be non-negative");
    }
    let rna = first.alphabet == crate::alphabet::Alphabet::Rna;
    let ppms: Vec<Motif> = motifs
        .iter()
        .map(|m| convert_type(m, MotifType::Ppm, None, false))
        .collect::<Result<_>>()?;
    let widths: Vec<usize> = ppms.iter().map(|m| m.width()).collect();
    let mut rng = rng_for(opts.seed, 0);
    let mut out = Implanted {
        sequences: Vec::with_capacity(sequences.len()),
        implants: Vec::new(),
        warnings: Vec::new(),
    };
    let narrowest = *widths.iter().min().unwrap();
    if sequences
        .iter()
        .any(|s| s.as_ref().chars().count() < narrowest)
    {
        out.warnings.push(
            "some sequences are shorter than the narrowest motif; no implants will be placed in them"
                .into(),
        );
    }
    let mut short_capacity = false;

    for (i, seq) in sequences.iter().enumerate() {
        let mut chars: Vec<char> = seq.as_ref().chars().collect();
        let len = chars.len();
        let n = match &opts.mode {
            ImplantMode::Count(n) => *n,
            ImplantMode::Rate(r) => {
                let mean = r * len as f64;
                if mean > 0.0 {
                    Poisson::new(mean)
                        .map(|p| p.sample(&mut rng) as usize)
                        .unwrap_or(0)
                } else {
                    0
                }
            }
            ImplantMode::Positions(p) => p[i].len(),
        };
        let mut placed: Vec<(usize, usize)> = Vec::new();
        for j in 0..n {
            let m = if motifs.len() == 1 {
                0
            } else {
                rng.random_range(0..motifs.len())
            };
            let w = widths[m];
            if w > len {
                continue;
            }
            let minus = opts.both_strands && rng.random_range(0..2) == 1;
            let start = if let ImplantMode::Positions(p) = &opts.mode {
                let start = p[i][j];
                if start < 1 || start + w - 1 > len {
                    out.warnings.push(format!(
                        "dropping out-of-bounds position for sequence {} (start={start}, \
                         width={w}, sequence length={len})",
                        i + 1
                    ));
                    continue;
                }
                if collides(start, start + w - 1, &placed, 0) {
                    out.warnings.push(format!(
                        "position {start} in sequence {} overlaps a previous implant; planting anyway",
                        i + 1
                    ));
                }
                start
            } else {
                let mut found = None;
                for _ in 0..opts.max_retries {
                    let u: f64 = (0..opts.centre_bias)
                        .map(|_| rng.random::<f64>())
                        .sum::<f64>()
                        / opts.centre_bias as f64;
                    let start =
                        (1 + (u * (len - w + 1) as f64).floor() as usize).clamp(1, len - w + 1);
                    if !collides(start, start + w - 1, &placed, opts.min_spacing) {
                        found = Some(start);
                        break;
                    }
                }
                match found {
                    Some(s) => s,
                    None => {
                        short_capacity = true;
                        continue;
                    }
                }
            };
            let mut site = sample_site(&ppms[m], &mut rng);
            if minus {
                site = reverse_complement_seq(&site, rna);
            }
            chars.splice(start - 1..start - 1 + w, site.chars());
            placed.push((start, start + w - 1));
            out.implants.push(Implant {
                sequence: i,
                motif: m,
                start,
                width: w,
                minus_strand: minus,
                planted: site,
            });
        }
        out.sequences.push(chars.into_iter().collect());
    }
    if short_capacity {
        out.warnings.push(
            "one or more implants could not be placed within max_retries attempts; consider \
             longer sequences, fewer implants, or a smaller min_spacing"
                .into(),
        );
    }
    Ok(out)
}
