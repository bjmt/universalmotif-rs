//! Conversions between motif types, information content and consensus
//! letters.

use crate::error::{Error, Result};
use crate::matrix::{sum, which_max};
use crate::motif::{Motif, MotifType};

/// Return a uniform background when `bkg` does not have one value per letter.
fn bkg_or_uniform(bkg: &[f64], n: usize) -> Vec<f64> {
    if bkg.len() == n {
        bkg.to_vec()
    } else {
        vec![1.0 / n as f64; n]
    }
}

/// Counts to probabilities. A positive `pseudocount` is spread evenly across
/// the letters.
pub fn pcm_to_ppm(pos: &[f64], pseudocount: f64) -> Vec<f64> {
    let spos = sum(pos);
    let num = pos.len() as f64;
    if pseudocount > 0.0 {
        pos.iter()
            .map(|&p| (p + pseudocount / num) / (spos + pseudocount))
            .collect()
    } else {
        pos.iter().map(|&p| p / spos).collect()
    }
}

/// Probabilities to counts. Each value is rounded, and any rounding error is
/// added to the largest count so the column sums to `nsites`. An `nsites` of 1
/// or less, or NaN, is replaced by 100.
pub fn ppm_to_pcm(pos: &[f64], nsites: f64) -> Vec<f64> {
    let nsites = if nsites.is_nan() || nsites <= 1.0 {
        100.0
    } else {
        nsites
    };
    let mut out: Vec<f64> = pos.iter().map(|&p| (p * nsites).round()).collect();
    let spos = sum(&out);
    if spos != nsites {
        let i = which_max(&out);
        out[i] += nsites - spos;
    }
    out
}

/// Probabilities to log2 odds against `bkg`. With a positive `pseudocount` the
/// column first goes through counts, using `nsites`, so the pseudocount can be
/// applied.
pub fn ppm_to_pwm(pos: &[f64], bkg: &[f64], pseudocount: f64, nsites: f64) -> Vec<f64> {
    let nsites = if nsites.is_nan() || nsites <= 1.0 {
        100.0
    } else {
        nsites
    };
    let bkg = bkg_or_uniform(bkg, pos.len());
    let pos = if pseudocount > 0.0 {
        pcm_to_ppm(&ppm_to_pcm(pos, nsites), pseudocount)
    } else {
        pos.to_vec()
    };
    pos.iter()
        .zip(&bkg)
        .map(|(&p, &b)| (p / b).log2())
        .collect()
}

/// Log2 odds back to probabilities. If `2^x` already sums to about 1 it is
/// returned as is; otherwise it is multiplied by `bkg` and, if needed,
/// renormalised.
pub fn pwm_to_ppm(pos: &[f64], bkg: &[f64]) -> Vec<f64> {
    let bkg = bkg_or_uniform(bkg, pos.len());
    let mut out: Vec<f64> = pos.iter().map(|&x| 2f64.powf(x)).collect();
    let s = sum(&out);
    if s > 0.99 && s < 1.01 {
        return out;
    }
    for (o, b) in out.iter_mut().zip(&bkg) {
        *o *= b;
    }
    let s = sum(&out);
    if s > 0.99 && s < 1.01 {
        return out;
    }
    for o in out.iter_mut() {
        *o /= s;
    }
    out
}

/// Probabilities to information content heights, either as Shannon IC, with
/// each probability scaled by the column's total IC, or as relative entropy
/// against `bkg`.
pub fn ppm_to_icm(pos: &[f64], bkg: &[f64], relative_entropy: bool) -> Vec<f64> {
    let n = pos.len();
    let bkg = if bkg.len() != n {
        vec![1.0 / n as f64; n]
    } else {
        let s = sum(bkg);
        if s > 1.01 || s < 0.99 {
            bkg.iter().map(|&b| b / s).collect()
        } else {
            bkg.to_vec()
        }
    };
    if relative_entropy {
        pos.iter()
            .zip(&bkg)
            .map(|(&p, &b)| {
                let tmp = p / b;
                let v = p * if tmp >= 0.0 { tmp.log2() } else { 0.0 };
                if v < 0.0 { 0.0 } else { v }
            })
            .collect()
    } else {
        let height_after = pos.iter().fold(0.0, |acc, &p| {
            acc + if p > 0.0 { -p * p.log2() } else { 0.0 }
        });
        let total_ic = (n as f64).log2() - height_after;
        pos.iter().map(|&p| p * total_ic).collect()
    }
}

/// Information content heights back to probabilities.
pub fn icm_to_ppm(pos: &[f64]) -> Vec<f64> {
    let total = sum(pos);
    pos.iter().map(|&p| p / total).collect()
}

/// Contribution weights to probabilities via `|x| / sum(|x|)`, the
/// TF-MoDISco convention. An all-zero column becomes uniform.
pub fn cwm_to_ppm(pos: &[f64]) -> Vec<f64> {
    let abs_sum = pos.iter().fold(0.0, |a, &v| a + v.abs());
    if abs_sum == 0.0 {
        vec![1.0 / pos.len() as f64; pos.len()]
    } else {
        pos.iter().map(|&v| v.abs() / abs_sum).collect()
    }
}

/// Total information content of a single column of the given type.
///
/// An `nsites` below 1 is replaced by 100. PCM columns get the pseudocount
/// directly; PPM columns go through counts first so the pseudocount can be
/// applied.
pub fn position_icscore(
    pos: &[f64],
    bkg: &[f64],
    motif_type: MotifType,
    pseudocount: f64,
    nsites: f64,
    relative_entropy: bool,
) -> f64 {
    let nsites = if nsites < 1.0 { 100.0 } else { nsites };
    let bkg = bkg_or_uniform(bkg, pos.len());
    let ppm = match motif_type {
        MotifType::Pcm => pcm_to_ppm(pos, pseudocount),
        MotifType::Ppm => pcm_to_ppm(&ppm_to_pcm(pos, nsites), pseudocount),
        MotifType::Pwm => pwm_to_ppm(pos, &bkg),
        MotifType::Icm => return sum(pos),
        MotifType::Cwm => cwm_to_ppm(pos),
    };
    sum(&ppm_to_icm(&ppm, &bkg, relative_entropy))
}

fn to_consensus_ppm(pos: &[f64], motif_type: MotifType, pseudocount: f64) -> Vec<f64> {
    match motif_type {
        MotifType::Pcm => pcm_to_ppm(pos, pseudocount),
        MotifType::Ppm => pos.to_vec(),
        MotifType::Pwm => pwm_to_ppm(pos, &[]),
        MotifType::Icm => icm_to_ppm(pos),
        MotifType::Cwm => cwm_to_ppm(pos),
    }
}

/// IUPAC consensus letter of a DNA or RNA column of four values in A, C, G, T/U
/// order.
pub fn consensus_letter(pos: &[f64], rna: bool, motif_type: MotifType, pseudocount: f64) -> char {
    assert_eq!(pos.len(), 4, "a DNA/RNA column needs four values");
    let p = to_consensus_ppm(pos, motif_type, pseudocount);

    if p[0] > 0.5 && p[1] <= 0.25 && p[2] <= 0.25 && p[3] <= 0.25 {
        return 'A';
    }
    if p[1] > 0.5 && p[0] <= 0.25 && p[2] <= 0.25 && p[3] <= 0.25 {
        return 'C';
    }
    if p[2] > 0.5 && p[0] <= 0.25 && p[1] <= 0.25 && p[3] <= 0.25 {
        return 'G';
    }
    if p[3] > 0.5 && p[0] <= 0.25 && p[1] <= 0.25 && p[2] <= 0.25 {
        return if rna { 'U' } else { 'T' };
    }

    if p[0] + p[1] > 0.75 {
        return 'M';
    }
    if p[0] + p[2] > 0.75 {
        return 'R';
    }
    if p[0] + p[3] > 0.75 {
        return 'W';
    }
    if p[1] + p[2] > 0.75 {
        return 'S';
    }
    if p[1] + p[3] > 0.75 {
        return 'Y';
    }
    if p[2] + p[3] > 0.75 {
        return 'K';
    }

    if p[0] > 0.25 && p[1] > 0.25 && p[3] > 0.25 {
        return 'H';
    }
    if p[1] > 0.25 && p[2] > 0.25 && p[3] > 0.25 {
        return 'B';
    }
    if p[0] > 0.25 && p[1] > 0.25 && p[2] > 0.25 {
        return 'V';
    }
    if p[0] > 0.25 && p[2] > 0.25 && p[3] > 0.25 {
        return 'D';
    }

    'N'
}

/// Consensus letter of an amino acid column of twenty values in
/// [`crate::alphabet::AMINO_ACIDS`] order, including the B, Z, J and X
/// ambiguity codes.
pub fn consensus_letter_aa(pos: &[f64], motif_type: MotifType, pseudocount: f64) -> char {
    assert_eq!(pos.len(), 20, "an amino acid column needs twenty values");
    let p = to_consensus_ppm(pos, motif_type, pseudocount);

    if p[2] >= 0.4 && p[11] >= 0.4 {
        return 'B';
    }
    if p[3] >= 0.4 && p[13] >= 0.4 {
        return 'Z';
    }
    if p[7] >= 0.4 && p[9] >= 0.4 {
        return 'J';
    }

    let imax = which_max(&p);
    if p[imax] < 0.1 {
        return 'X';
    }
    let mut sorted = p.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    if sorted[19] == sorted[18] {
        return 'X';
    }

    crate::alphabet::AMINO_ACIDS[imax]
}

/// A DNA/RNA column for an IUPAC consensus letter. `U` is treated as `T`, and
/// `N`, `-`, `+` and `.` give a uniform column.
pub fn consensus_to_ppm(letter: char) -> Result<[f64; 4]> {
    Ok(match letter {
        'A' => [0.997, 0.001, 0.001, 0.001],
        'C' => [0.001, 0.997, 0.001, 0.001],
        'G' => [0.001, 0.001, 0.997, 0.001],
        'T' | 'U' => [0.001, 0.001, 0.001, 0.997],
        'R' => [0.499, 0.001, 0.499, 0.001],
        'Y' => [0.001, 0.499, 0.001, 0.499],
        'M' => [0.499, 0.499, 0.001, 0.001],
        'K' => [0.001, 0.001, 0.499, 0.499],
        'S' => [0.001, 0.499, 0.499, 0.001],
        'W' => [0.499, 0.001, 0.001, 0.499],
        'H' => [0.333, 0.333, 0.001, 0.333],
        'B' => [0.001, 0.333, 0.333, 0.333],
        'V' => [0.333, 0.333, 0.333, 0.001],
        'D' => [0.333, 0.001, 0.333, 0.333],
        'N' | '-' | '+' | '.' => [0.25, 0.25, 0.25, 0.25],
        _ => {
            return Err(Error::Argument(format!(
                "'{letter}' is not a DNA IUPAC symbol"
            )));
        }
    })
}

/// An amino acid column for a consensus letter. B, Z and J split their mass
/// between two residues, and the other non-standard symbols X, U, O, `*`, `-`,
/// `+` and `.` give a uniform column.
pub fn consensus_to_ppm_aa(letter: char) -> Result<[f64; 20]> {
    let mut out = [0.001; 20];
    if let Some(i) = crate::alphabet::AMINO_ACIDS
        .iter()
        .position(|&l| l == letter)
    {
        out[i] = 0.981;
        return Ok(out);
    }
    match letter {
        'B' => {
            out[2] = 0.491;
            out[11] = 0.491;
        }
        'Z' => {
            out[3] = 0.491;
            out[13] = 0.491;
        }
        'J' => {
            out[7] = 0.491;
            out[9] = 0.491;
        }
        'X' | 'U' | 'O' | '*' | '-' | '+' | '.' => return Ok([0.05; 20]),
        _ => {
            return Err(Error::Argument(format!(
                "'{letter}' is not an amino acid IUPAC symbol"
            )));
        }
    }
    Ok(out)
}

/// Convert a motif to another type.
///
/// The motif's own pseudocount is used unless `pseudocount` is given, in which
/// case it also replaces the motif's pseudocount. An empty `nsites` is treated
/// as 100. Converting to [`MotifType::Cwm`] is an error, because contribution
/// weights cannot be recovered from probabilities. Converting a motif to its
/// own type returns it unchanged.
///
/// The consensus and total IC are left as they were.
///
/// # Example
///
/// ```
/// use universalmotif::MotifType;
/// use universalmotif::convert::convert_type;
/// use universalmotif::create::{CreateOptions, from_strings};
///
/// let opts = CreateOptions { nsites: Some(10.0), ..Default::default() };
/// let ppm = from_strings(&["TATAAA"], None, &opts)?;
/// let pwm = convert_type(&ppm, MotifType::Pwm, Some(1.0), false)?;
/// assert_eq!(pwm.motif_type, MotifType::Pwm);
/// assert!(pwm.matrix.get(3, 0) > 0.0); // T is favoured at the first position
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn convert_type(
    motif: &Motif,
    to: MotifType,
    pseudocount: Option<f64>,
    relative_entropy: bool,
) -> Result<Motif> {
    let from = motif.motif_type;
    if from == to {
        return Ok(motif.clone());
    }
    if to == MotifType::Cwm {
        return Err(Error::Unsupported(format!(
            "CWM cannot be derived from probabilistic motif types ({from}); a CWM must \
             originate from a contribution-attribution method (e.g. TF-MoDISco, DeepLIFT)"
        )));
    }

    let mut parts = motif.clone().into_parts();
    if let Some(p) = pseudocount {
        parts.pseudocount = p;
    }
    let pseudo = parts.pseudocount;
    let bkg = if parts.bkg.contains(&0.0) {
        let scaled: Vec<f64> = parts.bkg.iter().map(|b| b * 1000.0).collect();
        pcm_to_ppm(&scaled, 1.0)
    } else {
        parts.bkg.clone()
    };
    let nsites = parts.nsites.unwrap_or(100.0);

    let m = &parts.matrix;
    let to_icm = |ppm: &[f64]| ppm_to_icm(ppm, &bkg, relative_entropy);
    let to_pwm = |ppm: &[f64]| ppm_to_pwm(ppm, &bkg, pseudo, nsites);
    let to_pcm = |ppm: &[f64]| ppm_to_pcm(ppm, nsites);

    let mut out = match (from, to) {
        (MotifType::Pcm, MotifType::Ppm) => m.map_columns(|c| pcm_to_ppm(c, pseudo)),
        (MotifType::Pcm, MotifType::Pwm) => m.map_columns(|c| to_pwm(&pcm_to_ppm(c, pseudo))),
        (MotifType::Pcm, MotifType::Icm) => m.map_columns(|c| to_icm(&pcm_to_ppm(c, pseudo))),

        (MotifType::Ppm, MotifType::Pcm) => m.map_columns(to_pcm),
        (MotifType::Ppm, MotifType::Pwm) => m.map_columns(to_pwm),
        (MotifType::Ppm, MotifType::Icm) => m.map_columns(to_icm),

        (MotifType::Pwm, MotifType::Pcm) => m.map_columns(|c| to_pcm(&pwm_to_ppm(c, &bkg))),
        (MotifType::Pwm, MotifType::Ppm) => m.map_columns(|c| pwm_to_ppm(c, &bkg)),
        (MotifType::Pwm, MotifType::Icm) => m.map_columns(|c| to_icm(&pwm_to_ppm(c, &bkg))),

        (MotifType::Icm, MotifType::Pcm) => m.map_columns(|c| to_pcm(&icm_to_ppm(c))),
        (MotifType::Icm, MotifType::Ppm) => m.map_columns(icm_to_ppm),
        (MotifType::Icm, MotifType::Pwm) => m.map_columns(|c| to_pwm(&icm_to_ppm(c))),

        (MotifType::Cwm, MotifType::Pcm) => m.map_columns(|c| to_pcm(&cwm_to_ppm(c))),
        (MotifType::Cwm, MotifType::Ppm) => m.map_columns(cwm_to_ppm),
        (MotifType::Cwm, MotifType::Pwm) => m.map_columns(|c| to_pwm(&cwm_to_ppm(c))),
        (MotifType::Cwm, MotifType::Icm) => m.map_columns(|c| to_icm(&cwm_to_ppm(c))),

        (_, MotifType::Cwm)
        | (MotifType::Pcm, MotifType::Pcm)
        | (MotifType::Ppm, MotifType::Ppm)
        | (MotifType::Pwm, MotifType::Pwm)
        | (MotifType::Icm, MotifType::Icm) => unreachable!(),
    };

    // Any NaN becomes 0 after every conversion.
    out = out.map_columns(|c| {
        c.iter()
            .map(|&v| if v.is_nan() { 0.0 } else { v })
            .collect()
    });

    parts.matrix = out;
    parts.motif_type = to;
    Motif::from_parts(parts)
}
