//! Creating motifs from matrices, consensus strings, sets of sites, or at
//! random.

use rand::rngs::StdRng;
use rand::{Rng, RngExt, SeedableRng};
use rand_distr::{Distribution, Gamma, Poisson};

use crate::alphabet::{AA_IUPAC, AMINO_ACIDS, Alphabet, DNA_IUPAC};
use crate::convert::{consensus_to_ppm, consensus_to_ppm_aa, convert_type, pcm_to_ppm};
use crate::error::{Error, Result, arg_err};
use crate::matrix::Matrix;
use crate::motif::{GapInfo, Motif, MotifBuilder, MotifType, Strand};

/// Options shared by the motif constructors in this module.
#[derive(Debug, Clone)]
pub struct CreateOptions {
    /// Type of the returned motif.
    pub motif_type: MotifType,
    /// The motif's name; `"motif"` by default.
    pub name: String,
    /// The pseudocount used in conversions; 0 by default, so that a motif
    /// keeps the frequencies it was made from. This differs from the
    /// pseudocount of 1 that [`crate::motif::Motif::builder`] starts with,
    /// which is meant for matrices from elsewhere that may hold zero counts.
    pub pseudocount: f64,
    /// The number of sites. Counts are scaled to it, and probabilities from
    /// consensus strings are rounded to it.
    pub nsites: Option<f64>,
    /// One background frequency per letter. It is uniform when `None`, except
    /// that random motifs draw one instead.
    pub bkg: Option<Vec<f64>>,
    /// An alternative name.
    pub altname: Option<String>,
    /// The transcription factor family.
    pub family: Option<String>,
    /// The organism the motif comes from.
    pub organism: Option<String>,
    /// The number of sites the background frequencies came from.
    pub bkgsites: Option<f64>,
    /// The strands the motif applies to; both by default.
    pub strand: Strand,
    /// A P-value to attach to the motif.
    pub pval: Option<f64>,
    /// A q-value to attach to the motif.
    pub qval: Option<f64>,
    /// An E-value to attach to the motif.
    pub eval: Option<f64>,
    /// Free-form `(name, value)` pairs.
    pub extrainfo: Vec<(String, String)>,
}

impl Default for CreateOptions {
    fn default() -> Self {
        CreateOptions {
            motif_type: MotifType::Ppm,
            name: "motif".into(),
            pseudocount: 0.0,
            nsites: None,
            bkg: None,
            altname: None,
            family: None,
            organism: None,
            bkgsites: None,
            strand: Strand::Both,
            pval: None,
            qval: None,
            eval: None,
            extrainfo: Vec::new(),
        }
    }
}

impl CreateOptions {
    fn builder(&self, matrix: Matrix, alphabet: Alphabet) -> MotifBuilder {
        Motif::builder(matrix, alphabet)
            .name(self.name.clone())
            .pseudocount(self.pseudocount)
            .nsites(self.nsites)
            .bkg(self.bkg.clone())
            .altname(self.altname.clone())
            .family(self.family.clone())
            .organism(self.organism.clone())
            .bkgsites(self.bkgsites)
            .strand(self.strand)
            .pval(self.pval)
            .qval(self.qval)
            .eval(self.eval)
            .extrainfo(self.extrainfo.clone())
            .gap(GapInfo::default())
    }

    fn check(&self) -> Result<()> {
        if let Some(n) = self.nsites
            && (!n.is_finite() || n < 1.0 || n != n.floor())
        {
            return arg_err("nsites must be a positive integer");
        }
        Ok(())
    }
}

/// Guess the alphabet of a matrix from its row letters: DNA or RNA for four
/// nucleotide letters, AA for the twenty amino acids, and otherwise a custom
/// alphabet of the letters.
pub fn guess_matrix_alphabet(letters: &[char]) -> Result<Alphabet> {
    let all_in = |set: &[char]| letters.iter().all(|l| set.contains(l));
    if letters.len() == 4 && all_in(&crate::alphabet::DNA) {
        Ok(Alphabet::Dna)
    } else if letters.len() == 4 && all_in(&crate::alphabet::RNA) {
        Ok(Alphabet::Rna)
    } else if letters.len() == 20 && all_in(&AMINO_ACIDS) {
        Ok(Alphabet::Aa)
    } else {
        Alphabet::custom(letters.iter().copied())
    }
}

/// Create a motif from a matrix.
///
/// If `row_letters` is given, the rows are first sorted by letter, and the
/// alphabet is guessed from the letters when `alphabet` is `None`. Without row
/// letters the alphabet must be given and the rows must already be in its
/// order. The matrix type is detected from the values unless `opts.motif_type`
/// is CWM, and the motif is then converted to `opts.motif_type`.
pub fn from_matrix(
    matrix: Matrix,
    row_letters: Option<&[char]>,
    alphabet: Option<Alphabet>,
    opts: &CreateOptions,
) -> Result<Motif> {
    let detected = (opts.motif_type == MotifType::Cwm).then_some(MotifType::Cwm);
    matrix_motif(matrix, row_letters, alphabet, detected, opts)
}

/// A motif from a matrix whose type is known rather than detected, such as
/// the weights of a convolution filter, which can be log-odds scores with no
/// negative values. It is returned as `opts.motif_type`, as for
/// [`from_matrix`].
pub fn from_typed_matrix(
    matrix: Matrix,
    row_letters: Option<&[char]>,
    alphabet: Option<Alphabet>,
    matrix_type: MotifType,
    opts: &CreateOptions,
) -> Result<Motif> {
    matrix_motif(matrix, row_letters, alphabet, Some(matrix_type), opts)
}

fn matrix_motif(
    matrix: Matrix,
    row_letters: Option<&[char]>,
    alphabet: Option<Alphabet>,
    detected: Option<MotifType>,
    opts: &CreateOptions,
) -> Result<Motif> {
    opts.check()?;
    if matrix.as_slice().iter().any(|v| v.is_nan()) {
        return arg_err("matrix cannot have NaN values");
    }

    let (matrix, letters) = match row_letters {
        Some(l) => {
            if l.len() != matrix.nrow() {
                return arg_err("number of row letters does not match matrix rows");
            }
            let mut order: Vec<usize> = (0..l.len()).collect();
            order.sort_by_key(|&i| l[i]);
            let sorted: Vec<char> = order.iter().map(|&i| l[i]).collect();
            if sorted.windows(2).any(|w| w[0] == w[1]) {
                return arg_err("row letters must be unique");
            }
            let mut data = Vec::with_capacity(matrix.as_slice().len());
            for c in matrix.columns() {
                data.extend(order.iter().map(|&i| c[i]));
            }
            (
                Matrix::from_vec(matrix.nrow(), matrix.ncol(), data)?,
                Some(sorted),
            )
        }
        None => (matrix, None),
    };

    let alphabet = match (alphabet, &letters) {
        (Some(a @ Alphabet::Custom(_)), Some(l)) => {
            if l.iter().any(|c| a.index_of(*c).is_none()) {
                return arg_err("row letters do not match the provided alphabet");
            }
            a
        }
        (Some(a), _) => a,
        (None, Some(l)) => guess_matrix_alphabet(l)?,
        (None, None) => {
            return arg_err("please provide an alphabet or a matrix with row letters");
        }
    };
    match alphabet {
        Alphabet::Dna | Alphabet::Rna if matrix.nrow() != 4 => {
            return arg_err("incorrect number of rows");
        }
        Alphabet::Aa if matrix.nrow() != 20 => return arg_err("incorrect number of rows"),
        Alphabet::Custom(ref l) if l.len() != matrix.nrow() => {
            return arg_err("alphabet length does not match number of rows");
        }
        _ => {}
    }

    let motif = opts
        .builder(matrix, alphabet)
        .motif_type(detected)
        .build()?;
    convert_type(&motif, opts.motif_type, None, false)
}

fn sorted_unique(seqs: &[&str]) -> Vec<char> {
    let mut l: Vec<char> = seqs.iter().flat_map(|s| s.chars()).collect();
    l.sort_unstable();
    l.dedup();
    l
}

/// Guess the alphabet of a set of strings the way [`from_strings`] does.
pub fn guess_strings_alphabet(seqs: &[&str]) -> Result<Alphabet> {
    guess_string_alphabet(&sorted_unique(seqs))
}

fn guess_string_alphabet(letters: &[char]) -> Result<Alphabet> {
    let all_in = |set: &str| letters.iter().all(|&l| set.contains(l));
    if all_in(DNA_IUPAC) {
        Ok(Alphabet::Dna)
    } else if all_in("ACGU") {
        Ok(Alphabet::Rna)
    } else if all_in(AA_IUPAC) {
        Ok(Alphabet::Aa)
    } else {
        Alphabet::custom(letters.iter().copied())
    }
}

fn check_string_letters(letters: &[char], alphabet: &Alphabet) -> Result<()> {
    let (allowed, what): (String, &str) = match alphabet {
        Alphabet::Dna => (DNA_IUPAC.into(), "DNA "),
        Alphabet::Rna => ("ACGU".into(), "RNA "),
        Alphabet::Aa => (AA_IUPAC.into(), "AA "),
        Alphabet::Custom(l) => (l.iter().collect(), ""),
    };
    let bad: String = letters.iter().filter(|&&l| !allowed.contains(l)).collect();
    if bad.is_empty() {
        Ok(())
    } else {
        arg_err(format!(
            "the input has letters outside the {what}alphabet: {bad}"
        ))
    }
}

/// Create a motif from strings.
///
/// A single string is read as a consensus sequence, in which IUPAC codes are
/// accepted for DNA and amino acids. Several strings of equal length are read
/// as aligned sites and counted; `nsites` then defaults to the number of
/// strings. The alphabet is guessed from the letters when `None`: DNA with
/// IUPAC codes, then RNA, then amino acids, and otherwise a custom alphabet of
/// the letters present.
///
/// # Example
///
/// ```
/// use universalmotif::MotifType;
/// use universalmotif::create::{CreateOptions, from_strings};
///
/// let opts = CreateOptions { motif_type: MotifType::Pcm, ..Default::default() };
/// let m = from_strings(&["TGACTCA", "TGAGTCA"], None, &opts)?;
/// assert_eq!(m.width(), 7);
/// assert_eq!(m.matrix.get(3, 0), 2.0); // both sites start with T
/// # Ok::<(), universalmotif::Error>(())
/// ```
pub fn from_strings(
    seqs: &[&str],
    alphabet: Option<Alphabet>,
    opts: &CreateOptions,
) -> Result<Motif> {
    opts.check()?;
    let Some(first) = seqs.first() else {
        return arg_err("at least one string is required");
    };
    let width = first.chars().count();
    if seqs.iter().any(|s| s.chars().count() != width) {
        return arg_err("all strings must be of equal length");
    }
    if width == 0 {
        return arg_err("zero-length input strings are not allowed");
    }

    let letters = sorted_unique(seqs);
    let alphabet = match alphabet {
        None => guess_string_alphabet(&letters)?,
        Some(a) => {
            check_string_letters(&letters, &a)?;
            a
        }
    };

    let (matrix, opts) = if seqs.len() == 1 {
        let s = *first;
        let cols: Vec<Vec<f64>> = match &alphabet {
            Alphabet::Dna | Alphabet::Rna => s
                .chars()
                .map(|c| consensus_to_ppm(c).map(|p| p.to_vec()))
                .collect::<Result<_>>()?,
            Alphabet::Aa => s
                .chars()
                .map(|c| consensus_to_ppm_aa(c).map(|p| p.to_vec()))
                .collect::<Result<_>>()?,
            Alphabet::Custom(_) => s
                .chars()
                .map(|c| {
                    let mut col = vec![0.0; alphabet.len()];
                    col[alphabet.index_of(c).expect("letters checked above")] = 1.0;
                    col
                })
                .collect(),
        };
        (Matrix::from_columns(&cols)?, opts.clone())
    } else {
        match &alphabet {
            Alphabet::Dna | Alphabet::Rna => {
                if letters.iter().any(|&l| !"ACGTU".contains(l)) {
                    return arg_err(format!(
                        "{alphabet} ambiguity letters are not accepted if input is more than one \
                         string; to count them as letters of their own, give a custom alphabet"
                    ));
                }
            }
            Alphabet::Aa => {
                if letters.iter().any(|&l| "UOBJZX".contains(l)) {
                    return arg_err(
                        "AA ambiguity letters are not accepted if input is more than one string; \
                         to count them as letters of their own, give a custom alphabet",
                    );
                }
            }
            Alphabet::Custom(_) => {}
        }
        let n = alphabet.len();
        let mut counts = vec![vec![0.0; n]; width];
        for s in seqs {
            for (j, c) in s.chars().enumerate() {
                if let Some(i) = alphabet.index_of(c) {
                    counts[j][i] += 1.0;
                }
            }
        }
        let cols: Vec<Vec<f64>> = counts.iter().map(|c| pcm_to_ppm(c, 0.0)).collect();
        let mut opts = opts.clone();
        opts.nsites = opts.nsites.or(Some(seqs.len() as f64));
        (Matrix::from_columns(&cols)?, opts)
    };

    let motif = opts
        .builder(matrix, alphabet)
        .motif_type(Some(MotifType::Ppm))
        .build()?;
    // PPMs are passed through counts before being returned, which cleans up
    // the small values used for consensus letters.
    let motif = if opts.motif_type == MotifType::Ppm {
        convert_type(&motif, MotifType::Pcm, None, false)?
    } else {
        motif
    };
    convert_type(&motif, opts.motif_type, None, false)
}

/// Create a random motif of the given width.
///
/// Each column is drawn from a Dirichlet distribution whose parameters are the
/// background frequencies. Without a background, one is drawn too. If
/// `nsites` is not given, a PPM gets a random value between 50 and 250, and a
/// PCM gets 100. The same `seed` always gives the same motif.
pub fn random(width: usize, alphabet: Alphabet, opts: &CreateOptions, seed: u64) -> Result<Motif> {
    random_with(width, alphabet, opts, &mut StdRng::seed_from_u64(seed))
}

fn random_with<R: Rng + ?Sized>(
    width: usize,
    alphabet: Alphabet,
    opts: &CreateOptions,
    rng: &mut R,
) -> Result<Motif> {
    opts.check()?;
    if width == 0 {
        return arg_err("width must be greater than 0");
    }
    let n = alphabet.len();
    let bkg = match &opts.bkg {
        Some(b) if b.len() != n => {
            return arg_err(format!("bkg must have {n} values, one per letter"));
        }
        Some(b) => b.clone(),
        None => {
            let pois =
                Poisson::new(1000.0 / n as f64).map_err(|e| Error::Argument(e.to_string()))?;
            let draws: Vec<f64> = (0..n).map(|_| pois.sample(rng) / 1000.0).collect();
            let s = crate::matrix::sum(&draws);
            draws.iter().map(|d| d / s).collect()
        }
    };

    let mut data = Vec::with_capacity(width * n);
    for _ in 0..width {
        let col: Vec<f64> = bkg
            .iter()
            .map(|&a| {
                if a > 0.0 {
                    Gamma::new(a, 1.0).map(|g| g.sample(rng)).unwrap_or(0.0)
                } else {
                    0.0
                }
            })
            .collect();
        let s = crate::matrix::sum(&col);
        data.extend(col.iter().map(|v| v / s));
    }
    let matrix = Matrix::from_vec(n, width, data)?;

    let mut opts = opts.clone();
    opts.bkg = Some(bkg);
    if opts.nsites.is_none() {
        opts.nsites = match opts.motif_type {
            MotifType::Ppm => Some(rng.random_range(50..=250) as f64),
            MotifType::Pcm => Some(100.0),
            _ => None,
        };
    }
    from_matrix(matrix, None, Some(alphabet), &opts)
}
