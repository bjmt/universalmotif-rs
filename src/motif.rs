//! The motif type, its builder and its validation rules.

use std::fmt;
use std::ops::Deref;

use crate::alphabet::Alphabet;
use crate::convert::{consensus_letter, consensus_letter_aa, position_icscore};
use crate::error::{Error, Result, arg_err};
use crate::matrix::{Matrix, sum, which_max};

/// What the values of a motif matrix represent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MotifType {
    /// Position count matrix.
    Pcm,
    /// Position probability matrix.
    Ppm,
    /// Position weight matrix of log2 odds.
    Pwm,
    /// Information content matrix.
    Icm,
    /// Contribution weight matrix, such as those from TF-MoDISco.
    Cwm,
}

impl MotifType {
    /// Every motif type, in the order PCM, PPM, PWM, ICM, CWM.
    pub const ALL: [MotifType; 5] = [
        MotifType::Pcm,
        MotifType::Ppm,
        MotifType::Pwm,
        MotifType::Icm,
        MotifType::Cwm,
    ];

    /// Parse a type name: `"PCM"`, `"PPM"`, `"PWM"`, `"ICM"` or `"CWM"`.
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "PCM" => Ok(MotifType::Pcm),
            "PPM" => Ok(MotifType::Ppm),
            "PWM" => Ok(MotifType::Pwm),
            "ICM" => Ok(MotifType::Icm),
            "CWM" => Ok(MotifType::Cwm),
            _ => arg_err(format!(
                "incorrect type '{s}' (must be one of PCM, PPM, PWM, ICM, CWM)"
            )),
        }
    }

    /// The type's name in upper case, such as `"PPM"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            MotifType::Pcm => "PCM",
            MotifType::Ppm => "PPM",
            MotifType::Pwm => "PWM",
            MotifType::Icm => "ICM",
            MotifType::Cwm => "CWM",
        }
    }
}

impl fmt::Display for MotifType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Which strand(s) a motif applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Strand {
    /// The forward strand only.
    Plus,
    /// The reverse strand only.
    Minus,
    /// Both strands.
    #[default]
    Both,
}

impl Strand {
    /// Parse `"+"`, `"-"` or `"+-"`; `"-+"` and `"*"` also mean both strands.
    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            "+-" | "-+" | "*" => Ok(Strand::Both),
            _ => arg_err(format!("incorrect strand '{s}' (must be one of +, -, +-)")),
        }
    }

    /// `"+"`, `"-"` or `"+-"`.
    pub fn as_str(&self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Both => "+-",
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Gaps inside a motif. Gap `i` sits after column `gaploc[i]`, counting from 1,
/// and spans `mingap[i]` to `maxgap[i]` positions.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GapInfo {
    /// Whether the motif has any gaps.
    pub isgapped: bool,
    /// The column after which each gap sits, counting from 1.
    pub gaploc: Vec<usize>,
    /// The fewest positions each gap can span.
    pub mingap: Vec<usize>,
    /// The most positions each gap can span.
    pub maxgap: Vec<usize>,
}

/// Every field of a motif, with no invariants enforced.
///
/// Use [`Motif::from_parts`] to turn this into a validated [`Motif`]. Unlike
/// [`MotifBuilder`], it applies no defaults and recomputes nothing, which
/// makes it the right tool for restoring a motif exactly as it was saved.
#[derive(Debug, Clone)]
pub struct MotifParts {
    /// The motif's name.
    pub name: String,
    /// An alternative name, such as the transcription factor's name when the
    /// main name is a database ID.
    pub altname: Option<String>,
    /// The transcription factor family.
    pub family: Option<String>,
    /// The organism the motif comes from.
    pub organism: Option<String>,
    /// The motif matrix, with one row per alphabet letter and one column per
    /// position.
    pub matrix: Matrix,
    /// The alphabet the rows of the matrix follow.
    pub alphabet: Alphabet,
    /// What the values of the matrix represent.
    pub motif_type: MotifType,
    /// The total information content of the motif, in bits.
    pub icscore: f64,
    /// The number of sites the motif was built from.
    pub nsites: Option<f64>,
    /// The pseudocount added when converting to probabilities or weights.
    pub pseudocount: f64,
    /// Background frequencies, one per alphabet letter in row order.
    pub bkg: Vec<f64>,
    /// The number of sites the background frequencies came from.
    pub bkgsites: Option<f64>,
    /// One consensus letter per column, or `None` for custom alphabets.
    pub consensus: Option<String>,
    /// The strands the motif applies to.
    pub strand: Strand,
    /// A P-value attached to the motif, such as from discovery.
    pub pval: Option<f64>,
    /// A q-value attached to the motif.
    pub qval: Option<f64>,
    /// An E-value attached to the motif.
    pub eval: Option<f64>,
    /// Free-form `(name, value)` pairs. Names may be empty.
    pub extrainfo: Vec<(String, String)>,
    /// The motif's gaps, if it has any.
    pub gap: GapInfo,
}

// Equality treats two NaN totals of IC, as CWMs have, as equal.
impl PartialEq for MotifParts {
    fn eq(&self, o: &Self) -> bool {
        let MotifParts {
            name,
            altname,
            family,
            organism,
            matrix,
            alphabet,
            motif_type,
            icscore,
            nsites,
            pseudocount,
            bkg,
            bkgsites,
            consensus,
            strand,
            pval,
            qval,
            eval,
            extrainfo,
            gap,
        } = self;
        let same_ic = icscore == &o.icscore || (icscore.is_nan() && o.icscore.is_nan());
        same_ic
            && name == &o.name
            && altname == &o.altname
            && family == &o.family
            && organism == &o.organism
            && matrix == &o.matrix
            && alphabet == &o.alphabet
            && motif_type == &o.motif_type
            && nsites == &o.nsites
            && pseudocount == &o.pseudocount
            && bkg == &o.bkg
            && bkgsites == &o.bkgsites
            && consensus == &o.consensus
            && strand == &o.strand
            && pval == &o.pval
            && qval == &o.qval
            && eval == &o.eval
            && extrainfo == &o.extrainfo
            && gap == &o.gap
    }
}

impl MotifParts {
    /// Check every validity rule and return a description of each failure.
    pub fn validate(&self) -> Vec<String> {
        let mut msg = Vec::new();
        let m = &self.matrix;
        let ncol = m.ncol();

        let g = &self.gap;
        if g.gaploc.len() != g.mingap.len() || g.gaploc.len() != g.maxgap.len() {
            msg.push("gaploc, mingap and maxgap should all be the same length".into());
        }
        for &loc in &g.gaploc {
            if loc == 0 {
                msg.push("position 0 gaps or less are not allowed".into());
            }
            if loc >= ncol {
                msg.push("gap location values should not exceed motif size - 1".into());
            }
        }
        for (&lo, &hi) in g.mingap.iter().zip(&g.maxgap) {
            if lo > hi {
                msg.push("mingap cannot be greater than the corresponding maxgap".into());
            }
            if hi == 0 {
                msg.push("maxgap must be greater than 0".into());
            }
        }

        if ncol == 0 {
            msg.push("motif must have at least one column".into());
        }
        if m.as_slice().iter().any(|v| v.is_nan()) {
            msg.push("motif cannot contain NaN values".into());
        }
        if !self.pseudocount.is_finite() || self.pseudocount < 0.0 {
            msg.push("pseudocount must be a finite, non-negative number".into());
        }

        let nlet = self.alphabet.len();
        match self.alphabet {
            Alphabet::Dna | Alphabet::Rna if m.nrow() != 4 => {
                msg.push("DNA/RNA motifs must have 4 rows".into())
            }
            Alphabet::Aa if m.nrow() != 20 => msg.push("AA motifs must have 20 rows".into()),
            Alphabet::Custom(_) if m.nrow() != nlet => {
                msg.push("alphabet length does not match number of rows in motif".into())
            }
            _ => {}
        }

        if self.bkg.len() != m.nrow() {
            msg.push("bkg must have one value per alphabet letter".into());
        } else {
            if self.bkg.iter().any(|&b| b < 0.0) {
                msg.push("bkg does not allow values less than 0".into());
            }
            if self.bkg.iter().any(|&b| b > 1.0) {
                msg.push("bkg does not allow values higher than 1".into());
            }
            let s = sum(&self.bkg);
            if !(0.99..=1.01).contains(&s) {
                msg.push("bkg probabilities must add up to 1".into());
            }
        }

        let colsums = m.col_sums();
        match self.motif_type {
            MotifType::Pcm => {
                if let Some(&first) = colsums.first() {
                    let tol = f64::max(0.01, first.abs() * 0.01);
                    if colsums.iter().skip(1).any(|&s| (s - first).abs() > tol) {
                        msg.push("PCM columns must all have the same sum".into());
                    }
                }
                if m.as_slice().iter().any(|&v| v > 0.0 && v < 1.0) {
                    msg.push("PCM values must be whole counts, not between 0 and 1".into());
                }
            }
            MotifType::Ppm => {
                if !colsums.iter().all(|&s| s > 0.99 && s < 1.01) {
                    msg.push("PPM columns must sum to 1".into());
                }
                if !m.as_slice().iter().all(|&v| v >= 0.0) {
                    msg.push("PPM values cannot be negative".into());
                }
            }
            MotifType::Icm => {
                if m.as_slice().iter().any(|&v| v < 0.0) {
                    msg.push("ICM values cannot be negative".into());
                }
            }
            MotifType::Pwm | MotifType::Cwm => {}
        }

        if let Some(c) = &self.consensus
            && c.chars().count() != ncol
        {
            msg.push(
                "consensus string must have the same number of letters as motif positions".into(),
            );
        }

        msg
    }

    /// Recompute the total IC and the consensus from the matrix, type,
    /// background, pseudocount and `nsites`, for use after changing any of
    /// them.
    pub fn recompute_derived(&mut self) {
        self.icscore = Motif::compute_icscore(
            &self.matrix,
            &self.bkg,
            self.motif_type,
            self.pseudocount,
            self.nsites,
        );
        self.consensus = Motif::compute_consensus(
            &self.matrix,
            &self.alphabet,
            self.motif_type,
            self.pseudocount,
        );
    }
}

/// A validated sequence motif.
///
/// Fields are read through [`Deref`] to [`MotifParts`], so `motif.name`,
/// `motif.matrix` and so on work directly. To change a motif, go through
/// [`Motif::into_parts`] and [`Motif::from_parts`], or one of the functions in
/// [`crate::convert`] and [`crate::ops`].
#[derive(Debug, Clone, PartialEq)]
pub struct Motif(MotifParts);

impl Deref for Motif {
    type Target = MotifParts;
    fn deref(&self) -> &MotifParts {
        &self.0
    }
}

impl Motif {
    /// Validate a set of parts, unchanged, as a motif.
    pub fn from_parts(parts: MotifParts) -> Result<Motif> {
        let msg = parts.validate();
        if msg.is_empty() {
            Ok(Motif(parts))
        } else {
            Err(Error::InvalidMotif(msg))
        }
    }

    /// Take the motif apart into its fields, to change them freely and validate
    /// the result again with [`Motif::from_parts`].
    pub fn into_parts(self) -> MotifParts {
        self.0
    }

    /// Start building a motif from a matrix whose rows follow the order of
    /// `alphabet`'s letters.
    pub fn builder(matrix: Matrix, alphabet: Alphabet) -> MotifBuilder {
        MotifBuilder::new(matrix, alphabet)
    }

    /// Motif width, the number of positions.
    pub fn width(&self) -> usize {
        self.0.matrix.ncol()
    }

    /// Per-column consensus letters, recomputed from the current matrix.
    /// Returns `None` for custom alphabets.
    pub(crate) fn compute_consensus(
        matrix: &Matrix,
        alphabet: &Alphabet,
        motif_type: MotifType,
        pseudocount: f64,
    ) -> Option<String> {
        match alphabet {
            Alphabet::Dna | Alphabet::Rna => {
                let rna = *alphabet == Alphabet::Rna;
                Some(
                    matrix
                        .columns()
                        .map(|c| consensus_letter(c, rna, motif_type, pseudocount))
                        .collect(),
                )
            }
            Alphabet::Aa => Some(
                matrix
                    .columns()
                    .map(|c| consensus_letter_aa(c, motif_type, pseudocount))
                    .collect(),
            ),
            Alphabet::Custom(_) => None,
        }
    }

    /// Total information content of a matrix, summed over columns. An empty
    /// `nsites` is passed on as 0, which [`position_icscore`] treats as 100.
    /// Contribution weights carry no information content, so a CWM gets NaN.
    pub(crate) fn compute_icscore(
        matrix: &Matrix,
        bkg: &[f64],
        motif_type: MotifType,
        pseudocount: f64,
        nsites: Option<f64>,
    ) -> f64 {
        if motif_type == MotifType::Cwm {
            return f64::NAN;
        }
        let nsites = nsites.unwrap_or(0.0);
        let per_col: Vec<f64> = matrix
            .columns()
            .map(|c| position_icscore(c, bkg, motif_type, pseudocount, nsites, false))
            .collect();
        sum(&per_col)
    }
}

/// Builds a [`Motif`], filling in defaults for anything not given.
///
/// - The type is guessed when not given: PCM if every value is 0 or at least
///   1, PPM if columns sum to within 0.02 of 1, ICM if nothing is negative,
///   and PWM otherwise. A PPM whose columns do not sum to 1 is renormalized.
/// - For a PCM, `nsites` defaults to the first column sum; if a different
///   `nsites` is given, every column is rescaled to it and rounded.
/// - The background defaults to uniform.
/// - The total IC and the consensus are computed from the matrix.
/// - An empty name becomes `"[no name]"`, and empty metadata strings or NaN
///   numbers are dropped.
#[derive(Debug, Clone)]
pub struct MotifBuilder {
    matrix: Matrix,
    alphabet: Alphabet,
    motif_type: Option<MotifType>,
    name: String,
    altname: Option<String>,
    family: Option<String>,
    organism: Option<String>,
    icscore: Option<f64>,
    nsites: Option<f64>,
    pseudocount: f64,
    bkg: Option<Vec<f64>>,
    bkgsites: Option<f64>,
    strand: Strand,
    pval: Option<f64>,
    qval: Option<f64>,
    eval: Option<f64>,
    extrainfo: Vec<(String, String)>,
    gap: GapInfo,
}

macro_rules! setter {
    ($(#[$m:meta])* $name:ident: $ty:ty) => {
        $(#[$m])*
        pub fn $name(mut self, value: $ty) -> Self {
            self.$name = value;
            self
        }
    };
}

impl MotifBuilder {
    /// Start building a motif from its matrix and alphabet, with the name `"new
    /// motif"`, a pseudocount of 1, both strands and nothing else set.
    pub fn new(matrix: Matrix, alphabet: Alphabet) -> Self {
        MotifBuilder {
            matrix,
            alphabet,
            motif_type: None,
            name: "new motif".into(),
            altname: None,
            family: None,
            organism: None,
            icscore: None,
            nsites: None,
            pseudocount: 1.0,
            bkg: None,
            bkgsites: None,
            strand: Strand::Both,
            pval: None,
            qval: None,
            eval: None,
            extrainfo: Vec::new(),
            gap: GapInfo::default(),
        }
    }

    /// The motif's name.
    pub fn name(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    setter!(
        /// The matrix type; guessed from the values when `None`.
        motif_type: Option<MotifType>
    );
    setter!(
        /// An alternative name, such as the transcription factor's name.
        altname: Option<String>
    );
    setter!(
        /// The transcription factor family.
        family: Option<String>
    );
    setter!(
        /// The organism the motif comes from.
        organism: Option<String>
    );
    setter!(
        /// The total IC; computed from the matrix when `None`.
        icscore: Option<f64>
    );
    setter!(
        /// The number of sites the motif was built from.
        nsites: Option<f64>
    );
    setter!(
        /// The pseudocount used in conversions; 1 by default.
        pseudocount: f64
    );
    setter!(
        /// Background frequencies, one per letter; uniform when `None`.
        bkg: Option<Vec<f64>>
    );
    setter!(
        /// The number of sites the background frequencies came from.
        bkgsites: Option<f64>
    );
    setter!(
        /// The strands the motif applies to; both by default.
        strand: Strand
    );
    setter!(
        /// A P-value attached to the motif.
        pval: Option<f64>
    );
    setter!(
        /// A q-value attached to the motif.
        qval: Option<f64>
    );
    setter!(
        /// An E-value attached to the motif.
        eval: Option<f64>
    );
    setter!(
        /// Free-form `(name, value)` pairs.
        extrainfo: Vec<(String, String)>
    );
    setter!(
        /// The motif's gaps.
        gap: GapInfo
    );

    /// Check the matrix against the alphabet and build the motif.
    ///
    /// Without a type, it is guessed from the values: counts if every value
    /// is 0 or at least 1, probabilities if the columns sum to within about
    /// 0.02 of 1, information content if no value is negative, and weights
    /// otherwise. Probabilities that do not sum to 1 are renormalized. A PCM
    /// is rescaled to `nsites` if that is given, and gives `nsites` its column
    /// sum if not. The total IC is computed unless given, and the consensus
    /// always is.
    pub fn build(self) -> Result<Motif> {
        let mut m = self.matrix;
        let nrow = m.nrow();

        if m.ncol() == 0 {
            return arg_err("motif must have at least one column");
        }
        if m.as_slice().iter().any(|v| v.is_nan()) {
            return arg_err("matrix cannot have NaN values");
        }
        match &self.alphabet {
            Alphabet::Dna | Alphabet::Rna if nrow != 4 => {
                return arg_err("DNA/RNA motifs must have 4 rows");
            }
            Alphabet::Aa if nrow != 20 => return arg_err("AA motifs must have 20 rows"),
            Alphabet::Custom(l) if l.len() != nrow => {
                return arg_err("alphabet length does not match matrix rows");
            }
            _ => {}
        }

        let colsums = m.col_sums();
        let vals = m.as_slice();
        let all_pos = vals.iter().all(|&v| v >= 0.0);
        let sums_one = colsums.iter().all(|&s| s > 0.99 && s < 1.01);

        let motif_type = match self.motif_type {
            Some(t) => {
                if t == MotifType::Ppm && !sums_one && all_pos {
                    for (j, &s) in colsums.iter().enumerate() {
                        m.col_mut(j).iter_mut().for_each(|v| *v /= s);
                    }
                }
                t
            }
            None => {
                // Probabilities written to two decimals can sum to anywhere
                // from about 0.98 to 1.02, so guessing allows for that.
                let near_one = colsums.iter().all(|&s| (s - 1.0).abs() <= 0.0201);
                if vals.iter().all(|&v| v >= 1.0 || v == 0.0) {
                    MotifType::Pcm
                } else if near_one && all_pos {
                    if !sums_one {
                        for (j, &s) in colsums.iter().enumerate() {
                            m.col_mut(j).iter_mut().for_each(|v| *v /= s);
                        }
                    }
                    MotifType::Ppm
                } else if all_pos {
                    MotifType::Icm
                } else {
                    MotifType::Pwm
                }
            }
        };

        let nsites = match self.nsites.filter(|n| !n.is_nan()) {
            None if motif_type == MotifType::Pcm => Some(sum(m.col(0))),
            None => None,
            Some(n) => {
                if motif_type == MotifType::Pcm && colsums.iter().any(|&s| s != n) {
                    for (j, &s) in colsums.iter().enumerate() {
                        let col = m.col_mut(j);
                        for v in col.iter_mut() {
                            *v = (*v / s * n).round();
                        }
                        let possum = sum(col);
                        if possum != n {
                            let i = which_max(col);
                            col[i] += n - possum;
                        }
                    }
                }
                Some(n)
            }
        };

        let bkg = match self.bkg {
            None => vec![1.0 / nrow as f64; nrow],
            Some(b) if b.len() != nrow => {
                return arg_err(format!(
                    "bkg has {} values but the alphabet has {} letters",
                    b.len(),
                    nrow
                ));
            }
            Some(b) => b,
        };

        let icscore = self
            .icscore
            .filter(|v| !v.is_nan() && motif_type != MotifType::Cwm)
            .unwrap_or_else(|| {
                Motif::compute_icscore(&m, &bkg, motif_type, self.pseudocount, nsites)
            });
        let consensus = Motif::compute_consensus(&m, &self.alphabet, motif_type, self.pseudocount);

        let nonempty = |s: Option<String>| s.filter(|s| !s.is_empty());
        let notnan = |x: Option<f64>| x.filter(|x| !x.is_nan());

        Motif::from_parts(MotifParts {
            name: if self.name.is_empty() {
                "[no name]".into()
            } else {
                self.name
            },
            altname: nonempty(self.altname),
            family: nonempty(self.family),
            organism: nonempty(self.organism),
            matrix: m,
            alphabet: self.alphabet,
            motif_type,
            icscore,
            nsites,
            pseudocount: self.pseudocount,
            bkg,
            bkgsites: notnan(self.bkgsites),
            consensus,
            strand: self.strand,
            pval: notnan(self.pval),
            qval: notnan(self.qval),
            eval: notnan(self.eval),
            extrainfo: self.extrainfo,
            gap: self.gap,
        })
    }
}

/// Format a number for printing: shortest round-trip form,
/// switching to scientific notation for very small or large values.
fn fmt_num(x: f64) -> String {
    if x.is_infinite() {
        return if x > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    let a = x.abs();
    if a != 0.0 && !(1e-4..1e15).contains(&a) {
        format!("{x:e}")
    } else {
        format!("{x}")
    }
}

fn round_to(x: f64, digits: i32) -> f64 {
    let f = 10f64.powi(digits);
    (x * f).round() / f
}

fn shorten(s: &str) -> String {
    if s.chars().count() > 40 {
        let mut out: String = s.chars().take(40).collect();
        out.push_str("...");
        out
    } else {
        s.to_string()
    }
}

/// A printout of the motif's metadata followed by its matrix.
impl fmt::Display for Motif {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let p = &self.0;
        writeln!(f)?;
        writeln!(f, "       Motif name:   {}", shorten(&p.name))?;
        if let Some(a) = &p.altname {
            writeln!(f, "   Alternate name:   {}", shorten(a))?;
        }
        if let Some(a) = &p.family {
            writeln!(f, "           Family:   {}", shorten(a))?;
        }
        if let Some(a) = &p.organism {
            writeln!(f, "         Organism:   {}", shorten(a))?;
        }
        writeln!(f, "         Alphabet:   {}", shorten(&p.alphabet.name()))?;
        writeln!(f, "             Type:   {}", p.motif_type)?;
        if p.alphabet.is_nucleotide() {
            writeln!(f, "          Strands:   {}", p.strand)?;
        }
        if !p.icscore.is_nan() {
            writeln!(
                f,
                "         Total IC:   {}",
                fmt_num(round_to(p.icscore, 2))
            )?;
        }
        writeln!(
            f,
            "      Pseudocount:   {}",
            fmt_num(round_to(p.pseudocount, 2))
        )?;
        if let Some(c) = &p.consensus {
            let c = if p.gap.isgapped {
                gapped_consensus(c, &p.gap.gaploc)
            } else {
                c.clone()
            };
            writeln!(f, "        Consensus:   {}", shorten(&c))?;
        }
        if let Some(n) = p.nsites {
            writeln!(f, "     Target sites:   {}", fmt_num(n))?;
        }
        if let Some(n) = p.bkgsites {
            writeln!(f, " Background sites:   {}", fmt_num(n))?;
        }
        if let Some(n) = p.pval {
            writeln!(f, "          P-value:   {}", fmt_num(n))?;
        }
        if let Some(n) = p.qval {
            writeln!(f, "          Q-value:   {}", fmt_num(n))?;
        }
        if let Some(n) = p.eval {
            writeln!(f, "          E-value:   {}", fmt_num(n))?;
        }
        if p.gap.isgapped {
            let locs: Vec<String> = p
                .gap
                .gaploc
                .iter()
                .map(|l| format!("{}-{}", l, l + 1))
                .collect();
            let sizes: Vec<String> = p
                .gap
                .mingap
                .iter()
                .zip(&p.gap.maxgap)
                .map(|(a, b)| format!("{a}-{b}"))
                .collect();
            writeln!(f, "    Gap locations:   {}", shorten(&locs.join(", ")))?;
            writeln!(f, "        Gap sizes:   {}", shorten(&sizes.join(", ")))?;
        }
        if !p.extrainfo.is_empty() {
            for (i, (k, v)) in p.extrainfo.iter().take(3).enumerate() {
                let s = if k.is_empty() {
                    v.clone()
                } else {
                    format!("[{k}] {v}")
                };
                let lead = if i == 0 {
                    "       Extra info:   "
                } else {
                    "                     "
                };
                writeln!(f, "{lead}{}", shorten(&s))?;
            }
            if p.extrainfo.len() > 3 {
                writeln!(f, "                     ...")?;
            }
        }
        writeln!(f)?;
        write_matrix(f, p)
    }
}

fn gapped_consensus(c: &str, gaploc: &[usize]) -> String {
    let chars: Vec<char> = c.chars().collect();
    let mut out = String::new();
    let mut start = 0;
    for &g in gaploc {
        let g = g.min(chars.len());
        out.extend(&chars[start..g]);
        out.push_str("..");
        start = g;
    }
    out.extend(&chars[start..]);
    out
}

fn write_matrix(f: &mut fmt::Formatter<'_>, p: &MotifParts) -> fmt::Result {
    let m = &p.matrix;
    let headers: Vec<String> = match &p.consensus {
        Some(c) => c.chars().map(String::from).collect(),
        None => (1..=m.ncol()).map(|i| i.to_string()).collect(),
    };
    let mut cols: Vec<Vec<String>> = Vec::with_capacity(m.ncol());
    for c in m.columns() {
        let rounded: Vec<f64> = c.iter().map(|&v| round_to(v, 2)).collect();
        let decimals = rounded
            .iter()
            .filter(|v| v.is_finite())
            .map(|&v| {
                if round_to(v, 0) == v {
                    0
                } else if round_to(v, 1) == v {
                    1
                } else {
                    2
                }
            })
            .max()
            .unwrap_or(0);
        cols.push(
            rounded
                .iter()
                .map(|&v| {
                    if v.is_finite() {
                        format!("{v:.decimals$}")
                    } else {
                        fmt_num(v)
                    }
                })
                .collect(),
        );
    }
    let widths: Vec<usize> = cols
        .iter()
        .zip(&headers)
        .map(|(c, h)| {
            c.iter()
                .map(|s| s.len())
                .chain([h.len()])
                .max()
                .unwrap_or(1)
        })
        .collect();
    let letters = p.alphabet.letters();
    write!(f, " ")?;
    for (h, w) in headers.iter().zip(&widths) {
        write!(f, " {h:>w$}")?;
    }
    writeln!(f)?;
    for (i, l) in letters.iter().enumerate() {
        write!(f, "{l}")?;
        for (c, w) in cols.iter().zip(&widths) {
            write!(f, " {:>w$}", c[i])?;
        }
        writeln!(f)?;
    }
    Ok(())
}
