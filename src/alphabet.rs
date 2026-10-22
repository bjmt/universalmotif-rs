//! Motif alphabets.

use crate::error::{Result, arg_err};

/// The four DNA letters, in row order.
pub const DNA: [char; 4] = ['A', 'C', 'G', 'T'];
/// The four RNA letters, in row order.
pub const RNA: [char; 4] = ['A', 'C', 'G', 'U'];
/// The twenty standard amino acids, in row order.
pub const AMINO_ACIDS: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

/// DNA IUPAC letters accepted in consensus strings.
pub const DNA_IUPAC: &str = "ACGTMRWSYKVHDBN";
/// Amino acid letters accepted in consensus strings, including ambiguity codes.
pub const AA_IUPAC: &str = "ACDEFGHIKLMNPQRSTVWYUOBJZX";

/// The alphabet of a motif. The rows of a motif matrix always follow the
/// order of [`Alphabet::letters`].
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Alphabet {
    Dna,
    Rna,
    Aa,
    /// Any other set of single-character letters, kept sorted and unique.
    Custom(Vec<char>),
}

impl Alphabet {
    /// Build a custom alphabet from arbitrary letters. The letters are sorted
    /// and deduplicated, and at least two distinct letters are required.
    pub fn custom(letters: impl IntoIterator<Item = char>) -> Result<Self> {
        let mut letters: Vec<char> = letters.into_iter().collect();
        letters.sort_unstable();
        letters.dedup();
        if letters.len() < 2 {
            return arg_err("a custom alphabet needs at least two distinct letters");
        }
        Ok(Alphabet::Custom(letters))
    }

    /// Parse an alphabet name: `"DNA"`, `"RNA"`, `"AA"`, or a string of
    /// letters making up a custom alphabet.
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "DNA" => Ok(Alphabet::Dna),
            "RNA" => Ok(Alphabet::Rna),
            "AA" => Ok(Alphabet::Aa),
            "custom" => {
                arg_err("alphabet 'custom' is not accepted; please provide the actual letters")
            }
            _ => Alphabet::custom(name.chars()),
        }
    }

    /// The letters of the alphabet, in row order.
    pub fn letters(&self) -> &[char] {
        match self {
            Alphabet::Dna => &DNA,
            Alphabet::Rna => &RNA,
            Alphabet::Aa => &AMINO_ACIDS,
            Alphabet::Custom(l) => l,
        }
    }

    pub fn len(&self) -> usize {
        self.letters().len()
    }

    pub fn is_empty(&self) -> bool {
        self.letters().is_empty()
    }

    /// Whether this is DNA or RNA.
    pub fn is_nucleotide(&self) -> bool {
        matches!(self, Alphabet::Dna | Alphabet::Rna)
    }

    /// The alphabet's name: `"DNA"`, `"RNA"`,
    /// `"AA"`, or the concatenated letters of a custom alphabet.
    pub fn name(&self) -> String {
        match self {
            Alphabet::Dna => "DNA".into(),
            Alphabet::Rna => "RNA".into(),
            Alphabet::Aa => "AA".into(),
            Alphabet::Custom(l) => l.iter().collect(),
        }
    }

    /// Row index of a letter, if it belongs to the alphabet.
    pub fn index_of(&self, letter: char) -> Option<usize> {
        self.letters().iter().position(|&l| l == letter)
    }
}

impl std::fmt::Display for Alphabet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.name())
    }
}
