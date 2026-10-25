//! Edge cases for motif construction and manipulation: alphabets, matrices,
//! the builder and validation, the constructors in `create`, type
//! conversions, reverse complements, trimming, rounding, gaps and merging.

use universalmotif::alphabet::{AMINO_ACIDS, DNA, RNA};
use universalmotif::convert::{
    consensus_letter, consensus_letter_aa, consensus_to_ppm, consensus_to_ppm_aa, convert_type,
    cwm_to_ppm, icm_to_ppm, pcm_to_ppm, position_icscore, ppm_to_icm, ppm_to_pcm, ppm_to_pwm,
    pwm_to_ppm,
};
use universalmotif::create::{
    CreateOptions, from_matrix, from_strings, from_typed_matrix, guess_matrix_alphabet,
    guess_strings_alphabet, random,
};
use universalmotif::merge::{
    MergeOptions, MergeSimilarOptions, merge_motifs, merge_motifs_detailed, merge_similar,
};
use universalmotif::ops::{
    TrimFrom, add_gap, normalize, reverse_complement, round_motif, switch_alph, trim_cwm,
    trim_motif, ungap,
};
use universalmotif::{Alphabet, Error, GapInfo, Matrix, Motif, MotifParts, MotifType, Strand};

// Helpers

fn arg_msg<T: std::fmt::Debug>(r: universalmotif::Result<T>) -> String {
    match r {
        Err(Error::Argument(m)) => m,
        other => panic!("expected an Argument error, got {other:?}"),
    }
}

fn unsupported_msg<T: std::fmt::Debug>(r: universalmotif::Result<T>) -> String {
    match r {
        Err(Error::Unsupported(m)) => m,
        other => panic!("expected an Unsupported error, got {other:?}"),
    }
}

fn invalid_msgs<T: std::fmt::Debug>(r: universalmotif::Result<T>) -> Vec<String> {
    match r {
        Err(Error::InvalidMotif(m)) => m,
        other => panic!("expected an InvalidMotif error, got {other:?}"),
    }
}

fn mat(cols: &[&[f64]]) -> Matrix {
    Matrix::from_columns(cols).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    if a.is_infinite() || b.is_infinite() {
        return a == b;
    }
    (a - b).abs() <= 1e-12 * f64::max(1.0, f64::max(a.abs(), b.abs()))
}

fn assert_vec_close(got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "lengths differ: {got:?} vs {want:?}");
    for (g, w) in got.iter().zip(want) {
        assert!(close(*g, *w), "{got:?} != {want:?}");
    }
}

fn assert_matrix_close(got: &Matrix, want: &Matrix) {
    assert_eq!((got.nrow(), got.ncol()), (want.nrow(), want.ncol()));
    assert_vec_close(got.as_slice(), want.as_slice());
}

/// Shannon information content of one probability column, in bits.
fn shannon_ic(p: &[f64]) -> f64 {
    let h: f64 = p.iter().filter(|&&x| x > 0.0).map(|&x| -x * x.log2()).sum();
    (p.len() as f64).log2() - h
}

/// A DNA PCM with ten sites: an A, then an even split of C and G.
fn as_pcm() -> Motif {
    Motif::builder(mat(&[&[10., 0., 0., 0.], &[0., 5., 5., 0.]]), Alphabet::Dna)
        .build()
        .unwrap()
}

/// A DNA PPM without a pseudocount: uniform, A, uniform, C, uniform, uniform.
fn trim_ppm() -> Motif {
    let u: &[f64] = &[0.25; 4];
    Motif::builder(
        mat(&[u, &[1., 0., 0., 0.], u, &[0., 1., 0., 0.], u, u]),
        Alphabet::Dna,
    )
    .pseudocount(0.0)
    .nsites(Some(100.0))
    .build()
    .unwrap()
}

fn sharp_ppm(pattern: &str, p: f64, name: &str, nsites: f64) -> Motif {
    let q = (1.0 - p) / 3.0;
    let cols: Vec<Vec<f64>> = pattern
        .chars()
        .map(|c| {
            let i = Alphabet::Dna.index_of(c).unwrap();
            let mut col = vec![q; 4];
            col[i] = p;
            col
        })
        .collect();
    Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
        .name(name)
        .motif_type(Some(MotifType::Ppm))
        .nsites(Some(nsites))
        .pseudocount(0.0)
        .build()
        .unwrap()
}

// Alphabets

#[test]
fn custom_alphabet_is_sorted_and_deduplicated() {
    let a = Alphabet::custom("CBAAB".chars()).unwrap();
    assert_eq!(a, Alphabet::Custom(vec!['A', 'B', 'C']));
    assert_eq!(a.letters(), &['A', 'B', 'C']);
    assert_eq!(a.len(), 3);
    assert!(!a.is_empty());
    assert!(!a.is_nucleotide());
    assert_eq!(a.name(), "ABC");
    assert_eq!(a.to_string(), "ABC");
    assert_eq!(a.index_of('C'), Some(2));
    assert_eq!(a.index_of('D'), None);
}

#[test]
fn custom_alphabet_needs_two_distinct_letters() {
    let want = "a custom alphabet needs at least two distinct letters";
    assert_eq!(arg_msg(Alphabet::custom("".chars())), want);
    assert_eq!(arg_msg(Alphabet::custom("X".chars())), want);
    assert_eq!(arg_msg(Alphabet::custom("XXXX".chars())), want);
}

#[test]
fn builtin_alphabets_report_their_letters() {
    assert_eq!(Alphabet::Dna.letters(), &DNA);
    assert_eq!(Alphabet::Rna.letters(), &RNA);
    assert_eq!(Alphabet::Aa.letters(), &AMINO_ACIDS);
    assert_eq!(Alphabet::Aa.len(), 20);
    for (a, nuc, name) in [
        (Alphabet::Dna, true, "DNA"),
        (Alphabet::Rna, true, "RNA"),
        (Alphabet::Aa, false, "AA"),
    ] {
        assert_eq!(a.is_nucleotide(), nuc);
        assert!(!a.is_empty());
        assert_eq!(a.to_string(), name);
        assert_eq!(name.parse::<Alphabet>().unwrap(), a);
    }
    assert_eq!(Alphabet::Rna.index_of('U'), Some(3));
    assert_eq!(Alphabet::Rna.index_of('T'), None);
    assert_eq!(Alphabet::Aa.index_of('Y'), Some(19));
}

#[test]
fn alphabet_parsing_rejects_custom_and_single_letters() {
    assert_eq!(
        arg_msg("custom".parse::<Alphabet>()),
        "alphabet 'custom' is not accepted; please provide the actual letters"
    );
    assert_eq!(
        arg_msg("Q".parse::<Alphabet>()),
        "a custom alphabet needs at least two distinct letters"
    );
    // Names are case sensitive, so lower case "dna" is a custom alphabet.
    assert_eq!(
        "dna".parse::<Alphabet>().unwrap(),
        Alphabet::Custom(vec!['a', 'd', 'n'])
    );
}

// Matrices

#[test]
fn matrix_construction_errors() {
    assert_eq!(
        arg_msg(Matrix::from_vec(0, 3, vec![])),
        "a motif matrix needs at least one row"
    );
    assert_eq!(
        arg_msg(Matrix::from_vec(2, 3, vec![0.0; 5])),
        "matrix data has 5 values, expected 6 (2 rows x 3 columns)"
    );
    assert_eq!(
        arg_msg(Matrix::from_columns::<Vec<f64>>(&[])),
        "a motif matrix needs at least one column"
    );
    assert_eq!(
        arg_msg(Matrix::from_columns(&[vec![1.0, 2.0], vec![3.0]])),
        "all matrix columns must have the same length"
    );
    assert_eq!(
        arg_msg(Matrix::from_columns(&[Vec::<f64>::new()])),
        "a motif matrix needs at least one row"
    );
}

#[test]
fn matrix_accessors_are_position_major() {
    let mut m = Matrix::from_vec(2, 3, vec![1., 2., 3., 4., 5., 6.]).unwrap();
    assert_eq!((m.nrow(), m.ncol()), (2, 3));
    assert_eq!(m.get(1, 0), 2.0);
    assert_eq!(m.get(0, 2), 5.0);
    assert_eq!(m.col(1), &[3., 4.]);
    assert_eq!(m.col_sums(), vec![3., 7., 11.]);
    m.col_mut(2)[1] = 10.0;
    assert_eq!(m.get(1, 2), 10.0);
    let cols: Vec<&[f64]> = m.columns().collect();
    assert_eq!(cols, vec![&[1., 2.][..], &[3., 4.], &[5., 10.]]);
    assert_eq!(m.columns().len(), 3);

    let doubled = m.map_columns(|c| c.iter().map(|v| v * 2.0).collect());
    assert_eq!(doubled.as_slice(), &[2., 4., 6., 8., 10., 20.]);

    let s = m.slice_columns(1..3);
    assert_eq!((s.nrow(), s.ncol()), (2, 2));
    assert_eq!(s.as_slice(), &[3., 4., 5., 10.]);
    let empty = m.slice_columns(1..1);
    assert_eq!((empty.nrow(), empty.ncol()), (2, 0));

    assert_eq!(m.reversed().as_slice(), &[10., 5., 4., 3., 2., 1.]);
    assert_eq!(m.reversed().reversed(), m);
    assert_eq!(m.clone().into_vec(), vec![1., 2., 3., 4., 5., 10.]);
}

#[test]
fn zero_column_matrix_is_allowed_but_not_as_a_motif() {
    let m = Matrix::from_vec(4, 0, vec![]).unwrap();
    assert_eq!(m.ncol(), 0);
    assert_eq!(m.columns().len(), 0);
    assert!(m.col_sums().is_empty());
    assert_eq!(
        arg_msg(Motif::builder(m, Alphabet::Dna).build()),
        "motif must have at least one column"
    );
}

// Type, strand and trim names

#[test]
fn motif_type_and_strand_names_round_trip() {
    for t in MotifType::ALL {
        assert_eq!(t.as_str().parse::<MotifType>().unwrap(), t);
        assert_eq!(t.to_string(), t.as_str());
    }
    assert_eq!(
        arg_msg("pwm".parse::<MotifType>()),
        "incorrect type 'pwm' (must be one of PCM, PPM, PWM, ICM, CWM)"
    );
    for s in [Strand::Plus, Strand::Minus, Strand::Both] {
        assert_eq!(s.as_str().parse::<Strand>().unwrap(), s);
        assert_eq!(s.to_string(), s.as_str());
    }
    assert_eq!("-+".parse::<Strand>().unwrap(), Strand::Both);
    assert_eq!("*".parse::<Strand>().unwrap(), Strand::Both);
    assert_eq!(Strand::default(), Strand::Both);
    assert_eq!(
        arg_msg("both".parse::<Strand>()),
        "incorrect strand 'both' (must be one of +, -, +-)"
    );
    for t in [TrimFrom::Both, TrimFrom::Left, TrimFrom::Right] {
        assert_eq!(t.as_str().parse::<TrimFrom>().unwrap(), t);
        assert_eq!(t.to_string(), t.as_str());
    }
    assert_eq!(TrimFrom::default(), TrimFrom::Both);
    assert_eq!(
        arg_msg("middle".parse::<TrimFrom>()),
        "incorrect trim_from 'middle' (must be one of both, left, right)"
    );
}

// Validation of motif parts

fn base_parts() -> MotifParts {
    as_pcm().into_parts()
}

#[test]
fn validate_reports_every_gap_problem() {
    let mut p = base_parts();
    p.gap = GapInfo {
        isgapped: true,
        gaploc: vec![0, 2],
        mingap: vec![1],
        maxgap: vec![0],
    };
    assert_eq!(
        p.validate(),
        vec![
            "gaploc, mingap and maxgap should all be the same length",
            "position 0 gaps or less are not allowed",
            "gap location values should not exceed motif size - 1",
            "mingap cannot be greater than the corresponding maxgap",
            "maxgap must be greater than 0",
        ]
    );
    let err = Motif::from_parts(p).unwrap_err();
    assert!(
        err.to_string().starts_with(
            "invalid motif:\ngaploc, mingap and maxgap should all be the same length\n"
        )
    );
}

#[test]
fn validate_reports_matrix_problems() {
    let mut p = base_parts();
    p.matrix = Matrix::from_vec(4, 0, vec![]).unwrap();
    assert_eq!(
        p.validate(),
        vec![
            "motif must have at least one column",
            "consensus string must have the same number of letters as motif positions",
        ]
    );

    let mut p = base_parts();
    p.matrix = mat(&[&[10., 0., 0., 0.], &[0., f64::NAN, 5., 0.]]);
    p.pseudocount = -1.0;
    assert_eq!(
        p.validate(),
        vec![
            "motif cannot contain NaN values",
            "pseudocount must be a finite, non-negative number",
        ]
    );

    let mut p = base_parts();
    p.pseudocount = f64::INFINITY;
    assert_eq!(
        p.validate(),
        vec!["pseudocount must be a finite, non-negative number"]
    );

    let mut p = base_parts();
    p.consensus = Some("A".into());
    assert_eq!(
        p.validate(),
        vec!["consensus string must have the same number of letters as motif positions"]
    );
}

#[test]
fn validate_checks_rows_against_the_alphabet() {
    let mut p = base_parts();
    p.alphabet = Alphabet::Aa;
    assert_eq!(p.validate(), vec!["AA motifs must have 20 rows"]);

    let mut p = base_parts();
    p.alphabet = Alphabet::custom("XYZ".chars()).unwrap();
    assert_eq!(
        p.validate(),
        vec!["alphabet length does not match number of rows in motif"]
    );

    let mut p = base_parts();
    p.alphabet = Alphabet::Rna;
    p.matrix = mat(&[&[10., 0., 0.], &[0., 5., 5.]]);
    p.bkg = vec![0.5, 0.25, 0.25];
    assert_eq!(p.validate(), vec!["DNA/RNA motifs must have 4 rows"]);
}

#[test]
fn validate_checks_the_background() {
    let mut p = base_parts();
    p.bkg = vec![0.5, 0.5, 0.0];
    assert_eq!(
        p.validate(),
        vec!["bkg must have one value per alphabet letter"]
    );

    let mut p = base_parts();
    p.bkg = vec![-0.5, 0.5, 0.5, 0.5];
    assert_eq!(p.validate(), vec!["bkg does not allow values less than 0"]);

    let mut p = base_parts();
    p.bkg = vec![1.5, -0.5, 0.0, 0.0];
    assert_eq!(
        p.validate(),
        vec![
            "bkg does not allow values less than 0",
            "bkg does not allow values higher than 1",
        ]
    );

    let mut p = base_parts();
    p.bkg = vec![0.5; 4];
    assert_eq!(p.validate(), vec!["bkg probabilities must add up to 1"]);

    // The sum is allowed to be off by up to 0.01.
    let mut p = base_parts();
    p.bkg = vec![0.2524; 4];
    assert!(p.validate().is_empty());
}

#[test]
fn validate_checks_each_motif_type() {
    let mut p = base_parts();
    p.matrix = mat(&[&[10., 0., 0., 0.], &[0., 5., 4., 0.]]);
    assert_eq!(p.validate(), vec!["PCM columns must all have the same sum"]);

    let mut p = base_parts();
    p.matrix = mat(&[&[0.5, 0.5, 0., 0.], &[0., 0.5, 0.5, 0.]]);
    assert_eq!(
        p.validate(),
        vec!["PCM values must be whole counts, not between 0 and 1"]
    );

    let mut p = base_parts();
    p.matrix = mat(&[&[0.; 4], &[0.; 4]]);
    assert_eq!(
        p.validate(),
        vec!["PCM columns must have at least one count"]
    );

    let mut p = base_parts();
    p.motif_type = MotifType::Ppm;
    p.matrix = mat(&[&[0.5, 0.5, 0., 0.], &[1.2, -0.2, 0., 0.]]);
    assert_eq!(p.validate(), vec!["PPM values cannot be negative"]);

    let mut p = base_parts();
    p.motif_type = MotifType::Ppm;
    p.matrix = mat(&[&[0.5, 0.4, 0., 0.], &[0., 1., 0., 0.]]);
    assert_eq!(p.validate(), vec!["PPM columns must sum to 1"]);

    let mut p = base_parts();
    p.motif_type = MotifType::Icm;
    p.matrix = mat(&[&[0.5, -0.1, 0., 0.], &[0., 1., 0., 0.]]);
    assert_eq!(p.validate(), vec!["ICM values cannot be negative"]);

    // Weights are free to be anything but NaN.
    for t in [MotifType::Pwm, MotifType::Cwm] {
        let mut p = base_parts();
        p.motif_type = t;
        p.matrix = mat(&[&[-3., 7., f64::NEG_INFINITY, 0.], &[0., 0., 0., 0.]]);
        assert!(p.validate().is_empty());
    }
}

#[test]
fn recompute_derived_follows_matrix_changes() {
    let mut p = base_parts();
    p.matrix = mat(&[&[0., 0., 0., 10.], &[0., 0., 10., 0.]]);
    p.pseudocount = 0.0;
    p.recompute_derived();
    assert_eq!(p.consensus.as_deref(), Some("TG"));
    assert_eq!(p.icscore, 4.0);

    p.alphabet = Alphabet::custom("WXYZ".chars()).unwrap();
    p.motif_type = MotifType::Cwm;
    p.recompute_derived();
    assert_eq!(p.consensus, None);
    assert!(p.icscore.is_nan());
    // Parts holding NaN totals still compare equal.
    assert_eq!(p.clone(), p);
}

// The builder

#[test]
fn builder_rejects_unusable_matrices() {
    let nan = mat(&[&[1., f64::NAN, 0., 0.]]);
    assert_eq!(
        arg_msg(Motif::builder(nan, Alphabet::Dna).build()),
        "matrix cannot have NaN values"
    );
    let three = mat(&[&[1., 0., 0.]]);
    assert_eq!(
        arg_msg(Motif::builder(three.clone(), Alphabet::Rna).build()),
        "DNA/RNA motifs must have 4 rows"
    );
    assert_eq!(
        arg_msg(Motif::builder(three, Alphabet::Aa).build()),
        "AA motifs must have 20 rows"
    );
    let four = mat(&[&[1., 0., 0., 0.]]);
    assert_eq!(
        arg_msg(Motif::builder(four.clone(), Alphabet::custom("XY".chars()).unwrap()).build()),
        "alphabet length does not match matrix rows"
    );
    assert_eq!(
        arg_msg(
            Motif::builder(four, Alphabet::Dna)
                .bkg(Some(vec![0.5, 0.25, 0.25]))
                .build()
        ),
        "bkg has 3 values but the alphabet has 4 letters"
    );
}

#[test]
fn builder_guesses_the_type() {
    let guess = |cols: &[&[f64]]| {
        Motif::builder(mat(cols), Alphabet::Dna)
            .build()
            .unwrap()
            .motif_type
    };
    assert_eq!(guess(&[&[2., 0., 0., 0.]]), MotifType::Pcm);
    assert_eq!(guess(&[&[0.5, 0.5, 0., 0.]]), MotifType::Ppm);
    assert_eq!(guess(&[&[0.5, 0.5, 0.5, 0.5]]), MotifType::Icm);
    assert_eq!(guess(&[&[-1., 1., 0., 0.]]), MotifType::Pwm);

    // Probabilities rounded to two decimals are renormalized.
    let m = Motif::builder(mat(&[&[0.5, 0.485, 0.0, 0.0]]), Alphabet::Dna)
        .build()
        .unwrap();
    assert_eq!(m.motif_type, MotifType::Ppm);
    assert_vec_close(m.matrix.col(0), &[0.5 / 0.985, 0.485 / 0.985, 0.0, 0.0]);

    // Uneven whole-number columns cannot be guessed without nsites.
    let err = arg_msg(
        Motif::builder(mat(&[&[2., 0., 0., 0.], &[3., 0., 0., 0.]]), Alphabet::Dna).build(),
    );
    assert!(err.starts_with("the matrix holds whole numbers, so it was read as counts"));
}

#[test]
fn builder_renormalizes_a_declared_ppm() {
    let m = Motif::builder(mat(&[&[2., 2., 0., 0.]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .build()
        .unwrap();
    assert_eq!(m.matrix.col(0), &[0.5, 0.5, 0.0, 0.0]);

    // Negative values are not renormalized, so validation catches them.
    let r = Motif::builder(mat(&[&[2., -1., 0., 0.]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .build();
    assert_eq!(invalid_msgs(r), vec!["PPM values cannot be negative"]);
}

#[test]
fn builder_rescales_counts_to_nsites() {
    let m = Motif::builder(mat(&[&[1., 2., 3., 4.]]), Alphabet::Dna)
        .nsites(Some(20.0))
        .build()
        .unwrap();
    assert_eq!(m.matrix.col(0), &[2., 4., 6., 8.]);
    assert_eq!(m.nsites, Some(20.0));

    // 10/3 rounds to 3 three times, and the missing count goes to the first
    // largest value.
    let xyz = Alphabet::custom("XYZ".chars()).unwrap();
    let m = Motif::builder(mat(&[&[1., 1., 1.], &[2., 2., 0.]]), xyz)
        .nsites(Some(10.0))
        .build()
        .unwrap();
    assert_eq!(m.matrix.col(0), &[4., 3., 3.]);
    assert_eq!(m.matrix.col(1), &[5., 5., 0.]);
    assert_eq!(m.consensus, None);

    // A NaN nsites is ignored, so the column sum is used.
    let m = Motif::builder(mat(&[&[1., 2., 3., 4.]]), Alphabet::Dna)
        .nsites(Some(f64::NAN))
        .build()
        .unwrap();
    assert_eq!(m.nsites, Some(10.0));
    assert_eq!(m.matrix.col(0), &[1., 2., 3., 4.]);
}

#[test]
fn builder_fills_defaults_and_drops_empty_metadata() {
    let m = Motif::builder(mat(&[&[1., 0., 0., 0.]]), Alphabet::Dna)
        .name("")
        .altname(Some(String::new()))
        .family(Some(String::new()))
        .organism(Some("human".into()))
        .pval(Some(f64::NAN))
        .qval(Some(0.5))
        .eval(Some(f64::NAN))
        .bkgsites(Some(f64::NAN))
        .build()
        .unwrap();
    assert_eq!(m.name, "[no name]");
    assert_eq!(m.altname, None);
    assert_eq!(m.family, None);
    assert_eq!(m.organism.as_deref(), Some("human"));
    assert_eq!(
        (m.pval, m.qval, m.eval, m.bkgsites),
        (None, Some(0.5), None, None)
    );
    assert_eq!(m.bkg, vec![0.25; 4]);
    assert_eq!(m.pseudocount, 1.0);
    assert_eq!(m.strand, Strand::Both);
    assert_eq!(m.width(), 1);

    let named = Motif::builder(mat(&[&[1., 0., 0., 0.]]), Alphabet::Dna)
        .build()
        .unwrap();
    assert_eq!(named.name, "new motif");
}

#[test]
fn builder_total_ic_can_be_given_except_for_cwms() {
    let m = as_pcm();
    // Pseudocount 1 spread over four letters, out of 11.
    let want = shannon_ic(&[10.25 / 11., 0.25 / 11., 0.25 / 11., 0.25 / 11.])
        + shannon_ic(&[0.25 / 11., 5.25 / 11., 5.25 / 11., 0.25 / 11.]);
    assert!((m.icscore - want).abs() < 1e-12, "{} vs {want}", m.icscore);
    assert_eq!(m.consensus.as_deref(), Some("AS"));

    let given = Motif::builder(mat(&[&[10., 0., 0., 0.]]), Alphabet::Dna)
        .icscore(Some(5.0))
        .build()
        .unwrap();
    assert_eq!(given.icscore, 5.0);
    let nan = Motif::builder(mat(&[&[10., 0., 0., 0.]]), Alphabet::Dna)
        .pseudocount(0.0)
        .icscore(Some(f64::NAN))
        .build()
        .unwrap();
    assert_eq!(nan.icscore, 2.0);
    let cwm = Motif::builder(mat(&[&[1., -1., 0., 0.]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Cwm))
        .icscore(Some(5.0))
        .build()
        .unwrap();
    assert!(cwm.icscore.is_nan());
    // |x| / sum|x| gives an even A/C split.
    assert_eq!(cwm.consensus.as_deref(), Some("M"));
}

// Printing

#[test]
fn display_of_a_plain_pcm() {
    let want = "
       Motif name:   new motif
         Alphabet:   DNA
             Type:   PCM
          Strands:   +-
         Total IC:   2.27
      Pseudocount:   1
        Consensus:   AS
     Target sites:   10

   A S
A 10 0
C  0 5
G  0 5
T  0 0
";
    assert_eq!(as_pcm().to_string(), want);
}

#[test]
fn display_of_every_field() {
    let m = Motif::builder(
        mat(&[
            &[0.25, 0.25, 0.5, 0.0],
            &[1., 0., 0., 0.],
            &[0., 0., 0., 1.],
        ]),
        Alphabet::Dna,
    )
    .name("a".repeat(45))
    .altname(Some("alt".into()))
    .family(Some("fam".into()))
    .organism(Some("org".into()))
    .strand(Strand::Plus)
    .icscore(Some(3.14659))
    .pseudocount(0.5)
    .nsites(Some(10.0))
    .bkgsites(Some(1e15))
    .pval(Some(1e-5))
    .qval(Some(f64::INFINITY))
    .eval(Some(0.5))
    .gap(GapInfo {
        isgapped: true,
        gaploc: vec![1],
        mingap: vec![2],
        maxgap: vec![4],
    })
    .extrainfo(vec![
        (String::new(), "plain".into()),
        ("k".into(), "v".into()),
        ("k2".into(), "v2".into()),
        ("k3".into(), "v3".into()),
    ])
    .build()
    .unwrap();
    let want = format!(
        "
       Motif name:   {}...
   Alternate name:   alt
           Family:   fam
         Organism:   org
         Alphabet:   DNA
             Type:   PPM
          Strands:   +
         Total IC:   3.15
      Pseudocount:   0.5
        Consensus:   N..AT
     Target sites:   10
 Background sites:   1e15
          P-value:   1e-5
          Q-value:   Inf
          E-value:   0.5
    Gap locations:   1-2
        Gap sizes:   2-4
       Extra info:   plain
                     [k] v
                     [k2] v2
                     ...

     N A T
A 0.25 1 0
C 0.25 0 0
G 0.50 0 0
T 0.00 0 1
",
        "a".repeat(40)
    );
    assert_eq!(m.to_string(), want);
}

#[test]
fn display_of_a_custom_alphabet_pwm_with_infinite_weights() {
    let xy = Alphabet::custom("YX".chars()).unwrap();
    let m = Motif::builder(mat(&[&[f64::NEG_INFINITY, 0.5], &[1., -1.]]), xy)
        .motif_type(Some(MotifType::Pwm))
        .build()
        .unwrap();
    // Column 0 becomes [0, 1] and column 1 becomes [0.8, 0.2].
    let want_ic = 1.0 + shannon_ic(&[0.8, 0.2]);
    assert!((m.icscore - want_ic).abs() < 1e-12);
    let want = "
       Motif name:   new motif
         Alphabet:   XY
             Type:   PWM
         Total IC:   1.28
      Pseudocount:   1

     1  2
X -Inf  1
Y  0.5 -1
";
    assert_eq!(m.to_string(), want);

    let cwm = Motif::builder(mat(&[&[1., -1., 0., 0.]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Cwm))
        .build()
        .unwrap();
    assert!(!cwm.to_string().contains("Total IC"));
}

// Creating motifs from matrices

#[test]
fn from_matrix_sorts_rows_by_letter_and_guesses_the_alphabet() {
    let pcm = CreateOptions {
        motif_type: MotifType::Pcm,
        ..Default::default()
    };
    let m = from_matrix(
        mat(&[&[1., 2., 3., 4.]]),
        Some(&['T', 'G', 'C', 'A']),
        None,
        &pcm,
    )
    .unwrap();
    assert_eq!(m.alphabet, Alphabet::Dna);
    assert_eq!(m.matrix.col(0), &[4., 3., 2., 1.]);
    assert_eq!(m.nsites, Some(10.0));
    assert_eq!(m.name, "motif");

    let m = from_matrix(
        mat(&[&[1., 2., 3., 4.]]),
        Some(&['U', 'G', 'C', 'A']),
        None,
        &pcm,
    )
    .unwrap();
    assert_eq!(m.alphabet, Alphabet::Rna);
    assert_eq!(m.consensus.as_deref(), Some("N"));

    let rev: Vec<char> = AMINO_ACIDS.iter().rev().copied().collect();
    let col: Vec<f64> = (1..=20).map(f64::from).collect();
    let m = from_matrix(
        Matrix::from_columns(&[col]).unwrap(),
        Some(&rev),
        None,
        &pcm,
    )
    .unwrap();
    assert_eq!(m.alphabet, Alphabet::Aa);
    let want: Vec<f64> = (1..=20).rev().map(f64::from).collect();
    assert_eq!(m.matrix.col(0), &want[..]);

    // Four letters that are not nucleotides make a custom alphabet.
    let m = from_matrix(
        mat(&[&[1., 2., 3., 4.]]),
        Some(&['Q', 'W', 'E', 'R']),
        None,
        &pcm,
    )
    .unwrap();
    assert_eq!(m.alphabet, Alphabet::Custom(vec!['E', 'Q', 'R', 'W']));
    assert_eq!(m.matrix.col(0), &[3., 1., 4., 2.]);

    // A given custom alphabet is checked against the row letters.
    let xyz = Alphabet::custom("XYZ".chars()).unwrap();
    let m = from_matrix(
        mat(&[&[1., 2., 3.]]),
        Some(&['Z', 'X', 'Y']),
        Some(xyz),
        &pcm,
    )
    .unwrap();
    assert_eq!(m.matrix.col(0), &[2., 3., 1.]);
}

#[test]
fn guess_matrix_alphabet_needs_exact_sets() {
    assert_eq!(
        guess_matrix_alphabet(&['A', 'C', 'G', 'T']).unwrap(),
        Alphabet::Dna
    );
    assert_eq!(
        guess_matrix_alphabet(&['A', 'C', 'G', 'U']).unwrap(),
        Alphabet::Rna
    );
    assert_eq!(
        guess_matrix_alphabet(&['A', 'C', 'G']).unwrap(),
        Alphabet::Custom(vec!['A', 'C', 'G'])
    );
    assert_eq!(guess_matrix_alphabet(&AMINO_ACIDS).unwrap(), Alphabet::Aa);
    assert_eq!(
        arg_msg(guess_matrix_alphabet(&['A'])),
        "a custom alphabet needs at least two distinct letters"
    );
}

#[test]
fn from_matrix_errors() {
    let opts = CreateOptions::default();
    let m4 = || mat(&[&[1., 2., 3., 4.]]);
    assert_eq!(
        arg_msg(from_matrix(
            mat(&[&[1., f64::NAN, 3., 4.]]),
            None,
            Some(Alphabet::Dna),
            &opts
        )),
        "matrix cannot have NaN values"
    );
    assert_eq!(
        arg_msg(from_matrix(m4(), Some(&['A', 'C', 'G']), None, &opts)),
        "number of row letters does not match matrix rows"
    );
    assert_eq!(
        arg_msg(from_matrix(m4(), Some(&['A', 'C', 'A', 'T']), None, &opts)),
        "row letters must be unique"
    );
    assert_eq!(
        arg_msg(from_matrix(
            mat(&[&[1., 2., 3.]]),
            Some(&['X', 'Y', 'W']),
            Some(Alphabet::custom("XYZ".chars()).unwrap()),
            &opts
        )),
        "row letters do not match the provided alphabet"
    );
    assert_eq!(
        arg_msg(from_matrix(m4(), None, None, &opts)),
        "please provide an alphabet or a matrix with row letters"
    );
    assert_eq!(
        arg_msg(from_matrix(
            mat(&[&[1., 2., 3.]]),
            None,
            Some(Alphabet::Dna),
            &opts
        )),
        "incorrect number of rows"
    );
    assert_eq!(
        arg_msg(from_matrix(m4(), None, Some(Alphabet::Aa), &opts)),
        "incorrect number of rows"
    );
    assert_eq!(
        arg_msg(from_matrix(
            m4(),
            None,
            Some(Alphabet::custom("XY".chars()).unwrap()),
            &opts
        )),
        "alphabet length does not match number of rows"
    );
    for n in [0.5, 2.5, f64::INFINITY, f64::NAN, 0.0] {
        let bad = CreateOptions {
            nsites: Some(n),
            ..Default::default()
        };
        assert_eq!(
            arg_msg(from_matrix(m4(), None, Some(Alphabet::Dna), &bad)),
            "nsites must be a positive integer",
            "nsites {n}"
        );
    }
}

#[test]
fn from_typed_matrix_skips_type_detection() {
    let opts = CreateOptions::default();
    // [1, 0, 0, 0] would be read as a count, but as log2 odds against a
    // uniform background it gives 2, 1, 1, 1 times 0.25, renormalized.
    let detected =
        from_matrix(mat(&[&[1., 0., 0., 0.]]), None, Some(Alphabet::Dna), &opts).unwrap();
    assert_eq!(detected.matrix.col(0), &[1., 0., 0., 0.]);
    let typed = from_typed_matrix(
        mat(&[&[1., 0., 0., 0.]]),
        None,
        Some(Alphabet::Dna),
        MotifType::Pwm,
        &opts,
    )
    .unwrap();
    assert_eq!(typed.motif_type, MotifType::Ppm);
    assert_vec_close(typed.matrix.col(0), &[0.4, 0.2, 0.2, 0.2]);

    assert_eq!(
        arg_msg(from_typed_matrix(
            mat(&[&[1., 0., 0., 0.]]),
            None,
            None,
            MotifType::Pwm,
            &opts
        )),
        "please provide an alphabet or a matrix with row letters"
    );
}

#[test]
fn from_matrix_keeps_cwms_as_given() {
    let opts = CreateOptions {
        motif_type: MotifType::Cwm,
        ..Default::default()
    };
    let m = from_matrix(
        mat(&[&[0.5, -0.25, 0., 0.], &[0., 0., 0., 2.]]),
        None,
        Some(Alphabet::Dna),
        &opts,
    )
    .unwrap();
    assert_eq!(m.motif_type, MotifType::Cwm);
    assert_eq!(m.matrix.col(0), &[0.5, -0.25, 0., 0.]);
    assert!(m.icscore.is_nan());
    assert_eq!(m.consensus.as_deref(), Some("MT"));
}

#[test]
fn create_options_metadata_reaches_the_motif() {
    let opts = CreateOptions {
        motif_type: MotifType::Pcm,
        name: "m1".into(),
        nsites: Some(4.0),
        bkg: Some(vec![0.3, 0.2, 0.2, 0.3]),
        altname: Some("alt".into()),
        family: Some("bZIP".into()),
        organism: Some("yeast".into()),
        bkgsites: Some(1000.0),
        strand: Strand::Plus,
        pval: Some(0.01),
        qval: Some(0.02),
        eval: Some(0.03),
        extrainfo: vec![("source".into(), "test".into())],
        ..Default::default()
    };
    let m = from_matrix(mat(&[&[1., 1., 1., 1.]]), None, Some(Alphabet::Dna), &opts).unwrap();
    assert_eq!(m.name, "m1");
    assert_eq!(m.altname.as_deref(), Some("alt"));
    assert_eq!(m.family.as_deref(), Some("bZIP"));
    assert_eq!(m.organism.as_deref(), Some("yeast"));
    assert_eq!(m.bkg, vec![0.3, 0.2, 0.2, 0.3]);
    assert_eq!(m.bkgsites, Some(1000.0));
    assert_eq!(m.strand, Strand::Plus);
    assert_eq!(
        (m.pval, m.qval, m.eval),
        (Some(0.01), Some(0.02), Some(0.03))
    );
    assert_eq!(
        m.extrainfo,
        vec![("source".to_string(), "test".to_string())]
    );
    assert_eq!(m.pseudocount, 0.0);
    assert!(!m.gap.isgapped);
}

// Creating motifs from strings

#[test]
fn guess_strings_alphabet_order() {
    assert_eq!(guess_strings_alphabet(&["ACGN"]).unwrap(), Alphabet::Dna);
    assert_eq!(
        guess_strings_alphabet(&["MKV", "RYS"]).unwrap(),
        Alphabet::Dna
    );
    assert_eq!(guess_strings_alphabet(&["ACGU"]).unwrap(), Alphabet::Rna);
    assert_eq!(guess_strings_alphabet(&["EFQ"]).unwrap(), Alphabet::Aa);
    // U with an N fits no nucleotide alphabet, but both are amino acid codes.
    assert_eq!(guess_strings_alphabet(&["ACGUN"]).unwrap(), Alphabet::Aa);
    assert_eq!(
        guess_strings_alphabet(&["xyz", "zzx"]).unwrap(),
        Alphabet::Custom(vec!['x', 'y', 'z'])
    );
    assert_eq!(
        arg_msg(guess_strings_alphabet(&["11"])),
        "a custom alphabet needs at least two distinct letters"
    );
}

#[test]
fn from_strings_errors() {
    let opts = CreateOptions::default();
    assert_eq!(
        arg_msg(from_strings(&[], None, &opts)),
        "at least one string is required"
    );
    assert_eq!(
        arg_msg(from_strings(&["ACG", "AC"], None, &opts)),
        "all strings must be of equal length"
    );
    assert_eq!(
        arg_msg(from_strings(&[""], None, &opts)),
        "zero-length input strings are not allowed"
    );
    assert_eq!(
        arg_msg(from_strings(&["ACZX"], Some(Alphabet::Dna), &opts)),
        "the input has letters outside the DNA alphabet: XZ"
    );
    assert_eq!(
        arg_msg(from_strings(&["ACGT"], Some(Alphabet::Rna), &opts)),
        "the input has letters outside the RNA alphabet: T"
    );
    assert_eq!(
        arg_msg(from_strings(&["AC1"], Some(Alphabet::Aa), &opts)),
        "the input has letters outside the AA alphabet: 1"
    );
    let custom = arg_msg(from_strings(
        &["xyw"],
        Some(Alphabet::custom("xyz".chars()).unwrap()),
        &opts,
    ));
    assert!(custom.starts_with("the input has letters outside the "));
    assert!(custom.ends_with(" alphabet: w"));
    let dna = arg_msg(from_strings(&["ACGN", "ACGT"], None, &opts));
    assert!(
        dna.starts_with("DNA ambiguity letters are not accepted if input is more than one string")
    );
    let rna = arg_msg(from_strings(&["ACGU", "ACGT"], Some(Alphabet::Rna), &opts));
    assert_eq!(rna, "the input has letters outside the RNA alphabet: T");
    let aa = arg_msg(from_strings(&["MKB", "MKL"], None, &opts));
    assert!(
        aa.starts_with("AA ambiguity letters are not accepted if input is more than one string")
    );
    let bad = CreateOptions {
        nsites: Some(1.5),
        ..Default::default()
    };
    assert_eq!(
        arg_msg(from_strings(&["ACGT"], None, &bad)),
        "nsites must be a positive integer"
    );
}

#[test]
fn from_strings_custom_alphabet_error_message() {
    let msg = arg_msg(from_strings(
        &["xyw"],
        Some(Alphabet::custom("xyz".chars()).unwrap()),
        &CreateOptions::default(),
    ));
    assert_eq!(msg, "the input has letters outside the alphabet: w");
}

#[test]
fn from_strings_consensus_columns() {
    let opts = CreateOptions {
        motif_type: MotifType::Pcm,
        nsites: Some(1000.0),
        ..Default::default()
    };
    let m = from_strings(&["ACGTN"], None, &opts).unwrap();
    assert_eq!(m.alphabet, Alphabet::Dna);
    assert_eq!(m.consensus.as_deref(), Some("ACGTN"));
    assert_eq!(m.matrix.col(0), &[997., 1., 1., 1.]);
    assert_eq!(m.matrix.col(3), &[1., 1., 1., 997.]);
    assert_eq!(m.matrix.col(4), &[250.; 4]);
    assert_eq!(m.nsites, Some(1000.0));

    let m = from_strings(&["ACGU"], Some(Alphabet::Rna), &opts).unwrap();
    assert_eq!(m.consensus.as_deref(), Some("ACGU"));
    assert_eq!(m.matrix.col(3), &[1., 1., 1., 997.]);

    let m = from_strings(&["MB"], Some(Alphabet::Aa), &opts).unwrap();
    assert_eq!(m.consensus.as_deref(), Some("MB"));
    let m_row = Alphabet::Aa.index_of('M').unwrap();
    assert_eq!(m.matrix.get(m_row, 0), 981.0);
    assert_eq!(m.matrix.col_sums(), vec![1000.0, 1000.0]);
    assert_eq!(m.matrix.get(2, 1), 491.0); // D
    assert_eq!(m.matrix.get(11, 1), 491.0); // N
}

#[test]
fn from_strings_custom_alphabets() {
    // A single string over a custom alphabet is one-hot.
    let m = from_strings(&["abca"], None, &CreateOptions::default()).unwrap();
    assert_eq!(m.alphabet, Alphabet::Custom(vec!['a', 'b', 'c']));
    assert_eq!(m.motif_type, MotifType::Ppm);
    assert_eq!(m.consensus, None);
    assert_eq!(
        m.matrix.as_slice(),
        &[1., 0., 0., 0., 1., 0., 0., 0., 1., 1., 0., 0.]
    );

    // Several strings are counted, and nsites becomes their number.
    let pcm = CreateOptions {
        motif_type: MotifType::Pcm,
        ..Default::default()
    };
    let m = from_strings(
        &["ab", "ba", "aa"],
        Some(Alphabet::custom("ab".chars()).unwrap()),
        &pcm,
    )
    .unwrap();
    assert_eq!(m.nsites, Some(3.0));
    assert_eq!(m.matrix.as_slice(), &[2., 1., 2., 1.]);

    // Ambiguity letters can be counted as letters of their own.
    let m = from_strings(
        &["AN", "NA"],
        Some(Alphabet::custom("AN".chars()).unwrap()),
        &pcm,
    )
    .unwrap();
    assert_eq!(m.matrix.as_slice(), &[1., 1., 1., 1.]);
}

#[test]
fn from_strings_sites_round_trip_through_types() {
    let sites = ["TGACTCA", "TGAGTCA", "TGACTAA", "AGACTCA"];
    let pcm = from_strings(
        &sites,
        None,
        &CreateOptions {
            motif_type: MotifType::Pcm,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(pcm.nsites, Some(4.0));
    assert_eq!(pcm.matrix.col(0), &[1., 0., 0., 3.]);
    assert_eq!(pcm.matrix.col(3), &[0., 3., 1., 0.]);
    assert_eq!(pcm.consensus.as_deref(), Some("TGACTCA"));
    for t in [MotifType::Ppm, MotifType::Pwm, MotifType::Icm] {
        let other = from_strings(
            &sites,
            None,
            &CreateOptions {
                motif_type: t,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(other.motif_type, t);
        let back = convert_type(&other, MotifType::Pcm, None, false).unwrap();
        assert_eq!(back.matrix, pcm.matrix, "{t}");
    }
}

// Random motifs

#[test]
fn random_motifs_are_reproducible() {
    let opts = CreateOptions::default();
    let a = random(6, Alphabet::Dna, &opts, 42).unwrap();
    let b = random(6, Alphabet::Dna, &opts, 42).unwrap();
    let c = random(6, Alphabet::Dna, &opts, 43).unwrap();
    assert_eq!(a, b);
    assert_ne!(a.matrix, c.matrix);
    assert_eq!(a.width(), 6);
    assert_eq!(a.motif_type, MotifType::Ppm);
    let n = a.nsites.unwrap();
    assert!((50.0..=250.0).contains(&n) && n == n.floor(), "nsites {n}");
    for s in a.matrix.col_sums() {
        assert!((s - 1.0).abs() < 1e-9);
    }
    assert!((a.bkg.iter().sum::<f64>() - 1.0).abs() < 1e-9);
}

#[test]
fn random_motif_options() {
    let pcm = CreateOptions {
        motif_type: MotifType::Pcm,
        ..Default::default()
    };
    let m = random(4, Alphabet::Rna, &pcm, 1).unwrap();
    assert_eq!(m.nsites, Some(100.0));
    assert_eq!(m.matrix.col_sums(), vec![100.0; 4]);

    let pwm = CreateOptions {
        motif_type: MotifType::Pwm,
        ..Default::default()
    };
    assert_eq!(random(4, Alphabet::Dna, &pwm, 1).unwrap().nsites, None);

    // A letter with no background weight never appears.
    let skewed = CreateOptions {
        bkg: Some(vec![0.5, 0.5, 0.0, 0.0]),
        nsites: Some(20.0),
        ..Default::default()
    };
    let m = random(5, Alphabet::Dna, &skewed, 7).unwrap();
    assert_eq!(m.bkg, vec![0.5, 0.5, 0.0, 0.0]);
    assert_eq!(m.nsites, Some(20.0));
    for c in m.matrix.columns() {
        assert_eq!((c[2], c[3]), (0.0, 0.0));
    }

    let aa = random(3, Alphabet::Aa, &CreateOptions::default(), 3).unwrap();
    assert_eq!(aa.matrix.nrow(), 20);
    assert_eq!(aa.consensus.as_ref().map(|c| c.len()), Some(3));

    let xyz = random(
        3,
        Alphabet::custom("xyz".chars()).unwrap(),
        &CreateOptions::default(),
        3,
    )
    .unwrap();
    assert_eq!(xyz.matrix.nrow(), 3);
    assert_eq!(xyz.consensus, None);
}

#[test]
fn random_motif_errors() {
    let opts = CreateOptions::default();
    assert_eq!(
        arg_msg(random(0, Alphabet::Dna, &opts, 1)),
        "width must be greater than 0"
    );
    let bkg = CreateOptions {
        bkg: Some(vec![0.5, 0.5]),
        ..Default::default()
    };
    assert_eq!(
        arg_msg(random(3, Alphabet::Dna, &bkg, 1)),
        "bkg must have 4 values, one per letter"
    );
    let nsites = CreateOptions {
        nsites: Some(-3.0),
        ..Default::default()
    };
    assert_eq!(
        arg_msg(random(3, Alphabet::Dna, &nsites, 1)),
        "nsites must be a positive integer"
    );
}

// Column conversions

#[test]
fn column_conversions_by_hand() {
    assert_eq!(pcm_to_ppm(&[3., 1., 0., 6.], 0.0), vec![0.3, 0.1, 0.0, 0.6]);
    assert_eq!(
        pcm_to_ppm(&[2., 0., 0., 0.], 4.0),
        vec![0.5, 1. / 6., 1. / 6., 1. / 6.]
    );

    assert_eq!(
        ppm_to_pcm(&[0.3, 0.1, 0.0, 0.6], 10.0),
        vec![3., 1., 0., 6.]
    );
    // An nsites of 1 or less, or NaN, means 100.
    assert_eq!(ppm_to_pcm(&[0.25; 4], 1.0), vec![25.; 4]);
    assert_eq!(ppm_to_pcm(&[0.25; 4], f64::NAN), vec![25.; 4]);
    // Rounding 0.5 up four times overshoots by 2, so the counts are floored
    // and the sites go to the earliest of the tied remainders.
    assert_eq!(ppm_to_pcm(&[0.25; 4], 2.0), vec![1., 1., 0., 0.]);
    // A column summing to more than 1 cannot be floored back to nsites; the
    // excess is taken from the largest count.
    assert_eq!(ppm_to_pcm(&[0.6, 0.6], 10.0), vec![4., 6.]);

    assert_eq!(
        ppm_to_pwm(&[0.5, 0.25, 0.125, 0.125], &[0.25; 4], 0.0, 10.0),
        vec![1., 0., -1., -1.]
    );
    // A background of the wrong length is replaced by a uniform one.
    assert_eq!(
        ppm_to_pwm(&[0.5, 0.5], &[0.1, 0.2, 0.7], 0.0, 10.0),
        vec![0., 0.]
    );
    // With a pseudocount the column goes through 8 counts first.
    let got = ppm_to_pwm(&[0.5, 0.25, 0.125, 0.125], &[0.25; 4], 1.0, 8.0);
    let want: Vec<f64> = [4.25, 2.25, 1.25, 1.25]
        .iter()
        .map(|c| (c / 9.0 / 0.25f64).log2())
        .collect();
    assert_vec_close(&got, &want);

    // Log2 probabilities sum to 1 as they are.
    assert_eq!(
        pwm_to_ppm(&[-2., -2., -2., -2.], &[0.1, 0.2, 0.3, 0.4]),
        vec![0.25; 4]
    );
    assert_eq!(
        pwm_to_ppm(&[1., 0., -1., -1.], &[0.25; 4]),
        vec![0.5, 0.25, 0.125, 0.125]
    );
    assert_vec_close(&pwm_to_ppm(&[1., 0., 0., 0.], &[]), &[0.4, 0.2, 0.2, 0.2]);

    assert_eq!(
        icm_to_ppm(&[0.5, 0.25, 0.25, 0.0]),
        vec![0.5, 0.25, 0.25, 0.0]
    );
    assert_eq!(cwm_to_ppm(&[0.0; 4]), vec![0.25; 4]);
    assert_eq!(
        cwm_to_ppm(&[1.0, -3.0, 0.0, 0.0]),
        vec![0.25, 0.75, 0.0, 0.0]
    );
}

#[test]
fn ppm_to_icm_backgrounds() {
    let p = [0.5, 0.25, 0.125, 0.125];
    // Entropy is 1.75 bits, leaving 0.25 bits to share out.
    assert_eq!(
        ppm_to_icm(&p, &[0.25; 4], false),
        vec![0.125, 0.0625, 0.03125, 0.03125]
    );
    // Relative entropy gives nothing to letters at or below the background.
    assert_eq!(ppm_to_icm(&p, &[0.25; 4], true), vec![0.5, 0.0, 0.0, 0.0]);
    // A background that does not sum to 1 is normalized, and one of the wrong
    // length is replaced by a uniform one.
    assert_eq!(
        ppm_to_icm(&[1., 0., 0., 0.], &[2.; 4], true),
        vec![2., 0., 0., 0.]
    );
    assert_eq!(
        ppm_to_icm(&[1., 0., 0., 0.], &[0.5], true),
        vec![2., 0., 0., 0.]
    );
    assert_eq!(
        ppm_to_icm(&[1., 0., 0., 0.], &[0.5, 0.5, 0.0, 0.0], true),
        vec![1., 0., 0., 0.]
    );
}

#[test]
fn position_icscore_by_type() {
    let bkg = [0.25; 4];
    // ICM columns are simply summed.
    assert_eq!(
        position_icscore(&[0.5, 0.25, 0., 0.], &bkg, MotifType::Icm, 1.0, 10.0, false),
        0.75
    );
    // An nsites below 1 is taken as 100, which leaves 0.25 at 25 counts.
    assert_eq!(
        position_icscore(&[0.25; 4], &bkg, MotifType::Ppm, 0.0, 0.0, false),
        0.0
    );
    assert_eq!(
        position_icscore(&[1., 0., 0., 0.], &bkg, MotifType::Ppm, 0.0, 0.5, false),
        2.0
    );
    assert_eq!(
        position_icscore(&[7., 0., 0., 0.], &bkg, MotifType::Pcm, 0.0, 7.0, false),
        2.0
    );
    assert_eq!(
        position_icscore(&[1., 0., 0., 0.], &bkg, MotifType::Pwm, 0.0, 7.0, false),
        2.0 - shannon_ic(&[0.25; 4]) - (2.0 - shannon_ic(&[0.4, 0.2, 0.2, 0.2]))
    );
    assert_eq!(
        position_icscore(&[-4., 0., 0., 0.], &bkg, MotifType::Cwm, 0.0, 7.0, false),
        2.0
    );
    assert_eq!(
        position_icscore(
            &[1., 0., 0., 0.],
            &[0.5, 0.5, 0., 0.],
            MotifType::Ppm,
            0.0,
            10.0,
            true
        ),
        1.0
    );
}

#[test]
fn consensus_letters_round_trip_through_columns() {
    for l in "ACGTMRWSYKHBVDN".chars() {
        let col = consensus_to_ppm(l).unwrap();
        assert_eq!(consensus_letter(&col, false, MotifType::Ppm, 0.0), l, "{l}");
    }
    assert_eq!(
        consensus_to_ppm('U').unwrap(),
        consensus_to_ppm('T').unwrap()
    );
    assert_eq!(
        consensus_letter(&consensus_to_ppm('U').unwrap(), true, MotifType::Ppm, 0.0),
        'U'
    );
    for l in ['-', '+', '.'] {
        assert_eq!(consensus_to_ppm(l).unwrap(), [0.25; 4]);
    }
    assert_eq!(
        arg_msg(consensus_to_ppm('Z')),
        "'Z' is not a DNA IUPAC symbol"
    );

    for &l in AMINO_ACIDS.iter().chain(&['B', 'Z', 'J', 'X']) {
        let col = consensus_to_ppm_aa(l).unwrap();
        assert_eq!(consensus_letter_aa(&col, MotifType::Ppm, 0.0), l, "{l}");
    }
    for l in ['U', 'O', '*', '-', '+', '.'] {
        assert_eq!(consensus_to_ppm_aa(l).unwrap(), [0.05; 20]);
    }
    assert_eq!(
        arg_msg(consensus_to_ppm_aa('1')),
        "'1' is not an amino acid IUPAC symbol"
    );
}

#[test]
fn consensus_letters_from_other_types() {
    assert_eq!(
        consensus_letter(&[0., 0., 9., 1.], false, MotifType::Pcm, 0.0),
        'G'
    );
    // A large pseudocount flattens the same counts.
    assert_eq!(
        consensus_letter(&[0., 0., 9., 1.], false, MotifType::Pcm, 100.0),
        'N'
    );
    assert_eq!(
        consensus_letter(&[2., 0., 0., 0.], false, MotifType::Pwm, 0.0),
        'A'
    );
    assert_eq!(
        consensus_letter(&[0., 0., 0., 1.5], false, MotifType::Icm, 0.0),
        'T'
    );
    assert_eq!(
        consensus_letter(&[0., -1., 0., 0.], true, MotifType::Cwm, 0.0),
        'C'
    );

    let mut aa = [0.0; 20];
    aa[0] = 0.09;
    aa[1] = 0.08;
    // Nothing reaches 0.1, so the column is unclear.
    assert_eq!(consensus_letter_aa(&aa, MotifType::Ppm, 0.0), 'X');
    // A tie for the top letter is unclear too.
    aa[0] = 0.3;
    aa[1] = 0.3;
    assert_eq!(consensus_letter_aa(&aa, MotifType::Ppm, 0.0), 'X');
    aa[1] = 0.2;
    assert_eq!(consensus_letter_aa(&aa, MotifType::Pcm, 0.0), 'A');
}

#[test]
#[should_panic(expected = "a DNA/RNA column needs four values")]
fn consensus_letter_needs_four_values() {
    consensus_letter(&[1.0, 0.0, 0.0], false, MotifType::Ppm, 0.0);
}

#[test]
#[should_panic(expected = "an amino acid column needs twenty values")]
fn consensus_letter_aa_needs_twenty_values() {
    consensus_letter_aa(&[0.25; 4], MotifType::Ppm, 0.0);
}

// Motif conversions

#[test]
fn convert_type_round_trips_exactly() {
    let ppm = Motif::builder(
        mat(&[&[0.5, 0.25, 0.125, 0.125], &[0.125, 0.125, 0.25, 0.5]]),
        Alphabet::Dna,
    )
    .pseudocount(0.0)
    .nsites(Some(8.0))
    .build()
    .unwrap();
    let pcm = convert_type(&ppm, MotifType::Pcm, None, false).unwrap();
    assert_eq!(pcm.matrix.col(0), &[4., 2., 1., 1.]);
    let pwm = convert_type(&ppm, MotifType::Pwm, None, false).unwrap();
    assert_eq!(pwm.matrix.col(0), &[1., 0., -1., -1.]);
    let icm = convert_type(&ppm, MotifType::Icm, None, false).unwrap();
    assert_eq!(icm.matrix.col(1), &[0.03125, 0.03125, 0.0625, 0.125]);
    // The derived fields are carried over unchanged.
    for m in [&pcm, &pwm, &icm] {
        assert_eq!(m.icscore, ppm.icscore);
        assert_eq!(m.consensus, ppm.consensus);
    }
    assert_eq!(ppm.icscore, 0.5);

    // Every pair of types, there and back.
    let all = [&ppm, &pcm, &pwm, &icm];
    for a in all {
        for b in all {
            let there = convert_type(a, b.motif_type, None, false).unwrap();
            assert_eq!(
                there.matrix, b.matrix,
                "{} to {}",
                a.motif_type, b.motif_type
            );
        }
    }
}

#[test]
fn convert_type_to_its_own_type_is_a_copy() {
    let m = as_pcm();
    assert_eq!(
        convert_type(&m, MotifType::Pcm, Some(5.0), true).unwrap(),
        m
    );
}

#[test]
fn convert_type_refuses_cwm_targets() {
    let msg = unsupported_msg(convert_type(&as_pcm(), MotifType::Cwm, None, false));
    assert!(
        msg.starts_with("CWM cannot be derived from probabilistic motif types (PCM)"),
        "{msg}"
    );
}

#[test]
fn convert_type_with_a_new_pseudocount() {
    let ppm = Motif::builder(mat(&[&[0.5, 0.25, 0.125, 0.125]]), Alphabet::Dna)
        .pseudocount(0.0)
        .nsites(Some(8.0))
        .build()
        .unwrap();
    let pwm = convert_type(&ppm, MotifType::Pwm, Some(1.0), false).unwrap();
    assert_eq!(pwm.pseudocount, 1.0);
    let probs = [4.25 / 9., 2.25 / 9., 1.25 / 9., 1.25 / 9.];
    let want: Vec<f64> = probs.iter().map(|p| (p / 0.25f64).log2()).collect();
    assert_vec_close(pwm.matrix.col(0), &want);
    assert!((pwm.icscore - shannon_ic(&probs)).abs() < 1e-12);
    assert_ne!(pwm.icscore, ppm.icscore);

    // Without nsites a PPM is taken to have 100 sites.
    let no_n = Motif::from_parts(MotifParts {
        nsites: None,
        ..ppm.into_parts()
    })
    .unwrap();
    let pwm = convert_type(&no_n, MotifType::Pwm, Some(1.0), false).unwrap();
    // 12.5 rounds up twice, and the extra site comes off the 50.
    let probs = [49.25 / 101., 25.25 / 101., 13.25 / 101., 13.25 / 101.];
    let want: Vec<f64> = probs.iter().map(|p| (p / 0.25f64).log2()).collect();
    assert_vec_close(pwm.matrix.col(0), &want);
}

#[test]
fn convert_type_softens_zero_backgrounds() {
    let ppm = Motif::builder(mat(&[&[0.5, 0.5, 0.0, 0.0]]), Alphabet::Dna)
        .pseudocount(0.0)
        .bkg(Some(vec![0.5, 0.5, 0.0, 0.0]))
        .build()
        .unwrap();
    let pwm = convert_type(&ppm, MotifType::Pwm, None, false).unwrap();
    // The background is turned into counts of 1000 with a pseudocount of 1.
    let b: [f64; 4] = [500.25 / 1001., 500.25 / 1001., 0.25 / 1001., 0.25 / 1001.];
    let want = [
        (0.5 / b[0]).log2(),
        (0.5 / b[1]).log2(),
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    assert_vec_close(pwm.matrix.col(0), &want);
    // The stored background is left alone.
    assert_eq!(pwm.bkg, vec![0.5, 0.5, 0.0, 0.0]);
}

#[test]
fn convert_type_from_cwm_recomputes_the_total_ic() {
    let cwm = Motif::builder(
        mat(&[&[2., -2., 0., 0.], &[0., 0., 0., -0.5]]),
        Alphabet::Dna,
    )
    .motif_type(Some(MotifType::Cwm))
    .pseudocount(0.0)
    .build()
    .unwrap();
    assert!(cwm.icscore.is_nan());
    let ppm = convert_type(&cwm, MotifType::Ppm, None, false).unwrap();
    assert_eq!(ppm.matrix.as_slice(), &[0.5, 0.5, 0., 0., 0., 0., 0., 1.]);
    assert_eq!(ppm.icscore, 3.0);
    for t in [MotifType::Pcm, MotifType::Pwm, MotifType::Icm] {
        let m = convert_type(&cwm, t, None, false).unwrap();
        assert!(!m.icscore.is_nan(), "{t}");
        let back = convert_type(&m, MotifType::Ppm, None, false).unwrap();
        assert_matrix_close(&back.matrix, &ppm.matrix);
    }
}

#[test]
fn convert_type_relative_entropy_is_not_invertible() {
    let ppm = Motif::builder(mat(&[&[0.5, 0.25, 0.125, 0.125]]), Alphabet::Dna)
        .pseudocount(0.0)
        .build()
        .unwrap();
    let icm = convert_type(&ppm, MotifType::Icm, None, true).unwrap();
    assert_eq!(icm.matrix.col(0), &[0.5, 0., 0., 0.]);
    let back = convert_type(&icm, MotifType::Ppm, None, false).unwrap();
    assert_eq!(back.matrix.col(0), &[1., 0., 0., 0.]);
}

#[test]
fn convert_type_amino_acids_round_trip() {
    let m = from_strings(
        &["MKVLA", "MKILA", "MRVLA"],
        None,
        &CreateOptions {
            motif_type: MotifType::Pcm,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(m.alphabet, Alphabet::Aa);
    assert_eq!(m.consensus.as_deref(), Some("MKVLA"));
    for t in [MotifType::Ppm, MotifType::Pwm, MotifType::Icm] {
        let other = convert_type(&m, t, None, false).unwrap();
        let back = convert_type(&other, MotifType::Pcm, None, false).unwrap();
        assert_eq!(back.matrix, m.matrix, "{t}");
    }
}

// Reverse complements

#[test]
fn reverse_complement_by_hand_and_twice() {
    let m = as_pcm();
    let rc = reverse_complement(&m, false).unwrap();
    assert_eq!(rc.matrix.as_slice(), &[0., 5., 5., 0., 0., 0., 0., 10.]);
    assert_eq!(rc.consensus.as_deref(), Some("ST"));
    assert_eq!(rc.nsites, m.nsites);
    assert!((rc.icscore - m.icscore).abs() < 1e-12);
    assert_eq!(reverse_complement(&rc, false).unwrap(), m);
}

#[test]
fn reverse_complement_moves_gaps() {
    let m = from_strings(
        &["ACGTTA"],
        None,
        &CreateOptions {
            name: "g".into(),
            altname: Some("alt".into()),
            pval: Some(0.1),
            extrainfo: vec![("k".into(), "v".into())],
            ..Default::default()
        },
    )
    .unwrap();
    let g = add_gap(&m, Some(&[1, 4]), &[1, 2], &[3, 5]).unwrap();
    let rc = reverse_complement(&g, false).unwrap();
    assert_eq!(
        rc.gap,
        GapInfo {
            isgapped: true,
            gaploc: vec![2, 5],
            mingap: vec![2, 1],
            maxgap: vec![5, 3],
        }
    );
    assert_eq!(rc.consensus.as_deref(), Some("TAACGT"));
    assert_eq!(rc.altname.as_deref(), Some("alt"));
    assert_eq!(rc.pval, Some(0.1));
    assert_eq!(rc.extrainfo, g.extrainfo);
    assert_eq!(reverse_complement(&rc, false).unwrap(), g);
}

#[test]
fn reverse_complement_of_other_alphabets() {
    let aa = from_strings(&["MFL"], None, &CreateOptions::default()).unwrap();
    assert_eq!(
        unsupported_msg(reverse_complement(&aa, false)),
        "only DNA/RNA motifs can be reverse complemented"
    );
    let xyz = Alphabet::custom("XYZ".chars()).unwrap();
    let m = Motif::builder(mat(&[&[1., 2., 3.], &[6., 0., 0.]]), xyz)
        .build()
        .unwrap();
    let rev = reverse_complement(&m, true).unwrap();
    assert_eq!(rev.matrix.as_slice(), &[0., 0., 6., 3., 2., 1.]);
    assert_eq!(rev.alphabet, m.alphabet);

    let rna = from_strings(&["ACGU"], None, &CreateOptions::default()).unwrap();
    let rc = reverse_complement(&rna, false).unwrap();
    assert_eq!(rc.consensus.as_deref(), Some("ACGU"));
}

// Trimming

#[test]
fn trim_motif_by_ends() {
    let m = trim_ppm();
    assert_eq!(m.icscore, 4.0);
    let both = trim_motif(&m, 0.25, TrimFrom::Both).unwrap().unwrap();
    assert_eq!(both.width(), 3);
    assert_eq!(both.consensus.as_deref(), Some("ANC"));
    assert_eq!(both.icscore, 4.0);
    let left = trim_motif(&m, 0.25, TrimFrom::Left).unwrap().unwrap();
    assert_eq!(left.consensus.as_deref(), Some("ANCNN"));
    let right = trim_motif(&m, 0.25, TrimFrom::Right).unwrap().unwrap();
    assert_eq!(right.consensus.as_deref(), Some("NANC"));
    // A cutoff of 0 keeps everything, and one above 2 bits keeps nothing.
    assert_eq!(trim_motif(&m, 0.0, TrimFrom::Both).unwrap().unwrap(), m);
    assert_eq!(trim_motif(&m, 2.5, TrimFrom::Both).unwrap(), None);
    assert_eq!(trim_motif(&m, 2.5, TrimFrom::Right).unwrap(), None);
}

#[test]
fn trim_motif_shifts_and_drops_gaps() {
    let g = add_gap(&trim_ppm(), Some(&[2, 4]), &[1, 3], &[2, 5]).unwrap();
    let t = trim_motif(&g, 0.25, TrimFrom::Both).unwrap().unwrap();
    assert_eq!(
        t.gap,
        GapInfo {
            isgapped: true,
            gaploc: vec![1],
            mingap: vec![1],
            maxgap: vec![2],
        }
    );
    // A gap left at the very edge is dropped, and with it the gapped flag.
    let g = add_gap(&trim_ppm(), Some(&[1]), &[1], &[2]).unwrap();
    let t = trim_motif(&g, 0.25, TrimFrom::Both).unwrap().unwrap();
    assert_eq!(t.gap, GapInfo::default());
}

fn trim_cwm_motif() -> Motif {
    Motif::builder(
        mat(&[
            &[0.01, 0., 0., 0.],
            &[1., -0.5, 0., 0.],
            &[0., 0., 0., 0.05],
            &[0., 2., 0., 0.],
            &[0.02, 0., 0., 0.],
        ]),
        Alphabet::Dna,
    )
    .motif_type(Some(MotifType::Cwm))
    .build()
    .unwrap()
}

#[test]
fn trim_cwm_keeps_the_middle() {
    let m = trim_cwm_motif();
    // Column sums of |x| are 0.01, 1.5, 0.05, 2 and 0.02, so the cutoff is 0.6.
    let both = trim_cwm(&m, 0.3, None, TrimFrom::Both).unwrap().unwrap();
    assert_eq!(both.matrix, m.matrix.slice_columns(1..4));
    assert!(both.icscore.is_nan());
    assert_eq!(both.consensus.as_deref(), Some("MTC"));
    let left = trim_cwm(&m, 0.3, None, TrimFrom::Left).unwrap().unwrap();
    assert_eq!(left.matrix, m.matrix.slice_columns(1..5));
    let right = trim_cwm(&m, 0.3, None, TrimFrom::Right).unwrap().unwrap();
    assert_eq!(right.matrix, m.matrix.slice_columns(0..4));
    // An absolute cutoff overrides the relative one.
    let abs = trim_cwm(&m, 0.0, Some(1.6), TrimFrom::Both)
        .unwrap()
        .unwrap();
    assert_eq!(abs.matrix, m.matrix.slice_columns(3..4));
    assert_eq!(trim_cwm(&m, 0.3, Some(5.0), TrimFrom::Both).unwrap(), None);
    // The largest column always reaches a relative cutoff of 1.
    let top = trim_cwm(&m, 1.0, None, TrimFrom::Both).unwrap().unwrap();
    assert_eq!(top.width(), 1);
}

#[test]
fn trim_cwm_errors_and_empty_motifs() {
    assert_eq!(
        unsupported_msg(trim_cwm(&as_pcm(), 0.3, None, TrimFrom::Both)),
        "trim_cwm requires CWM motifs; got type PCM"
    );
    let m = trim_cwm_motif();
    for t in [f64::NAN, -0.1] {
        assert_eq!(
            arg_msg(trim_cwm(&m, t, None, TrimFrom::Both)),
            "trim_threshold must be a single non-negative number"
        );
    }
    for a in [f64::NAN, -1.0] {
        assert_eq!(
            arg_msg(trim_cwm(&m, 0.3, Some(a), TrimFrom::Both)),
            "abs_threshold must be a single non-negative number"
        );
    }
    let zero = Motif::builder(mat(&[&[0.; 4], &[0.; 4]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Cwm))
        .build()
        .unwrap();
    assert_eq!(trim_cwm(&zero, 0.3, None, TrimFrom::Both).unwrap(), None);
    // With an absolute cutoff of 0 even an all-zero CWM survives.
    let kept = trim_cwm(&zero, 0.3, Some(0.0), TrimFrom::Both)
        .unwrap()
        .unwrap();
    assert_eq!(kept.width(), 2);
}

// Rounding, normalization, gaps and alphabets

#[test]
fn round_motif_by_hand() {
    let pcm = Motif::builder(mat(&[&[1., 9., 90., 0.]]), Alphabet::Dna)
        .pseudocount(0.0)
        .build()
        .unwrap();
    // 0.01 is dropped and 0.09 and 0.9 are rescaled to 9.09 and 90.91 counts.
    let r = round_motif(&pcm, 0.05).unwrap();
    assert_eq!(r.motif_type, MotifType::Pcm);
    assert_eq!(r.matrix.col(0), &[0., 9., 91., 0.]);
    assert_eq!(r.consensus.as_deref(), Some("G"));

    let ppm = Motif::builder(mat(&[&[0.05, 0.15, 0.3, 0.5]]), Alphabet::Dna)
        .pseudocount(0.0)
        .build()
        .unwrap();
    let r = round_motif(&ppm, 0.1).unwrap();
    assert_vec_close(r.matrix.col(0), &[0.0, 0.15 / 0.95, 0.3 / 0.95, 0.5 / 0.95]);
    // A tolerance of 0 changes nothing.
    assert_matrix_close(&round_motif(&ppm, 0.0).unwrap().matrix, &ppm.matrix);

    for t in [-0.1, 1.5, f64::NAN] {
        assert_eq!(
            arg_msg(round_motif(&ppm, t)),
            "pct_tolerance must be a numeric value between 0 and 1"
        );
    }
}

#[test]
fn normalize_reapplies_the_pseudocount() {
    let ppm = Motif::builder(mat(&[&[1., 0., 0., 0.]]), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .pseudocount(0.0)
        .nsites(Some(10.0))
        .build()
        .unwrap();
    let n = normalize(&ppm).unwrap();
    assert_eq!(n.pseudocount, 1.0);
    assert_vec_close(
        n.matrix.col(0),
        &[10.25 / 11., 0.25 / 11., 0.25 / 11., 0.25 / 11.],
    );
    // Normalizing again keeps the pseudocount and the counts behind it.
    let again = normalize(&n).unwrap();
    assert_eq!(again.pseudocount, 1.0);

    let pcm = as_pcm();
    let n = normalize(&pcm).unwrap();
    assert_eq!(n.matrix, pcm.matrix);
    assert_eq!(n.pseudocount, pcm.pseudocount);
}

#[test]
fn add_gap_options_and_errors() {
    let m = trim_ppm();
    // Without gaploc the gap goes in the middle.
    let mid = add_gap(&m, None, &[1], &[3]).unwrap();
    assert_eq!(
        mid.gap,
        GapInfo {
            isgapped: true,
            gaploc: vec![3],
            mingap: vec![1],
            maxgap: vec![3],
        }
    );
    // Shorter lists are recycled.
    let two = add_gap(&m, Some(&[1, 3]), &[0], &[2]).unwrap();
    assert_eq!(two.gap.mingap, vec![0, 0]);
    assert_eq!(two.gap.maxgap, vec![2, 2]);
    assert!(two.to_string().contains("Consensus:   N..AN..CNN"));

    assert_eq!(
        arg_msg(add_gap(&m, Some(&[]), &[1], &[2])),
        "gaploc, mingap and maxgap must not be empty"
    );
    assert_eq!(
        arg_msg(add_gap(&m, None, &[], &[2])),
        "gaploc, mingap and maxgap must not be empty"
    );
    assert_eq!(
        arg_msg(add_gap(&m, None, &[3], &[2])),
        "mingap cannot be greater than the corresponding maxgap"
    );
    assert_eq!(
        arg_msg(add_gap(&m, None, &[0], &[0])),
        "maxgap must be greater than 0"
    );
    assert_eq!(
        invalid_msgs(add_gap(&m, Some(&[6]), &[1], &[2])),
        vec!["gap location values should not exceed motif size - 1"]
    );
    assert_eq!(
        invalid_msgs(add_gap(&m, Some(&[0]), &[1], &[2])),
        vec!["position 0 gaps or less are not allowed"]
    );
}

#[test]
fn ungap_keeps_or_deletes_gaps() {
    let g = add_gap(&trim_ppm(), Some(&[2]), &[1], &[4]).unwrap();
    let off = ungap(&g, false);
    assert!(!off.gap.isgapped);
    assert_eq!(
        (
            off.gap.gaploc.clone(),
            off.gap.mingap.clone(),
            off.gap.maxgap.clone()
        ),
        (vec![2], vec![1], vec![4])
    );
    assert!(off.to_string().contains("Consensus:   NANCNN"));
    let gone = ungap(&g, true);
    assert_eq!(gone.gap, GapInfo::default());
    assert_eq!(gone, trim_ppm());
}

#[test]
fn switch_alph_swaps_t_and_u() {
    let dna = from_strings(&["TTAGN"], None, &CreateOptions::default()).unwrap();
    let rna = switch_alph(&dna).unwrap();
    assert_eq!(rna.alphabet, Alphabet::Rna);
    assert_eq!(rna.consensus.as_deref(), Some("UUAGN"));
    assert_eq!(rna.matrix, dna.matrix);
    assert_eq!(switch_alph(&rna).unwrap(), dna);

    let aa = from_strings(&["MFL"], None, &CreateOptions::default()).unwrap();
    assert_eq!(
        arg_msg(switch_alph(&aa)),
        "motif 'motif' has alphabet 'AA'; will not convert"
    );
}

// Merging

#[test]
fn merge_motifs_errors_and_single_motifs() {
    let opts = MergeOptions::default();
    assert_eq!(arg_msg(merge_motifs(&[], &opts)), "no motifs given");
    let dna = from_strings(&["ACGTA"], None, &CreateOptions::default()).unwrap();
    let rna = from_strings(&["ACGUA"], None, &CreateOptions::default()).unwrap();
    assert_eq!(
        arg_msg(merge_motifs(&[dna, rna], &opts)),
        "all motifs must share the same alphabet"
    );

    // One motif is just converted to a PPM, and renamed if asked.
    let pcm = as_pcm();
    let one = merge_motifs_detailed(
        std::slice::from_ref(&pcm),
        &MergeOptions {
            name: Some("solo".into()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(one.used, vec![0]);
    assert_eq!(one.motif.name, "solo");
    assert_eq!(
        one.motif.matrix,
        convert_type(&pcm, MotifType::Ppm, None, false)
            .unwrap()
            .matrix
    );
    let unnamed = merge_motifs(std::slice::from_ref(&pcm), &opts).unwrap();
    assert_eq!(unnamed.name, pcm.name);
}

#[test]
fn merge_motifs_averages_aligned_columns() {
    let a = sharp_ppm("ACGTA", 0.9, "a", 30.0);
    let b = sharp_ppm("ACGTA", 0.6, "b", 10.0);
    let merged = merge_motifs_detailed(&[b.clone(), a.clone()], &MergeOptions::default()).unwrap();
    assert_eq!(merged.used, vec![0, 1]);
    let m = merged.motif;
    // The names follow input order, even though a is the anchor.
    assert_eq!(m.name, "merged_b+a");
    assert_eq!(m.nsites, Some(40.0));
    assert_eq!(m.motif_type, MotifType::Ppm);
    assert_eq!(m.consensus.as_deref(), Some("ACGTA"));
    assert_matrix_close(&m.matrix, &sharp_ppm("ACGTA", 0.75, "x", 1.0).matrix);

    let weighted = merge_motifs(
        &[a, b],
        &MergeOptions {
            weighted: true,
            name: Some("w".into()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(weighted.name, "w");
    // (0.9 * 30 + 0.6 * 10) / 40
    assert_matrix_close(
        &weighted.matrix,
        &sharp_ppm("ACGTA", 0.825, "x", 1.0).matrix,
    );
}

#[test]
fn merge_motifs_with_reverse_complements() {
    let a = sharp_ppm("AACGTTG", 0.9, "a", 10.0);
    let rc = reverse_complement(&a, false).unwrap();
    let merged = merge_motifs(&[a.clone(), rc.clone()], &MergeOptions::default()).unwrap();
    assert_matrix_close(&merged.matrix, &a.matrix);
    assert_eq!(merged.nsites, Some(20.0));

    // Without reverse complements the second motif cannot line up in full.
    let fwd = merge_motifs(
        &[a.clone(), rc],
        &MergeOptions {
            rc: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_ne!(fwd.consensus, a.consensus);
    assert_eq!(fwd.width(), 8, "{fwd}");
}

#[test]
fn merge_motifs_extends_partial_overlaps_and_skips_short_ones() {
    let a = sharp_ppm("ACGTAC", 0.97, "a", 10.0);
    let b = sharp_ppm("GTACGG", 0.9, "b", 10.0);
    let merged = merge_motifs_detailed(
        &[a.clone(), b],
        &MergeOptions {
            min_overlap: 4,
            rc: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(merged.used, vec![0, 1]);
    assert_eq!(merged.motif.width(), 8);
    assert_eq!(merged.motif.consensus.as_deref(), Some("ACGTACGG"));
    // The columns covered by one motif only are that motif's.
    assert_vec_close(merged.motif.matrix.col(0), a.matrix.col(0));
    assert_vec_close(
        merged.motif.matrix.col(7),
        &[0.1 / 3., 0.1 / 3., 0.9, 0.1 / 3.],
    );

    let short = sharp_ppm("GTA", 0.9, "c", 10.0);
    let merged = merge_motifs_detailed(
        &[a.clone(), short],
        &MergeOptions {
            min_overlap: 5,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(merged.used, vec![0]);
    assert_eq!(merged.motif.name, "merged_a");
    assert_matrix_close(&merged.motif.matrix, &a.matrix);
    assert_eq!(merged.motif.nsites, Some(10.0));
}

#[test]
fn merge_motifs_counts_missing_nsites_as_one() {
    let pwm = |name: &str| {
        Motif::builder(
            mat(&[
                &[1., -1., -1., -1.],
                &[-1., 1., -1., -1.],
                &[-1., -1., 1., -1.],
                &[-1., -1., -1., 1.],
                &[1., -1., -1., -1.],
            ]),
            Alphabet::Dna,
        )
        .name(name)
        .build()
        .unwrap()
    };
    let (a, b) = (pwm("a"), pwm("b"));
    assert_eq!(a.nsites, None);
    let m = merge_motifs(&[a, b], &MergeOptions::default()).unwrap();
    assert_eq!(m.nsites, Some(2.0));
    // 2 against 0.5 three times, out of 3.5.
    assert_vec_close(
        m.matrix.col(0),
        &[2. / 3.5, 0.5 / 3.5, 0.5 / 3.5, 0.5 / 3.5],
    );
}

#[test]
fn merge_similar_edge_cases() {
    let opts = MergeSimilarOptions::default();
    let none = merge_similar(&[], &opts).unwrap();
    assert!(none.clusters.is_empty() && none.motifs.is_empty());

    let a = sharp_ppm("ACGTAC", 0.9, "a", 10.0);
    let one = merge_similar(std::slice::from_ref(&a), &opts).unwrap();
    assert_eq!(one.clusters, vec![0]);
    assert_eq!(one.motifs, vec![a.clone()]);

    let rna = from_strings(&["ACGUAC"], None, &CreateOptions::default()).unwrap();
    assert_eq!(
        arg_msg(merge_similar(&[a.clone(), rna], &opts)),
        "all motifs must share the same alphabet"
    );
    for q in [0.0, 1.5, f64::NAN] {
        let bad = MergeSimilarOptions {
            qvalue: q,
            ..Default::default()
        };
        assert_eq!(
            arg_msg(merge_similar(std::slice::from_ref(&a), &bad)),
            "qvalue must be in (0, 1]"
        );
    }
}
