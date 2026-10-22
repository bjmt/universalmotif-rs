//! Reader cases that universalmotif gets wrong (see DIVERGENCES.md).

use universalmotif::io::{cisbp, homer, jaspar, matrix, meme, transfac};
use universalmotif::{Alphabet, MotifType};

#[test]
fn jaspar_rows_in_any_order() {
    let sorted = jaspar::read(">m\nA [ 1 2 ]\nC [ 3 4 ]\nG [ 5 6 ]\nT [ 7 8 ]\n").unwrap();
    let shuffled = jaspar::read(">m\nT [ 7 8 ]\nG [ 5 6 ]\nA [ 1 2 ]\nC [ 3 4 ]\n").unwrap();
    assert_eq!(sorted, shuffled);
}

#[test]
fn jaspar_custom_and_unlabelled_protein() {
    let m = jaspar::read(">m\nX [ 1 2 ]\nY [ 3 0 ]\nZ [ 0 2 ]\n").unwrap();
    assert_eq!(m[0].alphabet, Alphabet::custom("XYZ".chars()).unwrap());
    let rows: String = (0..20)
        .map(|i| format!("{} {}\n", i % 3, (i + 1) % 2))
        .collect();
    let m = jaspar::read(&format!(">p\n{rows}")).unwrap();
    assert_eq!(m[0].alphabet, Alphabet::Aa);
}

#[test]
fn transfac_keeps_whole_field_values() {
    let text = "ID  m1\nOS  Homo sapiens\nHC  C2H2 zinc finger\nP0 A C G T\n\
                01 1 1 1 1 N\n02 4 0 0 0 A\nXX\n//\n";
    let m = transfac::read(text).unwrap();
    assert_eq!(m[0].organism.as_deref(), Some("Homo sapiens"));
    assert_eq!(m[0].family.as_deref(), Some("C2H2 zinc finger"));
}

#[test]
fn transfac_name_does_not_depend_on_field_order() {
    let a = transfac::read("AC M1\nID first\nP0 A C G T\n01 1 0 0 0 A\n//\n").unwrap();
    let b = transfac::read("ID first\nAC M1\nP0 A C G T\n01 1 0 0 0 A\n//\n").unwrap();
    assert_eq!(a[0].name, "M1");
    assert_eq!(b[0].name, "M1");
    assert_eq!(b[0].altname.as_deref(), Some("first"));
}

#[test]
fn transfac_column_order_from_header() {
    let a = transfac::read("ID m\nP0 A C G T\n01 1 2 3 4 N\n//\n").unwrap();
    let b = transfac::read("ID m\nP0 T G C A\n01 4 3 2 1 N\n//\n").unwrap();
    assert_eq!(a[0].matrix, b[0].matrix);
}

#[test]
fn cisbp_keeps_whole_field_values() {
    let text = "TF\tT1\nTF Name\tZNF 123\nFamily\tC2H2 ZF\nSpecies\tHomo sapiens\n\
                Pos\tA\tC\tG\tT\n1\t0.1\t0.2\t0.3\t0.4\n";
    let m = cisbp::read(text).unwrap();
    assert_eq!(m[0].name, "T1");
    assert_eq!(m[0].altname.as_deref(), Some("ZNF 123"));
    assert_eq!(m[0].family.as_deref(), Some("C2H2 ZF"));
    assert_eq!(m[0].organism.as_deref(), Some("Homo sapiens"));
}

#[test]
fn homer_names_with_spaces() {
    let text = ">ACGT\tmy motif/extra\t5.2\n0.97\t0.01\t0.01\t0.01\n";
    let m = homer::read(text).unwrap();
    assert_eq!(m[0].name, "my motif");
}

#[test]
fn meme_blocks_with_large_positions() {
    let text = "MEME version 5\nALPHABET= ACGT\nstrands: + -\n\n\
        MOTIF A MEME-1\n\tMotif A in BLOCKS format\nBL   MOTIF A width=4 seqs=2\n\
        chr1                     (   12) ACGT  1\nchr2                     (12345) ACGA  1\n//\n\
        letter-probability matrix: alength= 4 w= 4 nsites= 2 E= 1e-3\n\
        1 0 0 0\n0 1 0 0\n0 0 1 0\n0.5 0 0 0.5\n";
    let sites = meme::read_sites(text).unwrap();
    let s: Vec<&str> = sites[0].sites.iter().map(|s| s.site.as_str()).collect();
    assert_eq!(s, ["ACGT", "ACGA"]);
    let m = meme::read(text, &Default::default()).unwrap();
    assert_eq!(m[0].motif_type, MotifType::Ppm);
    assert_eq!(m[0].altname.as_deref(), Some("MEME-1"));
}

#[test]
fn meme_custom_alphabet_background_by_letter() {
    let text = "MEME version 5\nALPHABET \"x\"\nA\nC\n1\nEND ALPHABET\n\
        Background letter frequencies\nA 0.5 C 0.3 1 0.2\n\
        MOTIF m\nletter-probability matrix: alength= 3 w= 1\n0.2 0.3 0.5\n";
    let m = meme::read(text, &Default::default()).unwrap();
    // Rows are stored in character order (1, A, C) while MEME lists letters
    // as A, C, 1.
    assert_eq!(m[0].alphabet.letters(), ['1', 'A', 'C']);
    assert_eq!(m[0].bkg, vec![0.2, 0.5, 0.3]);
    assert_eq!(m[0].matrix.col(0), [0.5, 0.2, 0.3]);
}

#[test]
fn matrix_comments_and_prefix_headers() {
    let text = "# comment\n>one\n0.1 0.2\n0.2 0.2\n0.3 0.3\n0.4 0.3\n>two\n1 0\n0 1\n0 0\n0 0\n";
    let opts = matrix::ReadOptions {
        headers: matrix::Headers::Prefix(">".into()),
        comment: Some("#".into()),
        sep: None,
        ..Default::default()
    };
    let m = matrix::read(text, &opts).unwrap();
    assert_eq!(m.len(), 2);
    assert_eq!(m[1].name, "two");
    assert_eq!(m[1].motif_type, MotifType::Pcm);
}

/// Writers universalmotif does not have: check that what they write reads
/// back unchanged, apart from the conversion to probabilities.
#[test]
fn cisbp_and_uniprobe_writers_round_trip() {
    use universalmotif::convert::convert_type;
    use universalmotif::io::uniprobe;
    let mut motifs = jaspar::read(include_str!("reference/files/jaspar.txt")).unwrap();
    motifs.extend(transfac::read(include_str!("reference/files/transfac_meta.txt")).unwrap());
    let ppms: Vec<_> = motifs
        .iter()
        .map(|m| convert_type(m, MotifType::Ppm, None, false).unwrap())
        .collect();

    for (back, m) in cisbp::read(&cisbp::write(&motifs).unwrap())
        .unwrap()
        .iter()
        .zip(&ppms)
    {
        assert_eq!(back.name, m.name);
        assert_eq!(back.altname, m.altname);
        assert_eq!(back.family, m.family);
        assert_eq!(back.organism, m.organism);
        assert_eq!(back.matrix, m.matrix);
    }
    for (back, m) in uniprobe::read(&uniprobe::write(&motifs).unwrap())
        .unwrap()
        .iter()
        .zip(&ppms)
    {
        assert_eq!(back.name, m.name);
        assert_eq!(back.matrix, m.matrix);
    }

    let rna = cisbp::read(include_str!("reference/files/cisbp_rna.txt")).unwrap();
    assert_eq!(
        cisbp::read(&cisbp::write(&rna).unwrap()).unwrap()[0].alphabet,
        Alphabet::Rna
    );
    assert!(uniprobe::write(&rna).is_err());
}
