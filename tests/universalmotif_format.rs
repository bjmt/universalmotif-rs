//! The R package's universalmotif text format, read from files the R package
//! wrote and written back.

use universalmotif::alphabet::Alphabet;
use universalmotif::io::universalmotif::{read, write};
use universalmotif::motif::{MotifType, Strand};

const MIXED: &str = include_str!("reference/files/universalmotif_mixed.yml");
const MINIMAL: &str = include_str!("reference/files/universalmotif_minimal.yml");

#[test]
fn reads_every_field_the_r_package_writes() {
    let (motifs, notes) = read(MIXED).unwrap();
    assert_eq!(motifs.len(), 5);

    // Quoted names, a quoted alternative name that reads as a boolean in
    // YAML, extra information, and matrix rows wrapped over three lines.
    let aa = &motifs[0];
    assert_eq!(aa.name, "toy: protein #1");
    assert_eq!(aa.altname.as_deref(), Some("yes"));
    assert_eq!(aa.family.as_deref(), Some("it's a family"));
    assert_eq!(aa.organism.as_deref(), Some("Homo sapiens"));
    assert_eq!(aa.alphabet, Alphabet::Aa);
    assert_eq!(aa.nsites, Some(3.0));
    assert_eq!(aa.pseudocount, 1.0);
    assert_eq!(
        aa.extrainfo,
        vec![
            ("source".into(), "made up".into()),
            ("note".into(), "a: b".into())
        ]
    );
    let l = Alphabet::Aa.index_of('L').unwrap();
    assert!((aa.matrix.get(l, 0) - 0.7625).abs() < 1e-12);
    assert_eq!(aa.consensus.as_deref(), Some("LDELKR"));

    let tbp = &motifs[1];
    assert_eq!(tbp.name, "MA0108.2 TBP");
    assert_eq!(
        (tbp.nsites, tbp.pval, tbp.qval, tbp.eval),
        (Some(50.0), Some(1e-5), Some(0.01), Some(3.0))
    );
    assert_eq!(tbp.bkg, vec![0.3, 0.2, 0.2, 0.3]);
    assert_eq!(tbp.consensus.as_deref(), Some("TATAAA"));

    let ss = &motifs[2];
    assert_eq!(ss.alphabet, Alphabet::custom("CEH".chars()).unwrap());
    assert_eq!(ss.width(), 8);

    let rna = &motifs[3];
    assert_eq!(
        (rna.alphabet.clone(), rna.motif_type),
        (Alphabet::Rna, MotifType::Pcm)
    );
    assert_eq!(rna.matrix.get(3, 0), 10.0);

    let dr = &motifs[4];
    assert!(dr.gap.isgapped);
    assert_eq!(
        (
            dr.gap.gaploc.clone(),
            dr.gap.mingap.clone(),
            dr.gap.maxgap.clone()
        ),
        (vec![3, 6], vec![0, 1], vec![1, 5])
    );

    // The higher-order matrices of the TBP motif are left out with a note.
    assert_eq!(notes.len(), 1);
    assert!(notes[0].contains("MA0108.2 TBP") && notes[0].contains("multifreq"));
}

#[test]
fn reads_the_minimal_output() {
    let (motifs, notes) = read(MINIMAL).unwrap();
    assert_eq!(motifs.len(), 1);
    assert_eq!(motifs[0].name, "MA0108.2 TBP");
    assert_eq!(motifs[0].strand, Strand::Both);
    assert_eq!(notes.len(), 1);
}

#[test]
fn round_trips_exactly() {
    let (motifs, _) = read(MIXED).unwrap();
    let text = write(&motifs, true).unwrap();
    assert!(text.starts_with("# universalmotif version 1.31.51\n"));
    let (back, notes) = read(&text).unwrap();
    assert!(notes.is_empty());
    assert_eq!(back, motifs);
    // Adding to a file leaves out the header.
    assert!(write(&motifs[..1], false).unwrap().starts_with("\n---\n"));
}

#[test]
fn writes_what_yaml_would_misread_in_quotes() {
    let (mut motifs, _) = read(MIXED).unwrap();
    let mut parts = motifs.remove(1).into_parts();
    parts.name = "-".into();
    parts.altname = Some("no".into());
    parts.strand = Strand::Minus;
    parts.extrainfo = vec![("y".into(), "two\nlines".into()), ("".into(), "1.5".into())];
    let m = universalmotif::motif::Motif::from_parts(parts).unwrap();
    let text = write(std::slice::from_ref(&m), true).unwrap();
    assert!(text.contains("name: '-'\n") && text.contains("strand: '-'\n"));
    let (back, _) = read(&text).unwrap();
    assert_eq!(back[0], m);
}

#[test]
fn refuses_files_from_old_versions() {
    let old = "# universalmotif version 1.1.50\n\nname: x\n";
    let e = read(old).unwrap_err().to_string();
    assert!(e.contains("before 1.1.67"), "{e}");
    assert!(read("").is_err());
}

#[test]
fn quotes_what_yaml_1_1_reads_as_numbers_or_keys() {
    let (motifs, _) = read(MIXED).unwrap();
    for name in [
        "0x1F",
        "0o17",
        "1_000",
        "1:30",
        "012",
        "<<",
        "=",
        "a\u{2028}b",
    ] {
        let mut parts = motifs[0].clone().into_parts();
        parts.name = name.into();
        let m = universalmotif::motif::Motif::from_parts(parts).unwrap();
        let text = write(std::slice::from_ref(&m), true).unwrap();
        let line = text.lines().find(|l| l.starts_with("name:")).unwrap();
        assert!(line.contains('\'') || line.contains('"'), "{line}");
        assert_eq!(read(&text).unwrap().0[0].name, name);
    }
}

#[test]
fn refuses_gap_sizes_too_large() {
    let (motifs, _) = read(MIXED).unwrap();
    let gapped = motifs
        .iter()
        .find(|m| m.gap.isgapped)
        .expect("a gapped motif");
    let text = write(std::slice::from_ref(gapped), true).unwrap();
    assert!(text.contains("maxgap:\n- 1\n- 5\n"), "{text}");
    let bad = text.replace("maxgap:\n- 1\n- 5\n", "maxgap:\n- 1\n- 1.0e30\n");
    let e = read(&bad).unwrap_err().to_string();
    assert!(e.contains("below 2^32"), "{e}");
}
