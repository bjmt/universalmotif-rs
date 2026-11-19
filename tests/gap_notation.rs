//! Gaps written into the alternative name of MEME, JASPAR and TRANSFAC files.

use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::io::{gap_notation, jaspar, meme, parse_gap_notation, transfac};
use universalmotif::motif::{GapInfo, Motif};
use universalmotif::ops::add_gap;

fn dr() -> Motif {
    let m = from_strings(
        &["AGGTCAAGGTCA"],
        None,
        &CreateOptions {
            name: "DR".into(),
            ..Default::default()
        },
    )
    .unwrap();
    add_gap(&m, Some(&[3, 6]), &[0, 1], &[1, 5]).unwrap()
}

#[test]
fn notation_round_trips() {
    let m = dr();
    let n = gap_notation(&m).unwrap();
    assert_eq!(n, "AGG[0-1]TCA[1-5]AGGTCA");
    assert_eq!(
        parse_gap_notation(&n, 12),
        Some(GapInfo {
            isgapped: true,
            gaploc: vec![3, 6],
            mingap: vec![0, 1],
            maxgap: vec![1, 5],
        })
    );
    // Only the exact notation, with as many letters as columns, counts.
    for not in [
        "AGG[0-1]TCA[1-5]AGGTC",
        "AGGTCA",
        "[1-5]AGGTCA",
        "AGG[1-]TCA",
        "AGG[1-5]TCA ",
        "MA0001.1",
    ] {
        assert_eq!(parse_gap_notation(not, 12), None, "{not}");
    }
    assert_eq!(gap_notation(&universalmotif::ops::ungap(&m, false)), None);
}

#[test]
fn formats_keep_gaps_in_the_alternative_name() {
    let m = dr();
    let texts = [
        meme::write(std::slice::from_ref(&m), &Default::default()).unwrap(),
        jaspar::write(std::slice::from_ref(&m)).unwrap(),
        transfac::write(std::slice::from_ref(&m), &Default::default()).unwrap(),
    ];
    assert!(texts[0].contains("MOTIF DR AGG[0-1]TCA[1-5]AGGTCA\n"));
    assert!(texts[1].starts_with(">DR AGG[0-1]TCA[1-5]AGGTCA\n"));
    assert!(texts[2].contains("NA AGG[0-1]TCA[1-5]AGGTCA\n"));
    let back = [
        meme::read(&texts[0], &Default::default()).unwrap(),
        jaspar::read(&texts[1]).unwrap(),
        transfac::read(&texts[2]).unwrap(),
    ];
    for b in back {
        assert_eq!(b[0].gap, m.gap);
        assert_eq!(b[0].altname.as_deref(), Some("AGG[0-1]TCA[1-5]AGGTCA"));
    }
}

#[test]
fn ordinary_alternative_names_are_left_alone() {
    let text = ">MA0001.1 AGL3\nA [ 1 2 3 ]\nC [ 1 1 1 ]\nG [ 1 1 1 ]\nT [ 1 0 0 ]\n";
    let m = &jaspar::read(text).unwrap()[0];
    assert!(!m.gap.isgapped);
    assert_eq!(m.altname.as_deref(), Some("AGL3"));
    // A notation that does not fit the motif's width is only a name.
    let text = ">x AC[1-2]G\nA [ 1 2 3 4 ]\nC [ 1 1 1 1 ]\nG [ 1 1 1 1 ]\nT [ 1 0 0 0 ]\n";
    assert!(!jaspar::read(text).unwrap()[0].gap.isgapped);
}
