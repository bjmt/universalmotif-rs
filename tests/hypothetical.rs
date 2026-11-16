//! Matching motifs against hypothetical contribution scores.

use universalmotif::Alphabet;
use universalmotif::contributions::{
    ContribScanOptions, ContribThreshold, scan_contributions, scan_hypothetical,
};
use universalmotif::create::{CreateOptions, from_strings, from_typed_matrix};
use universalmotif::matrix::Matrix;
use universalmotif::motif::{MotifType, Strand};
use universalmotif::ops::add_gap;

const ACGT: [u8; 4] = *b"ACGT";

/// Hypothetical scores for `seq` in which the model wants `wanted` (a letter
/// per position, or b'.' for no preference) with weight `w`, and the
/// projected scores they imply, both in the layout (4, length).
fn scores(seq: &[u8], wanted: &[u8], w: f32) -> (Vec<f32>, Vec<f32>) {
    let len = seq.len();
    let mut hyp = vec![0.0f32; 4 * len];
    for (p, &c) in wanted.iter().enumerate() {
        for (a, &l) in ACGT.iter().enumerate() {
            hyp[a * len + p] = if c == l {
                w
            } else if c == b'.' {
                0.01
            } else {
                -w / 3.0
            };
        }
    }
    let mut proj = vec![0.0f32; 4 * len];
    for (p, &c) in seq.iter().enumerate() {
        let a = ACGT.iter().position(|&l| l == c).unwrap();
        proj[a * len + p] = hyp[a * len + p];
    }
    (hyp, proj)
}

fn ppm(s: &str) -> universalmotif::Motif {
    from_strings(&[s], None, &CreateOptions::default()).unwrap()
}

#[test]
fn a_site_with_a_wrong_letter_still_matches_what_the_model_wants() {
    // TGACTCA with its middle C changed to A; the model still wants C there,
    // and gives the site weight through its other letters.
    let seq = b"AAAATGAATCAAAA";
    let wanted = b"....TGACTCA...";
    let (hyp, proj) = scores(seq, wanted, 0.5);
    let len = seq.len();
    let loose = ContribThreshold {
        min_similarity: 0.5,
        min_importance: 0.1,
    };
    let opts = ContribScanOptions::default();
    let hits =
        scan_hypothetical(&[ppm("TGACTCA")], &hyp, &proj, &[len], len, &[loose], &opts).unwrap();
    let best = hits
        .iter()
        .max_by(|a, b| a.score.total_cmp(&b.score))
        .unwrap();
    assert_eq!((best.start, best.end), (4, 11));
    assert!(best.similarity > 0.9, "{}", best.similarity);
    // The observed scores see a hole at the wrong letter.
    let as_cwm = CreateOptions {
        motif_type: MotifType::Cwm,
        ..Default::default()
    };
    let cwm = from_typed_matrix(
        ppm("TGACTCA").matrix.clone(),
        None,
        Some(Alphabet::Dna),
        MotifType::Cwm,
        &as_cwm,
    )
    .unwrap();
    let observed = scan_contributions(&[cwm], &proj, &[len], len, &[loose], &opts).unwrap();
    let best_observed = observed.iter().map(|h| h.similarity).fold(0.0, f64::max);
    assert!(best_observed < best.similarity);
}

#[test]
fn the_reverse_strand_and_thresholds() {
    // The model wants TGATAA's reverse complement, TTATCA.
    let seq = b"CCCTTATCACCC";
    let wanted = b"...TTATCA...";
    let (hyp, proj) = scores(seq, wanted, 0.5);
    let len = seq.len();
    let m = ppm("TGATAA");
    let t = ContribThreshold {
        min_similarity: 0.5,
        min_importance: 0.1,
    };
    let hits = scan_hypothetical(
        std::slice::from_ref(&m),
        &hyp,
        &proj,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    )
    .unwrap();
    let best = hits
        .iter()
        .max_by(|a, b| a.score.total_cmp(&b.score))
        .unwrap();
    assert_eq!((best.start, best.strand), (3, Strand::Minus));
    // For a PPM, min_importance is a mean per position: the site's is 0.5.
    let high = ContribThreshold {
        min_similarity: 0.5,
        min_importance: 0.6,
    };
    assert!(
        scan_hypothetical(
            std::slice::from_ref(&m),
            &hyp,
            &proj,
            &[len],
            len,
            &[high],
            &ContribScanOptions::default()
        )
        .unwrap()
        .is_empty()
    );
    // Without rc, only the forward strand.
    let fwd = ContribScanOptions {
        rc: false,
        ..Default::default()
    };
    assert!(
        scan_hypothetical(&[m], &hyp, &proj, &[len], len, &[t], &fwd)
            .unwrap()
            .iter()
            .all(|h| h.strand == Strand::Plus)
    );
}

#[test]
fn cwms_are_matched_relative_to_their_favoured_letters() {
    let seq = b"AAGATAAGAA";
    let wanted = b"..GATAAG..";
    let (hyp, proj) = scores(seq, wanted, 0.4);
    let len = seq.len();
    // A hypothetical CWM of GATAAG: 0.4 for the letter wanted, -0.13 for the others.
    let cols: Vec<[f64; 4]> = b"GATAAG"
        .iter()
        .map(|&c| ACGT.map(|l| if l == c { 0.4 } else { -0.4 / 3.0 }))
        .collect();
    let opts = CreateOptions {
        name: "GATA".into(),
        motif_type: MotifType::Cwm,
        ..Default::default()
    };
    let cwm = from_typed_matrix(
        Matrix::from_columns(&cols).unwrap(),
        None,
        Some(Alphabet::Dna),
        MotifType::Cwm,
        &opts,
    )
    .unwrap();
    // The site's importance, 6 x 0.4, is the CWM's own: a fraction of 1 is
    // just reached, and more is not.
    let reach = ContribThreshold {
        min_similarity: 0.9,
        min_importance: 0.99,
    };
    let miss = ContribThreshold {
        min_similarity: 0.9,
        min_importance: 1.01,
    };
    let o = ContribScanOptions::default();
    let hits = scan_hypothetical(
        std::slice::from_ref(&cwm),
        &hyp,
        &proj,
        &[len],
        len,
        &[reach],
        &o,
    )
    .unwrap();
    assert_eq!(hits.len(), 1);
    assert!((hits[0].similarity - 1.0).abs() < 1e-5);
    assert!(
        scan_hypothetical(&[cwm], &hyp, &proj, &[len], len, &[miss], &o)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn gapped_motifs_and_bad_input() {
    // TGA, a 3-letter spacer, then TCA; only the parts carry preferences.
    let seq = b"CCTGACCCTCACC";
    let wanted = b"..TGA...TCA..";
    let (hyp, proj) = scores(seq, wanted, 0.5);
    let len = seq.len();
    let m = add_gap(&ppm("TGATCA"), Some(&[3]), &[2], &[4]).unwrap();
    let t = ContribThreshold {
        min_similarity: 0.5,
        min_importance: 0.1,
    };
    let hits = scan_hypothetical(
        &[m],
        &hyp,
        &proj,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    )
    .unwrap();
    let best = hits
        .iter()
        .max_by(|a, b| a.score.total_cmp(&b.score))
        .unwrap();
    assert_eq!((best.start, best.end, best.gaps.clone()), (2, 11, vec![3]));
    let short = vec![0.0f32; 4 * len - 1];
    let err = scan_hypothetical(
        &[ppm("TGA")],
        &short,
        &proj,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    );
    assert!(
        err.unwrap_err()
            .to_string()
            .contains("hypothetical scores do not fill")
    );
    let protein = from_strings(&["MKV"], Some(Alphabet::Aa), &CreateOptions::default()).unwrap();
    let err = scan_hypothetical(
        &[protein],
        &hyp,
        &proj,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    );
    assert!(
        err.unwrap_err()
            .to_string()
            .contains("not a DNA or RNA motif")
    );
    // Observed scores still need CWMs.
    let err = scan_contributions(
        &[ppm("TGA")],
        &proj,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    );
    assert!(
        err.unwrap_err()
            .to_string()
            .contains("are matched with CWMs")
    );
}
