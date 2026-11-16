//! Comparing, laying out and merging gapped motifs with their gaps.

use universalmotif::Alphabet;
use universalmotif::compare::{CompareOptions, anchor_layout, compare_motifs};
use universalmotif::create::{CreateOptions, from_strings, from_typed_matrix};
use universalmotif::matrix::Matrix;
use universalmotif::merge::{
    MergeOptions, MergeSimilarOptions, merge_motifs_detailed, merge_similar,
};
use universalmotif::motif::{Motif, MotifType};
use universalmotif::ops::{add_gap, reverse_complement};
use universalmotif::similarity::{CwmCompareOptions, compare_cwms};

fn ppm(s: &str, name: &str) -> Motif {
    let opts = CreateOptions {
        name: name.into(),
        nsites: Some(20.0),
        ..Default::default()
    };
    from_strings(&[s], None, &opts).unwrap()
}

/// A direct repeat of AGGTCA with a spacer of `lo` to `hi` letters.
fn dr(lo: usize, hi: usize, name: &str) -> Motif {
    add_gap(&ppm("AGGTCAAGGTCA", name), Some(&[6]), &[lo], &[hi]).unwrap()
}

#[test]
fn a_gapped_motif_aligns_to_a_fixed_spacing_at_the_right_size() {
    // The repeat with a spacer of exactly 4, as an ungapped motif.
    let fixed = from_strings(
        &["AGGTCANNNNAGGTCA"],
        None,
        &CreateOptions {
            name: "DR4".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let motifs = [dr(1, 6, "DR"), fixed, ppm("CACGTGAC", "E-box")];
    let with = compare_motifs(&motifs, &[0], &CompareOptions::default()).unwrap();
    assert_eq!(with[1].query_gaps, vec![4]);
    assert!(with[1].target_gaps.is_empty());
    assert_eq!(with[1].overlap, 16);
    assert!(
        with[1].score > 0.99 && with[1].pvalue < 1e-4,
        "{:?}",
        with[1]
    );
    // Run together, only one half can line up with the fixed spacing, which
    // is far weaker evidence even after the six sizes tried are counted.
    let without = compare_motifs(
        &motifs,
        &[0],
        &CompareOptions {
            use_gaps: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(without[1].query_gaps.is_empty());
    assert_eq!(without[1].overlap, 6);
    assert!(with[1].pvalue < without[1].pvalue / 100.0);
    // Unrelated motifs stay unrelated, and each motif matches itself.
    assert!(with[2].qvalue > 0.05);
    assert!(with[0].score > 0.999);
}

#[test]
fn gapped_motifs_against_each_other() {
    // Spacers of 1-3 and 3-5 meet at 3.
    let motifs = [dr(1, 3, "a"), dr(3, 5, "b")];
    let rows = compare_motifs(&motifs, &[0], &CompareOptions::default()).unwrap();
    assert_eq!(
        (rows[1].query_gaps.clone(), rows[1].target_gaps.clone()),
        (vec![3], vec![3])
    );
    assert!(rows[1].score > 0.999);
}

#[test]
fn layouts_place_gapped_motifs_at_their_best_size() {
    let fixed = from_strings(&["AGGTCANNNNNAGGTCA"], None, &CreateOptions::default()).unwrap();
    let layout = anchor_layout(&[fixed, dr(1, 6, "DR")], 0, &CompareOptions::default()).unwrap();
    assert_eq!(layout.gaps, vec![vec![], vec![5]]);
    assert_eq!(layout.width, 17);
    assert_eq!(layout.placements[1].left, 0);
}

#[test]
fn merging_keeps_gaps_spanning_all_of_them() {
    // The second, turned around, still merges with the first part by part.
    let b = reverse_complement(&dr(3, 5, "b"), false).unwrap();
    let out = merge_motifs_detailed(&[dr(1, 3, "a"), b], &MergeOptions::default()).unwrap();
    let m = &out.motif;
    assert!(m.gap.isgapped && out.warnings.is_empty());
    assert_eq!(
        (
            m.gap.gaploc.clone(),
            m.gap.mingap.clone(),
            m.gap.maxgap.clone()
        ),
        (vec![6], vec![1], vec![5])
    );
    assert_eq!(m.width(), 12);
    assert_eq!(out.used, vec![0, 1]);
    let consensus = m.consensus.clone().unwrap_or_default();
    assert!(
        consensus.starts_with("AGGTCA") || consensus.starts_with("TGACCT"),
        "{consensus}"
    );
    // Without gaps, the merged motif has none.
    let off = MergeOptions {
        use_gaps: false,
        ..Default::default()
    };
    assert!(
        !merge_motifs_detailed(&[dr(1, 3, "a"), dr(3, 5, "b")], &off)
            .unwrap()
            .motif
            .gap
            .isgapped
    );
}

#[test]
fn mixing_gapped_and_ungapped_motifs_fills_the_gaps_with_a_warning() {
    let fixed = from_strings(
        &["AGGTCANNNNAGGTCA"],
        None,
        &CreateOptions {
            name: "DR4".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let out = merge_motifs_detailed(&[dr(2, 6, "DR"), fixed], &MergeOptions::default()).unwrap();
    assert!(!out.motif.gap.isgapped);
    assert_eq!(out.motif.width(), 16);
    assert_eq!(out.warnings.len(), 1);
    assert!(out.warnings[0].contains("DR"), "{}", out.warnings[0]);
    // merge_similar passes the warning on, naming the cluster.
    let fixed2 = from_strings(
        &["AGGTCANNNNAGGTCA"],
        None,
        &CreateOptions {
            name: "DR4".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let merged = merge_similar(&[dr(2, 6, "DR"), fixed2], &MergeSimilarOptions::default()).unwrap();
    assert_eq!(merged.motifs.len(), 1);
    assert!(merged.warnings.iter().any(|w| w.starts_with("cluster 0:")));
}

#[test]
fn cwms_are_compared_with_their_gaps_too() {
    let cwm = |s: &str, name: &str| {
        let cols: Vec<[f64; 4]> = s
            .bytes()
            .map(|c| {
                let mut col = [0.0; 4];
                if let Some(i) = b"ACGT".iter().position(|&l| l == c) {
                    col[i] = 0.5;
                }
                col
            })
            .collect();
        let opts = CreateOptions {
            name: name.into(),
            motif_type: MotifType::Cwm,
            ..Default::default()
        };
        from_typed_matrix(
            Matrix::from_columns(&cols).unwrap(),
            None,
            Some(Alphabet::Dna),
            MotifType::Cwm,
            &opts,
        )
        .unwrap()
    };
    let gapped = add_gap(&cwm("AGGTCAAGGTCA", "DR"), Some(&[6]), &[1], &[6]).unwrap();
    let fixed = cwm("AGGTCANNNNAGGTCA", "DR4");
    let rows = compare_cwms(&[gapped, fixed], &[0], &CwmCompareOptions::default()).unwrap();
    assert_eq!(rows[1].query_gaps, vec![4]);
    assert!(rows[1].score > 0.999, "{}", rows[1].score);
}
