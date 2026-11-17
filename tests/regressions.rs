//! Bugs found in review, each checked so that it stays fixed.

use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::enrich::{EnrichOptions, enrich_motifs};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::gapped::{GappedFinderOptions, gapped_motif_finder};
use universalmotif::genome::Genome;
use universalmotif::kmer::WindowOptions;
use universalmotif::merge::{MergeOptions, merge_motifs_detailed};
use universalmotif::ops::add_gap;
use universalmotif::shuffle::{ShuffleOptions, create_sequences, shuffle_sequences};

#[test]
fn merging_keeps_motifs_narrower_than_min_overlap() {
    let m = from_strings(&["TGAC"], None, &CreateOptions::default()).unwrap();
    let merged = merge_motifs_detailed(&[m.clone(), m], &MergeOptions::default()).unwrap();
    assert_eq!(merged.used, vec![0, 1]);
}

#[test]
fn gaps_out_of_order_are_refused() {
    let m = from_strings(&["TGACGTCATG"], None, &CreateOptions::default()).unwrap();
    assert!(add_gap(&m, Some(&[5, 2]), &[1], &[3]).is_err());
    assert!(add_gap(&m, Some(&[2, 5]), &[1], &[3]).is_ok());
}

#[test]
fn blank_lines_inside_a_fasta_record_are_refused() {
    let dir = std::env::temp_dir().join(format!("um-regress-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("blank.fa");
    std::fs::write(&path, ">chr1\nACGT\n\nTTTT\n").unwrap();
    let opened = Genome::open(&path);
    let fetched = opened.and_then(|g| g.fetch("chr1", 0, 8));
    assert!(fetched.is_err(), "{fetched:?}");
    // Blank lines between records are fine.
    std::fs::write(&path, ">chr1\nACGT\n\n>chr2\nTTTT\n").unwrap();
    let _ = std::fs::remove_file(dir.join("blank.fa.fai"));
    let g = Genome::open(&path).unwrap();
    assert_eq!(g.fetch("chr2", 0, 4).unwrap(), "TTTT");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn windows_longer_than_a_sequence_shuffle_it_whole() {
    let seq = "AACCGGTTAACCGGTTACGTACGTAAAACCCCGGGGTTTTACACGTGT".to_string();
    let opts = ShuffleOptions {
        k: 1,
        window: Some(WindowOptions {
            size: 60.0,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    let out = shuffle_sequences(std::slice::from_ref(&seq), &opts).unwrap();
    assert_ne!(out[0], seq);
    let mut a: Vec<char> = out[0].chars().collect();
    let mut b: Vec<char> = seq.chars().collect();
    a.sort_unstable();
    b.sort_unstable();
    assert_eq!(a, b);
}

#[test]
fn one_letter_alphabets_need_one_frequency() {
    assert!(create_sequences(&['A'], 1, 10, Some(&[0.5, 0.5]), 1, 1).is_err());
    assert_eq!(
        create_sequences(&['A'], 1, 5, None, 1, 1).unwrap(),
        vec!["AAAAA"]
    );
}

#[test]
fn an_empty_sequence_does_not_stop_shuffled_backgrounds() {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 60, 100, None, 3, 1).unwrap();
    seqs.push(String::new());
    let finder = FinderOptions {
        min_width: 6,
        max_width: 6,
        ..Default::default()
    };
    assert!(motif_finder(&seqs, None, &finder).is_ok());
    let gapped = GappedFinderOptions {
        k: vec![5],
        nmotifs: 1,
        ..Default::default()
    };
    assert!(gapped_motif_finder(&seqs, None, &gapped).is_ok());
    let m = from_strings(&["TGACTCA"], None, &CreateOptions::default()).unwrap();
    assert!(enrich_motifs(&[m], &seqs, None, &EnrichOptions::default()).is_ok());
}

#[test]
fn stability_with_no_sequences_is_an_error() {
    use universalmotif::select::{StabilityOptions, motif_stability};
    let m = from_strings(&["ACGTACGT"], None, &CreateOptions::default()).unwrap();
    let seqs = vec!["ACGTACGTACGTACGTAAAA".to_string(); 10];
    let none: Vec<String> = Vec::new();
    let opts = (FinderOptions::default(), StabilityOptions::default());
    assert!(
        motif_stability(
            std::slice::from_ref(&m),
            &seqs,
            Some(&none),
            &opts.0,
            &opts.1
        )
        .is_err()
    );
    assert!(motif_stability(&[m], &none, None, &opts.0, &opts.1).is_err());
}

#[test]
fn gapped_comparisons_take_a_query_twice() {
    use universalmotif::compare::{CompareOptions, compare_motifs};
    let o = CreateOptions {
        nsites: Some(20.0),
        ..Default::default()
    };
    let a = from_strings(&["AGGTCAAGGTCA"], None, &o).unwrap();
    let b = from_strings(&["CACGTGAC"], None, &o).unwrap();
    let g = add_gap(&a, Some(&[6]), &[1], &[3]).unwrap();
    let rows = compare_motifs(&[g, b], &[0, 0], &CompareOptions::default()).unwrap();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0].score, rows[2].score);
}

#[test]
fn refining_a_motif_marked_gapped_without_gaps_treats_it_as_ungapped() {
    use universalmotif::motif::{GapInfo, Motif};
    use universalmotif::select::{RefineOptions, refine_motifs};
    let o = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let mut p = from_strings(&["TTGACTCATG"], None, &o)
        .unwrap()
        .into_parts();
    p.gap = GapInfo {
        isgapped: true,
        ..Default::default()
    };
    let m = Motif::from_parts(p).unwrap();
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 100, 100, None, 3, 1).unwrap();
    for s in seqs.iter_mut() {
        s.replace_range(40..50, "TTGACTCATG");
    }
    let refined = refine_motifs(&[m], &seqs, None, &RefineOptions::default()).unwrap();
    assert!(refined[0].sites >= 100);
}

#[test]
fn missing_contribution_scores_count_as_zero() {
    use universalmotif::contributions::{ContribScanOptions, ContribThreshold, scan_hypothetical};
    let gat = from_strings(&["GATTACA"], None, &CreateOptions::default()).unwrap();
    let len = 20;
    let mut scores = vec![0.01f32; 4 * len];
    scores[2 * len + 5] = f32::NAN;
    let t = ContribThreshold {
        min_similarity: 0.9,
        min_importance: 0.5,
    };
    let hits = scan_hypothetical(
        &[gat],
        &scores,
        &scores,
        &[len],
        len,
        &[t],
        &ContribScanOptions::default(),
    )
    .unwrap();
    assert!(
        hits.iter()
            .all(|h| h.similarity.is_finite() && h.importance.is_finite())
    );
    assert!(hits.is_empty());
}
