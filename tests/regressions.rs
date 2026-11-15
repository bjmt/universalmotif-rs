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
