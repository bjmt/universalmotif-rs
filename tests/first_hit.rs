//! `ScanOptions::first_hit` keeps the first hit a full scan reports for each
//! motif and sequence, and presence from `motif_matches` is the same whether
//! or not scanning stops early.

use std::collections::BTreeMap;

use universalmotif::Alphabet;
use universalmotif::create::{CreateOptions, from_strings, random};
use universalmotif::motif::Motif;
use universalmotif::ops::add_gap;
use universalmotif::packed::EncodedSequences;
use universalmotif::scan::{
    HitTable, MatchValue, NoOverlaps, ScanOptions, motif_matches, scan_encoded, scan_fasta,
    scan_table,
};
use universalmotif::shuffle::create_sequences;

fn dna_motifs() -> Vec<Motif> {
    let o = CreateOptions {
        bkg: Some(vec![0.25; 4]),
        ..Default::default()
    };
    let mut motifs: Vec<Motif> = (0..16)
        .map(|i| random(6 + i % 9, Alphabet::Dna, &o, 100 + i as u64).unwrap())
        .collect();
    // Gapped motifs, with one and two gaps.
    motifs.push(add_gap(&motifs[3], Some(&[4]), &[0], &[4]).unwrap());
    motifs.push(add_gap(&motifs[8], Some(&[3, 7]), &[1, 0], &[3, 2]).unwrap());
    // A palindrome, hit on both strands at the same start.
    let s = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    motifs.push(from_strings(&["CACGTG"], None, &s).unwrap());
    motifs
}

/// Short sequences, some without hits, and long ones scanned in several
/// chunks, with unknown letters splitting them into stretches.
fn dna_sequences() -> Vec<String> {
    let mut seqs = create_sequences(&['A', 'C', 'G', 'T'], 300, 400, None, 7, 1).unwrap();
    let mut long = create_sequences(&['A', 'C', 'G', 'T'], 2, 150_000, None, 8, 1).unwrap();
    long[0].replace_range(70_000..70_500, &"N".repeat(500));
    long[1] = long[1].to_ascii_lowercase();
    // A repeat few motifs hit fills the first chunks, so that many first
    // hits are in a later one.
    long.push("AC".repeat(70_000) + &long[0][..30_000]);
    seqs.extend(long);
    seqs.push("ACGT".into());
    seqs.push("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".into());
    seqs
}

/// The rows of a full scan that come first for their motif and sequence.
fn first_rows(t: &HitTable) -> Vec<usize> {
    (0..t.len())
        .filter(|&i| i == 0 || (t.motif[i], t.sequence[i]) != (t.motif[i - 1], t.sequence[i - 1]))
        .collect()
}

type Row = (u32, u32, i64, i64, bool, f64, f64, f64, Option<Vec<usize>>);

fn rows(t: &HitTable, which: &[usize]) -> Vec<Row> {
    which
        .iter()
        .map(|&i| {
            (
                t.motif[i],
                t.sequence[i],
                t.start[i],
                t.end[i],
                t.minus[i],
                t.score[i],
                t.score_pct[i],
                t.pvalue[i],
                t.gap_sizes(i),
            )
        })
        .collect()
}

fn all_rows(t: &HitTable) -> Vec<Row> {
    rows(t, &(0..t.len()).collect::<Vec<_>>())
}

#[test]
fn first_hit_is_the_first_row_of_a_full_scan() {
    let motifs = dna_motifs();
    let seqs = dna_sequences();
    let mut ties = 0;
    for (rc, correct_gaps) in [(true, false), (false, false), (true, true)] {
        let full_opts = ScanOptions {
            pvalue: 1e-3,
            rc,
            correct_gaps,
            nthreads: 2,
            ..Default::default()
        };
        let full = scan_table(&motifs, &seqs, &full_opts).unwrap();
        let want = rows(&full, &first_rows(&full));
        assert!(want.len() > 1000, "only {} first hits", want.len());
        assert!(
            want.iter().any(|r| r.2 >= 1 << 16),
            "no first hits past a chunk"
        );
        // Pairs whose first start has a hit on both strands.
        ties += first_rows(&full)
            .iter()
            .filter(|&&i| {
                i + 1 < full.len()
                    && (full.motif[i + 1], full.sequence[i + 1], full.start[i + 1])
                        == (full.motif[i], full.sequence[i], full.start[i])
                    && full.minus[i + 1]
            })
            .count();
        let mut tables = Vec::new();
        for nthreads in [1, 2] {
            let opts = ScanOptions {
                first_hit: true,
                nthreads,
                ..full_opts.clone()
            };
            let got = scan_table(&motifs, &seqs, &opts).unwrap();
            assert_eq!(
                all_rows(&got),
                want,
                "rc {rc}, correct_gaps {correct_gaps}, nthreads {nthreads}"
            );
            let enc = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
            assert_eq!(scan_encoded(&motifs, &enc, &opts).unwrap(), got);
            tables.push(got);
        }
        assert_eq!(tables[0], tables[1]);
    }
    assert!(ties > 0, "no first hits on both strands at one start");
}

#[test]
fn first_hit_in_other_alphabets_and_fasta_files() {
    let o = CreateOptions {
        bkg: Some(vec![0.05; 20]),
        ..Default::default()
    };
    let mut motifs: Vec<Motif> = (0..6)
        .map(|i| random(4 + i, Alphabet::Aa, &o, 200 + i as u64).unwrap())
        .collect();
    motifs.push(add_gap(&motifs[4], Some(&[3]), &[0], &[3]).unwrap());
    let letters = Alphabet::Aa.letters().to_vec();
    let mut seqs = create_sequences(&letters, 200, 300, None, 9, 1).unwrap();
    seqs.extend(create_sequences(&letters, 2, 150_000, None, 10, 1).unwrap());
    let base = ScanOptions {
        pvalue: 1e-3,
        nthreads: 2,
        ..Default::default()
    };
    let full = scan_table(&motifs, &seqs, &base).unwrap();
    let want = rows(&full, &first_rows(&full));
    assert!(want.len() > 200, "only {} first hits", want.len());
    let opts = ScanOptions {
        first_hit: true,
        ..base
    };
    let got = scan_table(&motifs, &seqs, &opts).unwrap();
    assert_eq!(all_rows(&got), want);

    let dna = dna_motifs();
    let dna_seqs = dna_sequences();
    let dir = std::env::temp_dir().join(format!("um-first-hit-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("seqs.fa");
    let records: Vec<(String, &String)> = dna_seqs
        .iter()
        .enumerate()
        .map(|(i, s)| (format!("s{i}"), s))
        .collect();
    std::fs::write(&path, universalmotif::io::fasta::write(&records, 80)).unwrap();
    let opts = ScanOptions {
        pvalue: 1e-3,
        first_hit: true,
        nthreads: 2,
        ..Default::default()
    };
    let want = scan_table(&dna, &dna_seqs, &opts).unwrap();
    let f = scan_fasta(&dna, &path, &opts, true).unwrap();
    assert_eq!(f.table, want);
    for i in 0..want.len() {
        assert_eq!(f.matched(i), want.matched(i, &dna_seqs));
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn first_hit_cannot_be_combined_with_no_overlaps() {
    let opts = ScanOptions {
        first_hit: true,
        no_overlaps: Some(NoOverlaps::default()),
        ..Default::default()
    };
    let err = scan_table(&dna_motifs(), &["ACGTACGTACGT"], &opts).unwrap_err();
    assert!(err.to_string().contains("no_overlaps"), "{err}");
}

/// Presence of each motif in each sequence from the hits of a full scan.
fn presence(t: &HitTable) -> Vec<(u32, u32)> {
    let mut keys: Vec<(u32, u32)> = (0..t.len()).map(|i| (t.sequence[i], t.motif[i])).collect();
    keys.sort_unstable();
    keys.dedup();
    keys
}

#[test]
fn presence_is_the_same_with_early_exit() {
    let motifs = dna_motifs();
    let seqs = dna_sequences();
    let enc = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
    let overlaps = [
        None,
        Some(NoOverlaps::default()),
        Some(NoOverlaps {
            ignore_strand: true,
            ..Default::default()
        }),
        // Hits of one motif can be dropped for another's here, so presence
        // comes from a full scan.
        Some(NoOverlaps {
            ignore_motif: true,
            ..Default::default()
        }),
    ];
    for (k, no_overlaps) in overlaps.into_iter().enumerate() {
        for rc in [true, false] {
            if k > 0 && !rc {
                continue;
            }
            let mut results = Vec::new();
            for nthreads in [1, 2] {
                let opts = ScanOptions {
                    pvalue: 1e-3,
                    rc,
                    no_overlaps,
                    nthreads,
                    ..Default::default()
                };
                let want = presence(&scan_encoded(&motifs, &enc, &opts).unwrap());
                let got = motif_matches(&motifs, &enc, &opts, MatchValue::Matches).unwrap();
                let keys: Vec<(u32, u32)> = got
                    .sequence
                    .iter()
                    .copied()
                    .zip(got.motif.iter().copied())
                    .collect();
                assert_eq!(keys, want, "{no_overlaps:?}, rc {rc}, nthreads {nthreads}");
                assert!(got.value.iter().all(|&v| v == 1.0));
                results.push(got);
            }
            assert_eq!(results[0], results[1]);
        }
    }
}

#[test]
fn matches_follow_first_hit() {
    let motifs = dna_motifs();
    let seqs = dna_sequences();
    let enc = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
    let opts = ScanOptions {
        pvalue: 1e-3,
        first_hit: true,
        nthreads: 2,
        ..Default::default()
    };
    let t = scan_encoded(&motifs, &enc, &opts).unwrap();
    let want: BTreeMap<(u32, u32), f64> = (0..t.len())
        .map(|i| ((t.sequence[i], t.motif[i]), t.score[i]))
        .collect();
    assert_eq!(want.len(), t.len());
    for value in [MatchValue::Matches, MatchValue::Counts, MatchValue::Scores] {
        let got = motif_matches(&motifs, &enc, &opts, value).unwrap();
        let keys: Vec<(u32, u32)> = got
            .sequence
            .iter()
            .copied()
            .zip(got.motif.iter().copied())
            .collect();
        assert_eq!(keys, want.keys().copied().collect::<Vec<_>>(), "{value:?}");
        for (k, v) in keys.iter().zip(&got.value) {
            let expected = match value {
                MatchValue::Scores => want[k],
                _ => 1.0,
            };
            assert_eq!(*v, expected, "{value:?} {k:?}");
        }
    }
}
