//! motif_matches summarises the same hits scan_table reports.

use std::collections::BTreeMap;

use universalmotif::Alphabet;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::ops::add_gap;
use universalmotif::packed::EncodedSequences;
use universalmotif::scan::{MatchValue, NoOverlaps, ScanOptions, motif_matches, scan_encoded};
use universalmotif::shuffle::create_sequences;

#[test]
fn matches_agree_with_the_hit_table() {
    let o = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let motifs = vec![
        from_strings(&["TGACTCA", "TGAGTCA"], None, &o).unwrap(),
        from_strings(&["CACGTG"], None, &o).unwrap(),
        add_gap(
            &from_strings(&["TTGACGTCAA"], None, &o).unwrap(),
            Some(&[5]),
            &[1],
            &[3],
        )
        .unwrap(),
    ];
    // More sequence than one block, so blocks are summarised separately.
    let seqs = create_sequences(&['A', 'C', 'G', 'T'], 1200, 4000, None, 3, 1).unwrap();
    let enc = EncodedSequences::new(&seqs, &Alphabet::Dna, 1).unwrap();
    for no_overlaps in [
        None,
        Some(NoOverlaps {
            ignore_strand: true,
            ..Default::default()
        }),
    ] {
        let opts = ScanOptions {
            pvalue: 1e-3,
            no_overlaps,
            nthreads: 2,
            ..Default::default()
        };
        let hits = scan_encoded(&motifs, &enc, &opts).unwrap();
        let mut want: BTreeMap<(u32, u32), (f64, f64)> = BTreeMap::new();
        for i in 0..hits.len() {
            let e = want
                .entry((hits.sequence[i] as u32, hits.motif[i] as u32))
                .or_insert((0.0, f64::NEG_INFINITY));
            e.0 += 1.0;
            e.1 = e.1.max(hits.score[i]);
        }
        for value in [MatchValue::Matches, MatchValue::Counts, MatchValue::Scores] {
            let got = motif_matches(&motifs, &enc, &opts, value).unwrap();
            assert_eq!((got.n_sequences, got.n_motifs), (1200, 3));
            let keys: Vec<(u32, u32)> = got
                .sequence
                .iter()
                .copied()
                .zip(got.motif.iter().copied())
                .collect();
            assert_eq!(keys, want.keys().copied().collect::<Vec<_>>(), "{value:?}");
            for (k, v) in keys.iter().zip(&got.value) {
                let (count, best) = want[k];
                let expected = match value {
                    MatchValue::Matches => 1.0,
                    MatchValue::Counts => count,
                    MatchValue::Scores => best,
                };
                assert!(
                    (v - expected).abs() < 1e-9,
                    "{value:?} {k:?}: {v} vs {expected}"
                );
            }
        }
    }
}
