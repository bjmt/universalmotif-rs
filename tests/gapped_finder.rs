//! Discovery of gapped motifs planted in random sequences.

use universalmotif::gapped::{GappedFinderOptions, gapped_motif_finder};

/// Random DNA from a simple generator, so that the test needs no data.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn dna(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize] as char)
            .collect()
    }
}

/// Sequences with `AGGTCA`, a spacer of 1 to 5 letters and `AGGTCA` again
/// planted in every other one.
fn direct_repeats(n: usize, seed: u64) -> Vec<String> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|i| {
            let mut s = rng.dna(150);
            if i % 2 == 0 {
                let gap = 1 + (rng.next() % 5) as usize;
                let site = format!("AGGTCA{}AGGTCA", rng.dna(gap));
                let at = (rng.next() % 100) as usize;
                s.replace_range(at..at + site.len(), &site);
            }
            s
        })
        .collect()
}

fn options() -> GappedFinderOptions {
    GappedFinderOptions {
        k: vec![5],
        nmotifs: 1,
        ..Default::default()
    }
}

#[test]
fn finds_a_planted_direct_repeat() {
    let seqs = direct_repeats(300, 7);
    let found = gapped_motif_finder(&seqs, None, &options()).unwrap();
    assert_eq!(found.len(), 1);
    let m = &found[0];
    assert!(m.motif.gap.isgapped);
    let (lo, hi) = (m.motif.gap.mingap[0], m.motif.gap.maxgap[0]);
    assert!(lo <= 2 && (4..=7).contains(&hi), "spacer {lo}-{hi}");
    let consensus = m.motif.consensus.as_deref().unwrap();
    assert!(
        consensus.contains("GGTC") || consensus.contains("GACC"),
        "{consensus}"
    );
    assert!(m.seqs_pos > 100 && m.ln_spacing_pvalue < (1e-3f64).ln());
}

#[test]
fn finds_it_without_holding_sequences_out() {
    let seqs = direct_repeats(300, 7);
    let opts = GappedFinderOptions {
        holdout: 0.0,
        ..options()
    };
    let found = gapped_motif_finder(&seqs, None, &opts).unwrap();
    assert_eq!(found.len(), 1);
    assert!(found[0].motif.gap.isgapped);
}

#[test]
fn finds_nothing_in_random_sequences() {
    let mut rng = Lcg(11);
    let seqs: Vec<String> = (0..300).map(|_| rng.dna(150)).collect();
    assert!(
        gapped_motif_finder(&seqs, None, &options())
            .unwrap()
            .is_empty()
    );
}

#[test]
fn checks_its_options() {
    let seqs = direct_repeats(20, 1);
    let bad = |o: GappedFinderOptions| gapped_motif_finder(&seqs, None, &o).is_err();
    assert!(bad(GappedFinderOptions {
        k: vec![7],
        ..Default::default()
    }));
    assert!(bad(GappedFinderOptions {
        min_gap: 5,
        max_gap: 2,
        ..Default::default()
    }));
    assert!(bad(GappedFinderOptions {
        hit_pvalue: 0.0,
        ..Default::default()
    }));
    assert!(bad(GappedFinderOptions {
        holdout: 0.95,
        ..Default::default()
    }));
    assert!(gapped_motif_finder::<String>(&[], None, &options()).is_err());
}

#[test]
fn spaced_pairs_count_each_sequence_once() {
    use universalmotif::finder::spaced_pairs;
    // AA then CC with no letters between, twice in the first sequence.
    let pos = ["AACCTTTTAACC", "GGGGGGGGGGGG"];
    let neg = ["TTTTTTTTTTTT"];
    let find = |rows: &[universalmotif::finder::SpacedPair], left: u64, right: u64| {
        rows.iter()
            .find(|p| p.left == left && p.right == right && p.gap == 0 && p.width == 1)
            .map(|p| (p.positives, p.negatives))
    };
    let (aa, cc, ca) = (0, 5, 4);
    let exact = spaced_pairs(&pos, &neg, 2, 0, 2, false, 1, false, 1).unwrap();
    assert_eq!(find(&exact, aa, cc), Some((1, 0)));
    assert_eq!(find(&exact, aa, ca), None);
    // One mismatch lets AA CA count where AA CC is.
    let loose = spaced_pairs(&pos, &neg, 2, 0, 2, false, 1, true, 1).unwrap();
    assert_eq!(find(&loose, aa, ca), Some((1, 0)));
    assert!(spaced_pairs(&pos, &neg, 7, 0, 2, false, 1, false, 1).is_err());
}
