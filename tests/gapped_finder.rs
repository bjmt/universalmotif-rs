//! Discovery of gapped motifs planted in random sequences.

use universalmotif::gapped::{GappedFinderOptions, Mismatch, gapped_motif_finder};

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
    assert!(m.seqs_pos > 100 && m.spacing_pvalue < 1e-3 && m.log10_spacing_pvalue < -3.0);
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
fn takes_the_longest_gaps_and_letters_outside_ascii() {
    let mut seqs = direct_repeats(300, 7);
    seqs[0].insert(0, 'é');
    let found = find(
        &seqs,
        GappedFinderOptions {
            max_gap: 31,
            ..options()
        },
    );
    assert_eq!(found.len(), 1);
    assert!(found[0].motif.gap.isgapped);
    // Spacers beyond 31 would give motifs too costly to compare.
    let too_long = GappedFinderOptions {
        max_gap: 32,
        ..options()
    };
    assert!(gapped_motif_finder(&seqs, None, &too_long).is_err());
}

fn find(seqs: &[String], opts: GappedFinderOptions) -> Vec<universalmotif::gapped::GappedMotif> {
    gapped_motif_finder(seqs, None, &opts).unwrap()
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

/// Sequences with `AGGTCA`, three letters and `TGTCGG` planted in every
/// sixth one, each site with one letter changed, so that no pair of 6-mers
/// is common enough to count exactly.
fn varied_sites(n: usize, seed: u64) -> Vec<String> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|i| {
            let mut s = rng.dna(150);
            if i % 6 == 0 {
                let mut site: Vec<u8> = format!("AGGTCA{}TGTCGG", rng.dna(3)).into_bytes();
                let mut at = (rng.next() % 12) as usize;
                if at >= 6 {
                    at += 3;
                }
                let old = site[at];
                site[at] = *b"ACGT"
                    .iter()
                    .filter(|&&c| c != old)
                    .nth((rng.next() % 3) as usize)
                    .unwrap();
                let pos = (rng.next() % 100) as usize;
                s.replace_range(pos..pos + site.len(), std::str::from_utf8(&site).unwrap());
            }
            s
        })
        .collect()
}

#[test]
fn retries_with_a_mismatch_only_when_nothing_is_found() {
    let opts = |mismatch| GappedFinderOptions {
        k: vec![6],
        nmotifs: 1,
        mismatch,
        ..Default::default()
    };
    let varied = varied_sites(300, 5);
    assert!(find(&varied, opts(Mismatch::Off)).is_empty());
    // The retry keeps motifs at half of stop_pvalue.
    let loose = find(
        &varied,
        GappedFinderOptions {
            stop_pvalue: 5e-4,
            ..opts(Mismatch::On)
        },
    );
    assert_eq!(loose.len(), 1);
    assert_eq!(loose[0].seed_pair, "AGGTCA[3-3]TGTCGG");
    assert_eq!(find(&varied, opts(Mismatch::Auto)), loose);
    // A search that finds something is not run again.
    let exact = direct_repeats(300, 7);
    assert_eq!(
        find(&exact, opts(Mismatch::Auto)),
        find(&exact, opts(Mismatch::Off))
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

#[test]
fn summarizes_the_spacing_of_its_parts() {
    let seqs = direct_repeats(300, 7);
    let found = gapped_motif_finder(&seqs, None, &options()).unwrap();
    let m = &found[0];
    let best = m.best_spacer.expect("a best spacer");
    let (lo, hi) = (m.motif.gap.mingap[0], m.motif.gap.maxgap[0]);
    assert!((lo..=hi).contains(&best), "best {best} outside {lo}-{hi}");
    assert!(
        m.spacer_pairs >= 5 && m.spacer_fold > 2.0,
        "{} pairs, fold {}",
        m.spacer_pairs,
        m.spacer_fold
    );
    assert!(m.strong_spacers.contains(&best));
    assert!(m.strong_spacers.iter().all(|&g| g <= 31));
}
