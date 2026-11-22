//! Discovery of pairs of motifs planted in random sequences.

use universalmotif::Strand;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::pair_finder::{PairFinderOptions, motif_pair_finder};

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

/// 1,500 sequences of 200 letters. `TGTCGG`, a spacer of 7 letters and
/// `TTGACT` are planted in every `every`-th one, `TGTCGG` alone in a third
/// of them and `TTGACT` alone in a sixth, so that both sites are common and
/// the pair rare.
fn planted(every: usize, seed: u64) -> Vec<String> {
    let mut rng = Lcg(seed);
    (0..1500)
        .map(|i| {
            let mut s = rng.dna(200);
            let mut put = |site: &str, rng: &mut Lcg| {
                let at = (rng.next() % (200 - site.len() as u64)) as usize;
                s.replace_range(at..at + site.len(), site);
            };
            if i % 3 == 1 {
                put("TGTCGG", &mut rng);
            }
            if i % 6 == 2 {
                put("TTGACT", &mut rng);
            }
            if every > 0 && i % every == 0 {
                let site = format!("TGTCGG{}TTGACT", rng.dna(7));
                put(&site, &mut rng);
            }
            s
        })
        .collect()
}

fn anchor() -> universalmotif::Motif {
    from_strings(
        &["TGTCGG"],
        None,
        &CreateOptions {
            name: "ARF".into(),
            nsites: Some(20.0),
            pseudocount: 1.0,
            ..Default::default()
        },
    )
    .unwrap()
}

fn options() -> PairFinderOptions {
    PairFinderOptions {
        nmotifs: 4,
        nthreads: 2,
        ..Default::default()
    }
}

#[test]
fn finds_a_planted_pair_in_a_twentieth_of_the_sequences() {
    let seqs = planted(20, 3);
    let res = motif_pair_finder(&seqs, None, &[anchor()], &options()).unwrap();
    assert_eq!(res.parts[0].name, "ARF");
    let top = res.pairs.first().expect("a pair");
    assert_eq!(res.parts[top.left].name, "ARF", "{top:?}");
    assert_eq!(top.left_strand, Strand::Plus);
    // The partner found may be read on either strand, which sets the
    // arrangement.
    let expected = if top.right_strand == Strand::Plus {
        "tandem"
    } else {
        "convergent"
    };
    assert_eq!(top.arrangement, expected);
    assert!((top.min_gap..=top.max_gap).contains(&7) && top.max_gap - top.min_gap <= 2);
    assert!(top.log10_pvalue < -5.0 && top.fold > 5.0, "{top:?}");
    assert!(!top.repeat);
    let consensus = top.motif.consensus.as_deref().unwrap();
    assert!(
        consensus.starts_with("TGTCGG") && consensus[6..].contains("TTGACT"),
        "{consensus}"
    );
    assert!(top.motif.gap.isgapped && top.motif.gap.gaploc == [6]);
    // Shifted copies of the sites are merged into one pair.
    assert_eq!(
        res.pairs
            .iter()
            .filter(|p| p.min_gap.abs_diff(7) <= 3)
            .count(),
        1
    );
}

#[test]
fn finds_it_without_an_anchor() {
    let seqs = planted(10, 4);
    let res = motif_pair_finder(&seqs, None, &[], &options()).unwrap();
    let top = res.pairs.first().expect("a pair");
    let consensus = top.motif.consensus.as_deref().unwrap();
    assert!(
        consensus.contains("CCGACA") || consensus.contains("TGTCGG"),
        "{consensus}"
    );
    assert!(
        consensus.contains("AGTCAA") || consensus.contains("TTGACT"),
        "{consensus}"
    );
}

#[test]
fn finds_nothing_in_random_sequences() {
    let seqs = planted(0, 5);
    let res = motif_pair_finder(&seqs, None, &[anchor()], &options()).unwrap();
    assert!(res.pairs.is_empty(), "{:?}", res.pairs);
}

#[test]
fn flags_tandem_repeats() {
    let mut rng = Lcg(6);
    let seqs: Vec<String> = (0..1500)
        .map(|i| {
            let mut s = rng.dna(200);
            if i % 10 == 0 {
                let unit = "CCCTAAA".repeat(4);
                let at = (rng.next() % 150) as usize;
                s.replace_range(at..at + unit.len(), &unit);
            }
            s
        })
        .collect();
    let opts = PairFinderOptions {
        mask: false,
        ..options()
    };
    let res = motif_pair_finder(&seqs, None, &[], &opts).unwrap();
    let same: Vec<_> = res.pairs.iter().filter(|p| p.left == p.right).collect();
    assert!(
        !same.is_empty() && same.iter().all(|p| p.repeat),
        "{:?}",
        res.pairs
    );
}

#[test]
fn checks_its_options() {
    let seqs = planted(20, 3);
    let bad = |opts: PairFinderOptions, motifs: &[universalmotif::Motif]| {
        motif_pair_finder(&seqs, None, motifs, &opts).is_err()
    };
    assert!(bad(
        PairFinderOptions {
            nmotifs: 0,
            ..options()
        },
        &[]
    ));
    assert!(bad(
        PairFinderOptions {
            max_gap: 32,
            ..options()
        },
        &[anchor()]
    ));
    assert!(bad(
        PairFinderOptions {
            max_window: 0,
            ..options()
        },
        &[anchor()]
    ));
    assert!(bad(
        PairFinderOptions {
            holdout: 0.0,
            ..options()
        },
        &[anchor()]
    ));
    assert!(bad(
        PairFinderOptions {
            stop_pvalue: 0.0,
            ..options()
        },
        &[anchor()]
    ));
    assert!(motif_pair_finder::<String>(&[], None, &[anchor()], &options()).is_err());
}
