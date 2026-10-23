//! Shuffling, random sequences, site sampling and implanting. These use a
//! different random number generator from the R package, so they are checked
//! by their properties rather than against R's output.

use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::kmer::{WindowOptions, count_kmers, kmers, sequence_letters};
use universalmotif::sample::{ImplantMode, ImplantOptions, implant_motifs, sample_sites};
use universalmotif::shuffle::{
    ShuffleMethod, ShuffleOptions, create_sequences, shuffle_motifs, shuffle_sequences,
};
use universalmotif::{Alphabet, Motif, MotifType};

const DNA: [char; 4] = ['A', 'C', 'G', 'T'];

fn random_seqs(n: usize, len: usize, seed: u64) -> Vec<String> {
    create_sequences(&DNA, n, len, None, seed, 1).unwrap()
}

fn sorted(s: &str) -> Vec<char> {
    let mut c: Vec<char> = s.chars().collect();
    c.sort_unstable();
    c
}

fn opts(k: usize, method: ShuffleMethod, seed: u64) -> ShuffleOptions {
    ShuffleOptions {
        k,
        method,
        seed,
        ..Default::default()
    }
}

#[test]
fn euler_shuffle_keeps_kmers() {
    let seqs = random_seqs(20, 120, 3);
    for k in 1..=4 {
        let out = shuffle_sequences(&seqs, &opts(k, ShuffleMethod::Euler, 9)).unwrap();
        assert_ne!(out, seqs, "k={k}: nothing was shuffled");
        for (a, b) in seqs.iter().zip(&out) {
            let alph = sequence_letters(&[a]);
            assert_eq!(count_kmers(a, k, &alph), count_kmers(b, k, &alph), "k={k}");
            if k > 1 {
                // The first and last k-1 letters stay put.
                assert_eq!(a[..k - 1], b[..k - 1], "k={k}");
                assert_eq!(a[a.len() - k + 1..], b[b.len() - k + 1..], "k={k}");
            }
        }
    }
}

#[test]
fn shuffles_keep_letters_and_are_seeded() {
    let seqs = random_seqs(10, 101, 4);
    for method in [ShuffleMethod::Euler, ShuffleMethod::Linear] {
        for k in [1, 2, 3] {
            let a = shuffle_sequences(&seqs, &opts(k, method, 1)).unwrap();
            let b = shuffle_sequences(&seqs, &opts(k, method, 1)).unwrap();
            let c = shuffle_sequences(&seqs, &opts(k, method, 2)).unwrap();
            assert_eq!(a, b, "{method:?} k={k}: same seed");
            assert_ne!(a, c, "{method:?} k={k}: different seed");
            for (s, t) in seqs.iter().zip(&a) {
                assert_eq!(sorted(s), sorted(t), "{method:?} k={k}");
            }
        }
    }
    // Different thread counts give the same result.
    let mut o = opts(2, ShuffleMethod::Euler, 5);
    let one = shuffle_sequences(&seqs, &o).unwrap();
    o.nthreads = 3;
    assert_eq!(one, shuffle_sequences(&seqs, &o).unwrap());
}

#[test]
fn linear_shuffle_moves_blocks() {
    let seqs = random_seqs(5, 62, 6);
    let k = 3;
    let out = shuffle_sequences(&seqs, &opts(k, ShuffleMethod::Linear, 8)).unwrap();
    for (a, b) in seqs.iter().zip(&out) {
        let blocks = |s: &str| {
            let mut v: Vec<String> = s.as_bytes()[..s.len() / k * k]
                .chunks(k)
                .map(|c| String::from_utf8(c.to_vec()).unwrap())
                .collect();
            v.sort();
            v
        };
        assert_eq!(blocks(a), blocks(b));
        assert_eq!(a[a.len() / k * k..], b[b.len() / k * k..]);
    }
}

#[test]
fn markov_shuffle_follows_kmer_frequencies() {
    // A strongly structured sequence: after A always comes C, and so on.
    let seq = "ACGT".repeat(500);
    let out = shuffle_sequences(&[&seq], &opts(2, ShuffleMethod::Markov, 3)).unwrap();
    let counts = count_kmers(&out[0], 2, &DNA);
    let kl = kmers(&DNA, 2);
    for (kmer, c) in kl.iter().zip(&counts) {
        if ["AC", "CG", "GT", "TA"].contains(&kmer.as_str()) {
            assert!(*c > 400, "{kmer}: {c}");
        } else {
            assert_eq!(*c, 0, "{kmer}");
        }
    }
    assert_eq!(out[0].len(), seq.len());
}

#[test]
fn windowed_shuffle_keeps_window_letters() {
    let seqs = random_seqs(4, 100, 12);
    let mut o = opts(1, ShuffleMethod::Euler, 2);
    o.window = Some(WindowOptions {
        size: 20.0,
        overlap: 0.0,
    });
    let out = shuffle_sequences(&seqs, &o).unwrap();
    assert_ne!(out, seqs);
    for (a, b) in seqs.iter().zip(&out) {
        for i in (0..100).step_by(20) {
            assert_eq!(sorted(&a[i..i + 20]), sorted(&b[i..i + 20]));
        }
    }
}

#[test]
fn created_sequences_follow_frequencies() {
    let seqs = create_sequences(&DNA, 50, 400, Some(&[0.1, 0.2, 0.3, 0.4]), 1, 2).unwrap();
    assert_eq!(seqs.len(), 50);
    assert!(seqs.iter().all(|s| s.len() == 400));
    let total: Vec<u64> = seqs
        .iter()
        .map(|s| count_kmers(s, 1, &DNA))
        .fold(vec![0; 4], |a, c| {
            a.iter().zip(&c).map(|(x, y)| x + y).collect()
        });
    for (c, f) in total.iter().zip([0.1, 0.2, 0.3, 0.4]) {
        let p = *c as f64 / 20000.0;
        assert!((p - f).abs() < 0.01, "{p} vs {f}");
    }
    assert_eq!(
        seqs,
        create_sequences(&DNA, 50, 400, Some(&[0.1, 0.2, 0.3, 0.4]), 1, 1).unwrap()
    );

    // Dinucleotide frequencies where only AA and TT occur after A and T.
    let mut f2 = vec![0.0; 16];
    f2[0] = 0.5; // AA
    f2[15] = 0.5; // TT
    let seqs = create_sequences(&DNA, 5, 50, Some(&f2), 3, 1).unwrap();
    for s in seqs {
        assert!(s.chars().all(|c| c == s.chars().next().unwrap()), "{s}");
    }
}

fn dna_motif(sites: &[&str], name: &str) -> Motif {
    from_strings(
        sites,
        Some(Alphabet::Dna),
        &CreateOptions {
            name: name.into(),
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn shuffled_motifs_keep_columns() {
    let a = dna_motif(&["ACGTAC", "ACGTTC", "AGGTAC"], "a");
    let b = dna_motif(&["TTTGCA", "TTAGCA"], "b");
    let c = universalmotif::convert::convert_type(
        &dna_motif(&["GGGCCC", "GGGCAC"], "c"),
        MotifType::Pwm,
        None,
        false,
    )
    .unwrap();
    let motifs = [a, b, c];
    let to_ppm =
        |m: &Motif| universalmotif::convert::convert_type(m, MotifType::Ppm, None, false).unwrap();
    let column_multiset = |ms: &[Motif]| {
        let mut v: Vec<Vec<u64>> = ms
            .iter()
            .flat_map(|m| {
                to_ppm(m)
                    .matrix
                    .columns()
                    .map(|c| c.iter().map(|x| (x * 1e9).round() as u64).collect())
                    .collect::<Vec<_>>()
            })
            .collect();
        v.sort();
        v
    };
    for k in [1, 2, 3] {
        let out = shuffle_motifs(&motifs, k, 4).unwrap();
        assert_eq!(out.len(), 3);
        for (o, m) in out.iter().zip(&motifs) {
            assert_eq!(o.width(), m.width());
            assert_eq!(o.motif_type, m.motif_type);
            assert_eq!(o.name, format!("{} [shuffled]", m.name));
        }
        assert_eq!(column_multiset(&out), column_multiset(&motifs), "k={k}");
        assert_eq!(out, shuffle_motifs(&motifs, k, 4).unwrap());
    }
}

#[test]
fn sampled_sites_follow_the_motif() {
    let m = dna_motif(&["AAGT", "ACGT", "ACGT", "ACGA"], "m");
    let sites = sample_sites(&m, 4000, 7).unwrap();
    assert!(sites.iter().all(|s| s.len() == 4));
    let ppm = universalmotif::convert::convert_type(&m, MotifType::Ppm, None, false).unwrap();
    for pos in 0..4 {
        for (l, letter) in DNA.iter().enumerate() {
            let p = sites
                .iter()
                .filter(|s| s.as_bytes()[pos] as char == *letter)
                .count() as f64
                / 4000.0;
            assert!((p - ppm.matrix.get(l, pos)).abs() < 0.03, "{pos} {letter}");
        }
    }
    assert_eq!(sites, sample_sites(&m, 4000, 7).unwrap());
}

#[test]
fn implants_are_placed_and_spaced() {
    let m = dna_motif(&["GGGGGGGG"], "g");
    let seqs: Vec<String> = vec!["A".repeat(200); 30];
    let o = ImplantOptions {
        mode: ImplantMode::Count(4),
        min_spacing: 5,
        seed: 3,
        ..Default::default()
    };
    let r = implant_motifs(std::slice::from_ref(&m), &seqs, &o).unwrap();
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    assert_eq!(r.implants.len(), 120);
    for (i, s) in r.sequences.iter().enumerate() {
        let mut here: Vec<_> = r.implants.iter().filter(|x| x.sequence == i).collect();
        here.sort_by_key(|x| x.start);
        for w in here.windows(2) {
            assert!(w[1].start > w[0].start + w[0].width - 1 + 5);
        }
        for x in &here {
            // Both strands of GGGGGGGG are G or C runs.
            let site = &s[x.start..x.start + x.width];
            assert_eq!(site, x.planted);
            assert!(site == "GGGGGGGG" || site == "CCCCCCCC");
        }
        assert_eq!(s.chars().filter(|&c| c != 'A').count(), 8 * here.len());
    }
    assert!(r.implants.iter().any(|x| x.minus_strand));
    assert!(r.implants.iter().any(|x| !x.minus_strand));

    // Fixed positions, forward strand only.
    let o = ImplantOptions {
        mode: ImplantMode::Positions(vec![vec![0, 49], vec![192], vec![300]]),
        both_strands: false,
        ..Default::default()
    };
    let r = implant_motifs(std::slice::from_ref(&m), &seqs[..3], &o).unwrap();
    assert_eq!(r.sequences[0][..8], *"GGGGGGGG");
    assert_eq!(r.sequences[0][49..57], *"GGGGGGGG");
    assert_eq!(r.sequences[1][192..], *"GGGGGGGG");
    assert_eq!(r.sequences[2], seqs[2]);
    assert_eq!(r.warnings.len(), 1);

    // A rate gives about rate * length implants per sequence.
    let o = ImplantOptions {
        mode: ImplantMode::Rate(0.01),
        seed: 11,
        ..Default::default()
    };
    let r = implant_motifs(std::slice::from_ref(&m), &seqs, &o).unwrap();
    let mean = r.implants.len() as f64 / 30.0;
    assert!((1.0..3.0).contains(&mean), "{mean}");
}
