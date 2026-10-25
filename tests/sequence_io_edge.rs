//! Edge cases and behaviour checks for the sequence utilities (sampling,
//! shuffling, complexity, k-mers, encoding, masking and deduplication) and for
//! the motif and FASTA readers and writers.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use universalmotif::complexity::{
    ComplexityMethod, ComplexityOptions, complexity, sequence_complexity,
};
use universalmotif::dedup::{Interval, dedup};
use universalmotif::io::{cisbp, fasta, homer, jaspar, matrix, meme, transfac, uniprobe};
use universalmotif::kmer::{
    BkgOptions, WindowOptions, count_kmers, get_bkg, kmers, sequence_letters, windows,
};
use universalmotif::mask::{mask_ranges, mask_seqs};
use universalmotif::packed::EncodedSequences;
use universalmotif::sample::{ImplantMode, ImplantOptions, implant_motifs, sample_sites};
use universalmotif::shuffle::{
    ShuffleMethod, ShuffleOptions, create_sequences, shuffle_motifs, shuffle_sequences,
};
use universalmotif::{Alphabet, Error, Matrix, Motif, MotifType, Strand};

// Helpers.

static COUNTER: AtomicUsize = AtomicUsize::new(0);

/// A fresh path in the temporary directory, unique to this process and call.
fn temp_path(name: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir().join(format!("um-seqio-edge-{}-{n}-{name}", std::process::id()))
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    e.write_all(bytes).unwrap();
    e.finish().unwrap()
}

fn arg_msg<T: std::fmt::Debug>(r: Result<T, Error>) -> String {
    match r {
        Err(Error::Argument(m)) => m,
        other => panic!("expected an argument error, got {other:?}"),
    }
}

fn unsupported_msg<T: std::fmt::Debug>(r: Result<T, Error>) -> String {
    match r {
        Err(Error::Unsupported(m)) => m,
        other => panic!("expected an unsupported error, got {other:?}"),
    }
}

fn motif(alphabet: Alphabet, cols: &[&[f64]], motif_type: MotifType, name: &str) -> Motif {
    Motif::builder(Matrix::from_columns(cols).unwrap(), alphabet)
        .name(name)
        .motif_type(Some(motif_type))
        .pseudocount(0.0)
        .build()
        .unwrap()
}

fn ppm(alphabet: Alphabet, cols: &[&[f64]], name: &str) -> Motif {
    motif(alphabet, cols, MotifType::Ppm, name)
}

/// A motif that always gives the same site.
fn one_hot(site: &str) -> Motif {
    let cols: Vec<Vec<f64>> = site
        .chars()
        .map(|c| {
            let mut v = vec![0.0; 4];
            v["ACGT".find(c).unwrap()] = 1.0;
            v
        })
        .collect();
    let cols: Vec<&[f64]> = cols.iter().map(Vec::as_slice).collect();
    ppm(Alphabet::Dna, &cols, site)
}

fn assert_matrix_close(a: &Matrix, b: &Matrix, tol: f64, ctx: &str) {
    assert_eq!((a.nrow(), a.ncol()), (b.nrow(), b.ncol()), "{ctx}: shape");
    for (x, y) in a.as_slice().iter().zip(b.as_slice()) {
        assert!((x - y).abs() <= tol, "{ctx}: {x} vs {y}");
    }
}

fn kmer_table(s: &str, k: usize) -> HashMap<Vec<char>, usize> {
    let chars: Vec<char> = s.chars().collect();
    let mut out = HashMap::new();
    for w in chars.windows(k) {
        *out.entry(w.to_vec()).or_insert(0) += 1;
    }
    out
}

fn sorted_letters(s: &str) -> Vec<char> {
    let mut v: Vec<char> = s.chars().collect();
    v.sort_unstable();
    v
}

/// A pseudo-random DNA string from a simple linear congruential generator.
fn random_dna(len: usize, seed: u64) -> String {
    let mut x = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 60) as usize & 3] as char
        })
        .collect()
}

// Sampling and implanting.

#[test]
fn sampled_sites_follow_the_motif_probabilities() {
    let m = ppm(
        Alphabet::Dna,
        &[
            &[0.7, 0.1, 0.1, 0.1],
            &[0.0, 0.0, 0.0, 1.0],
            &[0.25, 0.25, 0.25, 0.25],
        ],
        "m",
    );
    let n = 20000;
    let sites = sample_sites(&m, n, 7).unwrap();
    assert_eq!(sites.len(), n);
    assert!(sites.iter().all(|s| s.len() == 3));
    let freq = |pos: usize, letter: u8| {
        sites.iter().filter(|s| s.as_bytes()[pos] == letter).count() as f64 / n as f64
    };
    // Five standard deviations of a binomial proportion is below 0.02 here.
    for (letter, want) in [(b'A', 0.7), (b'C', 0.1), (b'G', 0.1), (b'T', 0.1)] {
        assert!(
            (freq(0, letter) - want).abs() < 0.02,
            "position 1, {letter}"
        );
    }
    assert_eq!(freq(1, b'T'), 1.0, "a probability of 1 is always drawn");
    for letter in *b"ACGT" {
        assert!(
            (freq(2, letter) - 0.25).abs() < 0.02,
            "position 3, {letter}"
        );
    }
    assert_eq!(sample_sites(&m, n, 7).unwrap(), sites, "same seed");
    assert_ne!(sample_sites(&m, 50, 8).unwrap(), sites[..50], "other seed");
    assert!(sample_sites(&m, 0, 1).unwrap().is_empty());
}

#[test]
fn sampled_sites_use_any_alphabet_and_type() {
    let custom = Alphabet::custom("xyz".chars()).unwrap();
    let m = ppm(custom, &[&[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0]], "c");
    assert_eq!(sample_sites(&m, 3, 1).unwrap(), ["yz", "yz", "yz"]);
    // A count matrix with no pseudocount samples its exact proportions.
    let pcm = motif(
        Alphabet::Dna,
        &[&[0.0, 0.0, 10.0, 0.0]],
        MotifType::Pcm,
        "p",
    );
    assert_eq!(sample_sites(&pcm, 4, 1).unwrap(), ["G"; 4]);
    let rna = ppm(Alphabet::Rna, &[&[0.0, 0.0, 0.0, 1.0]], "r");
    assert_eq!(sample_sites(&rna, 1, 1).unwrap(), ["U"]);
}

#[test]
fn implants_overwrite_only_their_sites() {
    let m = one_hot("AACG");
    let seqs: Vec<String> = (0..60).map(|_| "N".repeat(30)).collect();
    let out = implant_motifs(std::slice::from_ref(&m), &seqs, &ImplantOptions::default()).unwrap();
    assert!(out.warnings.is_empty(), "{:?}", out.warnings);
    assert_eq!(out.implants.len(), seqs.len());
    let (mut plus, mut minus) = (0, 0);
    for imp in &out.implants {
        let s = &out.sequences[imp.sequence];
        assert_eq!(s.len(), 30);
        assert_eq!(imp.width, 4);
        assert_eq!(&s[imp.start..imp.start + 4], imp.planted);
        // Everything else is untouched.
        assert_eq!(s.matches('N').count(), 26);
        if imp.minus_strand {
            assert_eq!(imp.planted, "CGTT");
            minus += 1;
        } else {
            assert_eq!(imp.planted, "AACG");
            plus += 1;
        }
    }
    assert!(
        plus > 10 && minus > 10,
        "both strands are used ({plus}, {minus})"
    );

    let opts = ImplantOptions {
        both_strands: false,
        ..Default::default()
    };
    let out = implant_motifs(&[m], &seqs, &opts).unwrap();
    assert!(out.implants.iter().all(|i| !i.minus_strand));
    assert!(out.sequences.iter().all(|s| s.contains("AACG")));
}

#[test]
fn implants_are_reproducible_by_seed() {
    let m = ppm(Alphabet::Dna, &[&[0.25; 4], &[0.25; 4], &[0.25; 4]], "u");
    let seqs = vec!["ACGTACGTACGTACGTACGT"; 8];
    let opts = ImplantOptions {
        mode: ImplantMode::Count(2),
        ..Default::default()
    };
    let a = implant_motifs(std::slice::from_ref(&m), &seqs, &opts).unwrap();
    let b = implant_motifs(std::slice::from_ref(&m), &seqs, &opts).unwrap();
    assert_eq!(a, b);
    let c = implant_motifs(&[m], &seqs, &ImplantOptions { seed: 2, ..opts }).unwrap();
    assert_ne!(a.implants, c.implants);
}

#[test]
fn implants_at_given_positions() {
    let m = one_hot("GGGG");
    let seqs = ["AAAAAAAAAAAA", "TTTTTTTTTTTT"];
    let opts = ImplantOptions {
        mode: ImplantMode::Positions(vec![vec![0, 10, 2], vec![8]]),
        both_strands: false,
        ..Default::default()
    };
    let out = implant_motifs(std::slice::from_ref(&m), &seqs, &opts).unwrap();
    assert_eq!(out.sequences, ["GGGGGGAAAAAA", "TTTTTTTTGGGG"]);
    let starts: Vec<(usize, usize)> = out.implants.iter().map(|i| (i.sequence, i.start)).collect();
    assert_eq!(starts, [(0, 0), (0, 2), (1, 8)]);
    assert_eq!(
        out.warnings,
        [
            "dropping out-of-bounds position for sequence 1 (start=10, width=4, sequence length=12)",
            "position 2 in sequence 1 overlaps a previous implant; planting anyway",
        ]
    );
    let bad = ImplantOptions {
        mode: ImplantMode::Positions(vec![vec![0]]),
        ..Default::default()
    };
    assert_eq!(
        arg_msg(implant_motifs(&[m], &seqs, &bad)),
        "give one list of positions per sequence"
    );
}

#[test]
fn implant_rate_gives_poisson_counts() {
    let m = one_hot("ACGT");
    let seqs = vec!["N".repeat(1000); 40];
    let opts = ImplantOptions {
        mode: ImplantMode::Rate(0.005),
        ..Default::default()
    };
    let out = implant_motifs(std::slice::from_ref(&m), &seqs, &opts).unwrap();
    // A mean of 5 per sequence, so 200 in all with a standard deviation of
    // about 14.
    let n = out.implants.len();
    assert!((140..=260).contains(&n), "{n} implants");
    let per_seq: Vec<usize> = (0..seqs.len())
        .map(|i| out.implants.iter().filter(|x| x.sequence == i).count())
        .collect();
    assert!(per_seq.iter().any(|&c| c != per_seq[0]), "counts vary");

    let zero = ImplantOptions {
        mode: ImplantMode::Rate(0.0),
        ..Default::default()
    };
    let out = implant_motifs(std::slice::from_ref(&m), &seqs, &zero).unwrap();
    assert!(out.implants.is_empty());
    assert_eq!(out.sequences, seqs);
    for r in [-1.0, f64::NAN] {
        let bad = ImplantOptions {
            mode: ImplantMode::Rate(r),
            ..Default::default()
        };
        assert_eq!(
            arg_msg(implant_motifs(std::slice::from_ref(&m), &seqs, &bad)),
            "rate must be non-negative"
        );
    }
}

#[test]
fn implants_keep_their_spacing() {
    let m = one_hot("ACGT");
    let seqs = vec!["N".repeat(200); 30];
    let opts = ImplantOptions {
        mode: ImplantMode::Count(6),
        min_spacing: 5,
        ..Default::default()
    };
    let out = implant_motifs(&[m], &seqs, &opts).unwrap();
    assert!(out.warnings.is_empty(), "{:?}", out.warnings);
    for i in 0..seqs.len() {
        let mut spans: Vec<(usize, usize)> = out
            .implants
            .iter()
            .filter(|x| x.sequence == i)
            .map(|x| (x.start, x.start + x.width))
            .collect();
        assert_eq!(spans.len(), 6);
        spans.sort_unstable();
        for w in spans.windows(2) {
            assert!(w[1].0 >= w[0].1 + 5, "sequence {i}: {spans:?}");
        }
    }
}

#[test]
fn centre_bias_pulls_implants_to_the_middle() {
    let m = one_hot("ACGT");
    let seqs = vec!["N".repeat(200); 400];
    let mean_offset = |bias: usize| {
        let opts = ImplantOptions {
            centre_bias: bias,
            ..Default::default()
        };
        let out = implant_motifs(std::slice::from_ref(&m), &seqs, &opts).unwrap();
        out.implants
            .iter()
            .map(|x| (x.start as f64 - 98.0).abs())
            .sum::<f64>()
            / out.implants.len() as f64
    };
    // Uniform starts lie about 49 from the middle on average, and means of ten
    // uniform draws about 14.
    let (flat, biased) = (mean_offset(1), mean_offset(10));
    assert!(flat > 38.0, "{flat}");
    assert!(biased < 22.0, "{biased}");
}

#[test]
fn implants_warn_when_they_do_not_fit() {
    let m = one_hot("ACGT");
    let out = implant_motifs(
        std::slice::from_ref(&m),
        &["ACG", "NNNNNNNNNN"],
        &ImplantOptions::default(),
    )
    .unwrap();
    assert_eq!(out.sequences[0], "ACG");
    assert!(out.implants.iter().all(|x| x.sequence == 1));
    assert_eq!(
        out.warnings,
        ["some sequences are shorter than the narrowest motif; no implants will be placed in them"]
    );
    let crowded = ImplantOptions {
        mode: ImplantMode::Count(5),
        min_spacing: 3,
        max_retries: 20,
        ..Default::default()
    };
    let out = implant_motifs(&[m], &["NNNNNNNNNN"], &crowded).unwrap();
    assert!(out.implants.len() < 5);
    assert_eq!(out.warnings.len(), 1);
    assert!(
        out.warnings[0].starts_with("one or more implants could not be placed within max_retries")
    );
}

#[test]
fn implants_choose_among_motifs_and_handle_rna() {
    let a = one_hot("AAAA");
    let c = one_hot("CCCCCC");
    let seqs = vec!["T".repeat(40); 50];
    let opts = ImplantOptions {
        both_strands: false,
        ..Default::default()
    };
    let out = implant_motifs(&[a, c], &seqs, &opts).unwrap();
    for imp in &out.implants {
        let want = if imp.motif == 0 { "AAAA" } else { "CCCCCC" };
        assert_eq!(imp.planted, want);
        assert_eq!(imp.width, want.len());
    }
    let used: Vec<usize> = out.implants.iter().map(|x| x.motif).collect();
    assert!(used.contains(&0) && used.contains(&1));

    let rna = ppm(
        Alphabet::Rna,
        &[&[1.0, 0.0, 0.0, 0.0], &[0.0, 1.0, 0.0, 0.0]],
        "r",
    );
    let out = implant_motifs(&[rna], &vec!["GGGGGG"; 30], &ImplantOptions::default()).unwrap();
    for imp in &out.implants {
        let want = if imp.minus_strand { "GU" } else { "AC" };
        assert_eq!(imp.planted, want);
    }
}

#[test]
fn implant_argument_errors() {
    let dna = one_hot("ACGT");
    let seqs = ["ACGTACGT"];
    let none: [Motif; 0] = [];
    assert_eq!(
        arg_msg(implant_motifs(&none, &seqs, &ImplantOptions::default())),
        "no motifs given"
    );
    let rna = ppm(Alphabet::Rna, &[&[0.25; 4]], "r");
    assert_eq!(
        arg_msg(implant_motifs(
            &[dna.clone(), rna],
            &seqs,
            &ImplantOptions::default()
        )),
        "all motifs must share the same alphabet"
    );
    let aa = ppm(Alphabet::Aa, &[&[0.05; 20]], "aa");
    assert_eq!(
        unsupported_msg(implant_motifs(&[aa], &seqs, &ImplantOptions::default())),
        "only DNA and RNA motifs can be implanted"
    );
    for opts in [
        ImplantOptions {
            centre_bias: 0,
            ..Default::default()
        },
        ImplantOptions {
            max_retries: 0,
            ..Default::default()
        },
    ] {
        assert_eq!(
            arg_msg(implant_motifs(std::slice::from_ref(&dna), &seqs, &opts)),
            "centre_bias and max_retries must be at least 1"
        );
    }
}

// Shuffling.

#[test]
fn euler_shuffles_keep_kmer_counts_and_ends() {
    let seq = random_dna(300, 3);
    for k in 2..=4 {
        let opts = ShuffleOptions {
            k,
            seed: 11,
            ..Default::default()
        };
        let out = shuffle_sequences(&[seq.as_str()], &opts).unwrap();
        let s = &out[0];
        assert_ne!(s, &seq, "k = {k}");
        assert_eq!(kmer_table(s, k), kmer_table(&seq, k), "k = {k}");
        assert_eq!(s[..k - 1], seq[..k - 1], "first letters, k = {k}");
        assert_eq!(
            s[s.len() - (k - 1)..],
            seq[seq.len() - (k - 1)..],
            "last letters, k = {k}"
        );
    }
}

#[test]
fn single_letter_shuffles_keep_composition() {
    let seq = random_dna(100, 9);
    let out = shuffle_sequences(&[seq.as_str()], &ShuffleOptions::default()).unwrap();
    assert_ne!(out[0], seq);
    assert_eq!(sorted_letters(&out[0]), sorted_letters(&seq));
    // A single letter or a one-letter alphabet has nothing to move.
    let opts = ShuffleOptions {
        k: 2,
        ..Default::default()
    };
    assert_eq!(
        shuffle_sequences(&["G", "TTTTTT"], &ShuffleOptions::default()).unwrap(),
        ["G", "TTTTTT"]
    );
    assert_eq!(shuffle_sequences(&["AAAAAA"], &opts).unwrap(), ["AAAAAA"]);
    assert_eq!(shuffle_sequences(&["AC"], &opts).unwrap(), ["AC"]);
}

#[test]
fn linear_shuffles_move_whole_blocks() {
    let seq = "AACCGGTTACGTA";
    let opts = ShuffleOptions {
        k: 3,
        method: ShuffleMethod::Linear,
        seed: 4,
        ..Default::default()
    };
    let out = shuffle_sequences(&[seq], &opts).unwrap();
    let s = &out[0];
    let blocks = |x: &str| {
        let mut b: Vec<String> = x.as_bytes()[..12]
            .chunks(3)
            .map(|c| String::from_utf8(c.to_vec()).unwrap())
            .collect();
        b.sort();
        b
    };
    assert_eq!(blocks(s), blocks(seq));
    assert_eq!(&s[12..], "A", "the leftover letter stays at the end");
    assert_eq!(sorted_letters(s), sorted_letters(seq));
}

#[test]
fn markov_shuffles_only_use_seen_transitions() {
    let opts = ShuffleOptions {
        k: 2,
        method: ShuffleMethod::Markov,
        seed: 5,
        ..Default::default()
    };
    let out = shuffle_sequences(&["ACACACACAC", "AAAACCCCGGGGTTTT"], &opts).unwrap();
    assert_eq!(out[0].len(), 10);
    assert!(out[0] == "ACACACACAC" || out[0] == "CACACACACA");
    // Every two-letter word of the second comes from the original, except
    // where the walk reaches T, which has no successor of its own.
    let seen = kmer_table("AAAACCCCGGGGTTTT", 2);
    let chars: Vec<char> = out[1].chars().collect();
    assert_eq!(chars.len(), 16);
    for w in chars.windows(2) {
        assert!(w[0] == 'T' || seen.contains_key(w), "{:?}", out[1]);
    }
    let long = random_dna(4000, 1);
    let opts = ShuffleOptions {
        k: 3,
        method: ShuffleMethod::Markov,
        ..Default::default()
    };
    let shuffled = &shuffle_sequences(&[long.as_str()], &opts).unwrap()[0];
    let (a, b) = (kmer_table(&long, 3), kmer_table(shuffled, 3));
    // Counts are kept on average, so each 3-mer, seen about 62 times, stays
    // close.
    for (kmer, &n) in &a {
        let m = *b.get(kmer).unwrap_or(&0) as f64;
        assert!((m - n as f64).abs() < 40.0, "{kmer:?}: {n} vs {m}");
    }
}

#[test]
fn shuffles_are_seeded_and_thread_independent() {
    let seqs: Vec<String> = (0..12).map(|i| random_dna(80, i)).collect();
    for method in [
        ShuffleMethod::Euler,
        ShuffleMethod::Markov,
        ShuffleMethod::Linear,
    ] {
        let opts = ShuffleOptions {
            k: 2,
            method,
            seed: 99,
            nthreads: 1,
            ..Default::default()
        };
        let one = shuffle_sequences(&seqs, &opts).unwrap();
        let many = shuffle_sequences(
            &seqs,
            &ShuffleOptions {
                nthreads: 3,
                ..opts.clone()
            },
        )
        .unwrap();
        assert_eq!(one, many, "{method}");
        let other = shuffle_sequences(&seqs, &ShuffleOptions { seed: 100, ..opts }).unwrap();
        assert_ne!(one, other, "{method}");
        // Each sequence gets its own generator, so it does not matter what
        // comes after it.
        let first = shuffle_sequences(
            &seqs[..1],
            &ShuffleOptions {
                seed: 99,
                k: 2,
                method,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(first[0], one[0], "{method}");
    }
}

#[test]
fn non_ascii_sequences_shuffle_like_ascii_ones() {
    let ascii = random_dna(60, 5);
    // A map that keeps the order of the letters, so the draws line up.
    let map = |c: char| match c {
        'A' => 'α',
        'C' => 'β',
        'G' => 'γ',
        _ => 'δ',
    };
    let greek: String = ascii.chars().map(map).collect();
    for method in [
        ShuffleMethod::Euler,
        ShuffleMethod::Markov,
        ShuffleMethod::Linear,
    ] {
        let opts = ShuffleOptions {
            k: 2,
            method,
            ..Default::default()
        };
        let a = shuffle_sequences(&[ascii.as_str()], &opts).unwrap();
        let g = shuffle_sequences(&[greek.as_str()], &opts).unwrap();
        assert_eq!(g[0], a[0].chars().map(map).collect::<String>(), "{method}");
    }
}

#[test]
fn windowed_shuffles_stay_inside_their_windows() {
    let seq = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTT";
    let opts = ShuffleOptions {
        window: Some(WindowOptions {
            size: 10.0,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    // Every window holds one letter, so nothing can change.
    assert_eq!(shuffle_sequences(&[seq], &opts).unwrap(), [seq]);
    let mixed = random_dna(40, 2);
    let out = &shuffle_sequences(&[mixed.as_str()], &opts).unwrap()[0];
    assert_ne!(out, &mixed);
    for start in (0..40).step_by(10) {
        assert_eq!(
            sorted_letters(&out[start..start + 10]),
            sorted_letters(&mixed[start..start + 10])
        );
    }
    // A fractional window is a share of each sequence's length.
    let frac = ShuffleOptions {
        k: 2,
        window: Some(WindowOptions {
            size: 0.25,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    let out = &shuffle_sequences(&[mixed.as_str()], &frac).unwrap()[0];
    for start in (0..40).step_by(10) {
        assert_eq!(
            kmer_table(&out[start..start + 10], 2),
            kmer_table(&mixed[start..start + 10], 2)
        );
    }
}

#[test]
fn markov_short_last_window_keeps_its_letters() {
    // Windows of 4 over 10 letters leave a last window of 2 ("GT"), shorter
    // than k = 3; it should be left as it is.
    let opts = ShuffleOptions {
        k: 3,
        method: ShuffleMethod::Markov,
        window: Some(WindowOptions {
            size: 4.0,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    let out = shuffle_sequences(&["ACGTACGTGT"], &opts).unwrap();
    assert_eq!(&out[0][8..], "GT");
}

#[test]
fn shuffle_argument_errors() {
    let seqs = ["ACGTACGT", "ACG"];
    let opts = |k, window| ShuffleOptions {
        k,
        window,
        ..Default::default()
    };
    assert_eq!(
        arg_msg(shuffle_sequences(&seqs, &opts(0, None))),
        "k must be at least 1"
    );
    assert_eq!(
        arg_msg(shuffle_sequences(&seqs, &opts(4, None))),
        "k (4) must not exceed the shortest sequence length (3)"
    );
    let w = |size, overlap| Some(WindowOptions { size, overlap });
    for bad in [w(0.0, 0.0), w(f64::NAN, 0.0), w(4.0, -1.0)] {
        assert_eq!(
            arg_msg(shuffle_sequences(&seqs, &opts(1, bad))),
            "the window size must be positive and the overlap non-negative"
        );
    }
    assert_eq!(
        arg_msg(shuffle_sequences(&seqs, &opts(1, w(3.0, 3.0)))),
        "the window overlap must be smaller than the window size"
    );
    assert_eq!(
        arg_msg(shuffle_sequences(&seqs, &opts(3, w(2.0, 0.0)))),
        "the window size must be at least k"
    );
}

#[test]
fn shuffle_methods_parse_and_print() {
    for m in [
        ShuffleMethod::Euler,
        ShuffleMethod::Markov,
        ShuffleMethod::Linear,
    ] {
        assert_eq!(m.to_string().parse::<ShuffleMethod>().unwrap(), m);
        assert_eq!(m.as_str(), m.to_string());
    }
    assert_eq!(ShuffleMethod::default(), ShuffleMethod::Euler);
    assert_eq!(
        arg_msg("Euler".parse::<ShuffleMethod>()),
        "unknown shuffling method 'Euler' (use euler, markov or linear)"
    );
}

#[test]
fn created_sequences_follow_their_frequencies() {
    let alph = ['A', 'C', 'G', 'T'];
    let seqs = create_sequences(&alph, 20, 500, None, 3, 1).unwrap();
    assert_eq!(seqs.len(), 20);
    assert!(seqs.iter().all(|s| s.len() == 500));
    let all: String = seqs.concat();
    for c in alph {
        let f = all.matches(c).count() as f64 / all.len() as f64;
        assert!((f - 0.25).abs() < 0.02, "{c}: {f}");
    }
    assert_eq!(create_sequences(&alph, 20, 500, None, 3, 4).unwrap(), seqs);
    assert_ne!(create_sequences(&alph, 20, 500, None, 4, 1).unwrap(), seqs);

    let skewed = [0.0, 0.9, 0.1, 0.0];
    let s = create_sequences(&alph, 1, 5000, Some(&skewed), 1, 1).unwrap();
    assert!(!s[0].contains('A') && !s[0].contains('T'));
    let c = s[0].matches('C').count() as f64 / 5000.0;
    assert!((c - 0.9).abs() < 0.03, "{c}");

    // Dinucleotides: only AC and CA, so the letters alternate.
    let mut di = [0.0; 16];
    di[1] = 0.5;
    di[4] = 0.5;
    let s = create_sequences(&alph, 10, 21, Some(&di), 1, 1).unwrap();
    for x in &s {
        assert!(x.starts_with("AC") || x.starts_with("CA"), "{x}");
        assert!(
            !x.contains("AA") && !x.contains("CC") && x.len() == 21,
            "{x}"
        );
    }
    // Sequences shorter than k still come out at the length asked for.
    // A context with no way forward falls back to the first letter.
    let mut dead_end = [0.0; 16];
    dead_end[1] = 1.0;
    let s = create_sequences(&alph, 1, 6, Some(&dead_end), 1, 1).unwrap();
    assert_eq!(s, ["ACACAC"]);
    let short = create_sequences(&alph, 3, 1, Some(&di), 1, 1).unwrap();
    assert!(short.iter().all(|x| x == "A" || x == "C"));
    assert_eq!(create_sequences(&alph, 0, 10, None, 1, 1).unwrap().len(), 0);
}

#[test]
fn created_sequences_check_their_arguments() {
    let alph = ['A', 'C'];
    assert_eq!(
        arg_msg(create_sequences(&[], 1, 5, None, 1, 1)),
        "the alphabet needs at least one letter"
    );
    assert_eq!(
        arg_msg(create_sequences(
            &alph,
            1,
            5,
            Some(&[0.5, 0.25, 0.25]),
            1,
            1
        )),
        "the number of frequencies must be a power of the alphabet size"
    );
    for bad in [[0.5, -0.5], [f64::NAN, 1.0], [f64::INFINITY, 1.0]] {
        assert_eq!(
            arg_msg(create_sequences(&alph, 1, 5, Some(&bad), 1, 1)),
            "frequencies must be finite and non-negative"
        );
    }
    assert_eq!(
        arg_msg(create_sequences(&alph, 1, 5, Some(&[0.0, 0.0]), 1, 1)),
        "at least one frequency must be above zero"
    );
}

#[test]
fn motif_shuffles_permute_columns() {
    let cols: Vec<Vec<f64>> = (0..5)
        .map(|i| {
            let mut v = [0.1; 4];
            v[i % 4] = 0.7;
            v[(i + 1) % 4] = 0.1 + 0.01 * i as f64;
            let s: f64 = v.iter().sum();
            v.iter().map(|x| x / s).collect()
        })
        .collect();
    let a = ppm(Alphabet::Dna, &[&cols[0], &cols[1], &cols[2]], "first");
    let pcm = motif(
        Alphabet::Dna,
        &[&[10.0, 0.0, 0.0, 0.0], &[0.0, 0.0, 0.0, 10.0]],
        MotifType::Pcm,
        "second",
    );
    let out = shuffle_motifs(&[a.clone(), pcm.clone()], 1, 3).unwrap();
    assert_eq!(out[0].name, "first [shuffled]");
    assert_eq!(out[1].name, "second [shuffled]");
    assert_eq!((out[0].width(), out[1].width()), (3, 2));
    assert_eq!(out[0].motif_type, MotifType::Ppm);
    assert_eq!(out[1].motif_type, MotifType::Pcm);
    // The columns, as probabilities, are those of the inputs in a new order.
    let as_ppm = |m: &Motif| {
        universalmotif::convert::convert_type(m, MotifType::Ppm, None, false)
            .unwrap()
            .matrix
            .columns()
            .map(|c| {
                c.iter()
                    .map(|x| (x * 1e6).round() as i64)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let mut before: Vec<_> = [&a, &pcm].iter().flat_map(|m| as_ppm(m)).collect();
    let mut after: Vec<_> = out.iter().flat_map(as_ppm).collect();
    before.sort();
    after.sort();
    assert_eq!(before, after);
    assert_eq!(
        shuffle_motifs(&[a.clone(), pcm.clone()], 1, 3).unwrap(),
        out
    );
}

#[test]
fn motif_block_shuffles_keep_the_ends() {
    // Eight columns in blocks of three leave one column at each end.
    let cols: Vec<Vec<f64>> = (0..8)
        .map(|i| {
            let p = 0.3 + 0.05 * i as f64;
            vec![p, (1.0 - p) / 3.0, (1.0 - p) / 3.0, (1.0 - p) / 3.0]
        })
        .collect();
    let refs: Vec<&[f64]> = cols.iter().map(Vec::as_slice).collect();
    let m = ppm(Alphabet::Dna, &refs, "m");
    for seed in 1..6 {
        let out = &shuffle_motifs(std::slice::from_ref(&m), 3, seed).unwrap()[0];
        let firsts: Vec<f64> = out.matrix.columns().map(|c| c[0]).collect();
        assert!((firsts[0] - 0.3).abs() < 1e-12);
        assert!((firsts[7] - 0.65).abs() < 1e-12);
        // The middle six move as two blocks of three consecutive columns.
        for b in [1, 4] {
            let step = firsts[b + 1] - firsts[b];
            assert!((step - 0.05).abs() < 1e-9, "seed {seed}: {firsts:?}");
            assert!((firsts[b + 2] - firsts[b + 1] - 0.05).abs() < 1e-9);
        }
    }
    // Two columns in blocks of three: the remainder of 2 splits one to each end.
    // Seven columns in blocks of three: the one left over stays first.
    let seven = ppm(Alphabet::Dna, &refs[..7], "seven");
    for seed in 1..4 {
        let out = &shuffle_motifs(std::slice::from_ref(&seven), 3, seed).unwrap()[0];
        assert_eq!(out.matrix.col(0), seven.matrix.col(0));
    }
    let two = ppm(Alphabet::Dna, &refs[..2], "two");
    assert_eq!(
        shuffle_motifs(std::slice::from_ref(&two), 3, 1).unwrap()[0].matrix,
        two.matrix
    );
}

#[test]
fn motif_shuffle_errors() {
    let none: [Motif; 0] = [];
    assert_eq!(arg_msg(shuffle_motifs(&none, 1, 1)), "no motifs given");
    let dna = ppm(Alphabet::Dna, &[&[0.25; 4]], "d");
    let rna = ppm(Alphabet::Rna, &[&[0.25; 4]], "r");
    assert_eq!(
        arg_msg(shuffle_motifs(&[dna.clone(), rna], 1, 1)),
        "all motifs must share the same alphabet"
    );
    assert_eq!(
        arg_msg(shuffle_motifs(&[dna], 0, 1)),
        "k must be at least 1"
    );
}

// Complexity.

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-12
}

#[test]
fn complexity_by_hand() {
    use ComplexityMethod::*;
    // Wootton-Federhen: log_n(L! / prod n_i!) / L, here log2(4! / (2! 2!)) / 4.
    assert!(approx(
        complexity("AACC", WoottonFederhen, None, 7),
        6f64.log2() / 4.0
    ));
    assert!(approx(
        complexity("AACC", WoottonFederhen, Some("ACGT"), 7),
        6f64.ln() / 4f64.ln() / 4.0
    ));
    // Its entropy approximation: two letters, each half the sequence.
    assert!(approx(
        complexity("AACC", WoottonFederhenFast, None, 7),
        1.0
    ));
    assert!(approx(
        complexity("AACC", WoottonFederhenFast, Some("ACGT"), 7),
        0.5
    ));
    // Trifonov over a two-letter alphabet: words of 1 to 4 letters, seen
    // 1, 1, 1, 1 of a possible 2, 3, 2, 1.
    assert!(approx(
        complexity("AAAA", Trifonov, Some("AC"), 7),
        1.0 / 12.0
    ));
    assert!(approx(
        complexity("AAAA", Trifonov, Some("AC"), 2),
        1.0 / 6.0
    ));
    assert!(approx(
        complexity("AAAA", TrifonovFast, Some("AC"), 7),
        4.0 / 8.0
    ));
    assert!(approx(
        complexity("AAAA", TrifonovFast, Some("AC"), 2),
        2.0 / 5.0
    ));
    assert!(approx(complexity("ACGT", Trifonov, None, 7), 1.0));
    // An empty alphabet string means the sequence's own letters.
    assert_eq!(
        complexity("AACG", Trifonov, Some(""), 7),
        complexity("AACG", Trifonov, None, 7)
    );
    // DUST: four AAA triplets give 4 * 3 / 2 over (6 - 2) - 1.
    assert!(approx(complexity("AAAAAA", Dust, None, 7), 2.0));
    // ACG and CGT twice each: (1 + 1) / (6 - 1).
    assert!(approx(complexity("ACGTACGT", Dust, None, 7), 0.4));
    assert!(approx(complexity("ACGTTGCA", Dust, None, 7), 0.0));
}

#[test]
fn complexity_methods_parse_and_print() {
    use ComplexityMethod::*;
    for m in [
        WoottonFederhen,
        WoottonFederhenFast,
        Trifonov,
        TrifonovFast,
        Dust,
    ] {
        assert_eq!(m.to_string().parse::<ComplexityMethod>().unwrap(), m);
        assert_eq!(m.as_str(), m.to_string());
    }
    assert_eq!(Dust.as_str(), "DUST");
    assert_eq!(ComplexityMethod::default(), WoottonFederhen);
    assert!(
        arg_msg("dust".parse::<ComplexityMethod>()).starts_with("unknown complexity method 'dust'")
    );
}

#[test]
fn windowed_complexity_matches_whole_windows() {
    let seq = random_dna(40, 8);
    for method in [
        ComplexityMethod::WoottonFederhen,
        ComplexityMethod::WoottonFederhenFast,
        ComplexityMethod::Trifonov,
        ComplexityMethod::TrifonovFast,
        ComplexityMethod::Dust,
    ] {
        let opts = ComplexityOptions {
            method,
            ..Default::default()
        };
        let w = sequence_complexity(&[seq.as_str()], &opts).unwrap();
        // Windows of 20 overlapping by 10.
        let spans: Vec<(usize, usize)> = w.iter().map(|x| (x.start, x.end)).collect();
        assert_eq!(spans, [(0, 20), (10, 30), (20, 40)], "{method}");
        let alph = sequence_letters(&[seq.as_str()])
            .into_iter()
            .collect::<String>();
        for x in &w {
            assert_eq!(x.sequence, 0);
            let want = complexity(&seq[x.start..x.end], method, Some(&alph), 7);
            assert_eq!(x.complexity.to_bits(), want.to_bits(), "{method}");
        }
        let threaded = sequence_complexity(
            &[seq.as_str()],
            &ComplexityOptions {
                nthreads: 2,
                ..opts
            },
        )
        .unwrap();
        assert_eq!(threaded, w);
    }
}

#[test]
fn windowed_complexity_options() {
    let seqs = ["ACGTACGTAACCGGTTACGTACGTAACCGGTTACGTACGT", "ACGTAC"];
    // A fractional window is a share of each sequence, with no overlap by
    // default when half of it rounds to 0.
    let opts = ComplexityOptions {
        window: 0.5,
        ..Default::default()
    };
    let w = sequence_complexity(&seqs, &opts).unwrap();
    let spans: Vec<(usize, usize, usize)> =
        w.iter().map(|x| (x.sequence, x.start, x.end)).collect();
    assert_eq!(spans, [(0, 0, 20), (0, 20, 40), (1, 0, 3), (1, 3, 6)]);
    // Windows longer than a sequence shrink to it.
    let w = sequence_complexity(&seqs, &ComplexityOptions::default()).unwrap();
    assert_eq!((w[3].sequence, w[3].start, w[3].end), (1, 0, 6));
    // A given overlap, and a fractional one.
    let opts = ComplexityOptions {
        window: 10.0,
        overlap: Some(5.0),
        ..Default::default()
    };
    let w = sequence_complexity(&seqs[..1], &opts).unwrap();
    assert_eq!(
        w.iter().map(|x| x.start).collect::<Vec<_>>(),
        [0, 5, 10, 15, 20, 25, 30]
    );
    let opts = ComplexityOptions {
        window: 10.0,
        overlap: Some(0.1),
        ..Default::default()
    };
    let w = sequence_complexity(&seqs[..1], &opts).unwrap();
    assert_eq!(w[1].start, 6, "an overlap of a tenth of 40 letters is 4");
    // Extra alphabet letters count towards the possible letters.
    let wide = ComplexityOptions {
        method: ComplexityMethod::WoottonFederhenFast,
        alphabet: Some("ACGTN".into()),
        window: 6.0,
        ..Default::default()
    };
    let w = sequence_complexity(&["AACCAA"], &wide).unwrap();
    assert!(approx(
        w[0].complexity,
        complexity(
            "AACCAA",
            ComplexityMethod::WoottonFederhenFast,
            Some("ACGNT"),
            7
        )
    ));
    let single = ComplexityOptions {
        alphabet: Some("C".into()),
        ..Default::default()
    };
    assert_eq!(sequence_complexity(&["AAAA"], &single).unwrap().len(), 1);
}

#[test]
fn complexity_counts_letters_not_bytes() {
    // Two letters, once each, over an alphabet of those two letters.
    assert!(approx(
        complexity("αβ", ComplexityMethod::WoottonFederhenFast, None, 7),
        1.0
    ));
    let opts = ComplexityOptions {
        method: ComplexityMethod::WoottonFederhenFast,
        window: 3.0,
        overlap: Some(0.0),
        ..Default::default()
    };
    let w = sequence_complexity(&["αβγαβγ"], &opts).unwrap();
    assert_eq!(w.len(), 2);
}

#[test]
fn windowed_complexity_errors() {
    assert_eq!(
        arg_msg(sequence_complexity(
            &["AAAA"],
            &ComplexityOptions::default()
        )),
        "complexity cannot be calculated with a single-letter alphabet (A)"
    );
    let zero = ComplexityOptions {
        window: 0.0,
        ..Default::default()
    };
    assert_eq!(
        arg_msg(sequence_complexity(&["ACGT"], &zero)),
        "the window size must be greater than 0"
    );
    let overlap = ComplexityOptions {
        window: 4.0,
        overlap: Some(4.0),
        ..Default::default()
    };
    assert_eq!(
        arg_msg(sequence_complexity(&["ACGTACGT"], &overlap)),
        "the window overlap must be smaller than the window size"
    );
}

// k-mers and backgrounds.

#[test]
fn kmers_are_listed_in_order() {
    assert_eq!(kmers(&['A', 'C'], 0), [""]);
    assert_eq!(kmers(&['A', 'C'], 2), ["AA", "AC", "CA", "CC"]);
    assert_eq!(kmers(&['T', 'A'], 2), ["TT", "TA", "AT", "AA"]);
    assert_eq!(kmers(&['A', 'C', 'G', 'T'], 3).len(), 64);
    assert_eq!(kmers(&['A', 'C', 'G', 'T'], 3)[27], "CGT");
    assert_eq!(
        sequence_letters(&["TGCA", "βα", "aN", "Aβ"]),
        ['A', 'C', 'G', 'N', 'T', 'a', 'α', 'β']
    );
}

#[test]
fn kmer_counts_skip_unknown_letters() {
    let dna = ['A', 'C', 'G', 'T'];
    // N and the lower-case letters break the runs.
    let c = count_kmers("ACNGTacGT", 2, &dna);
    let names = kmers(&dna, 2);
    let got: Vec<(&str, u64)> = names
        .iter()
        .zip(&c)
        .filter(|(_, n)| **n > 0)
        .map(|(k, n)| (k.as_str(), *n))
        .collect();
    assert_eq!(got, [("AC", 1), ("GT", 2)]);
    assert_eq!(count_kmers("ACGT", 0, &dna), [0]);
    assert_eq!(count_kmers("ACG", 4, &dna).iter().sum::<u64>(), 0);
    // A non-ASCII alphabet with an ASCII sequence, and the other way round.
    assert_eq!(count_kmers("ABAB", 1, &['A', 'B', 'é']), [2, 2, 0]);
    assert_eq!(count_kmers("AéAé", 2, &['A', 'é']), [0, 2, 1, 0]);
    // Every k-mer of a long sequence is counted once per occurrence.
    let seq = random_dna(500, 4);
    for k in 1..=4 {
        let counts = count_kmers(&seq, k, &dna);
        assert_eq!(counts.iter().sum::<u64>() as usize, 500 - k + 1);
        let table = kmer_table(&seq, k);
        for (kmer, n) in kmers(&dna, k).iter().zip(&counts) {
            let key: Vec<char> = kmer.chars().collect();
            assert_eq!(*table.get(&key).unwrap_or(&0) as u64, *n, "{kmer}");
        }
    }
}

#[test]
fn sliding_windows() {
    assert_eq!(windows(10, 4, 0, true), [(0, 4), (4, 8), (8, 10)]);
    assert_eq!(windows(10, 4, 0, false), [(0, 4), (4, 8)]);
    assert_eq!(windows(10, 4, 2, true), [(0, 4), (2, 6), (4, 8), (6, 10)]);
    assert_eq!(windows(8, 4, 0, false), [(0, 4), (4, 8)]);
    assert_eq!(windows(4, 4, 1, true), [(0, 4)]);
    assert_eq!(windows(3, 1, 0, true), [(0, 1), (1, 2), (2, 3)]);
    assert_eq!(windows(6, 4, 3, false), [(0, 4), (1, 5), (2, 6)]);
    assert!(windows(3, 4, 0, true).is_empty());
    assert!(windows(10, 4, 4, true).is_empty());
    assert!(windows(10, 0, 0, true).is_empty());
}

fn bkg_counts(rows: &[universalmotif::kmer::BkgRow]) -> Vec<(String, u64)> {
    rows.iter().map(|r| (r.kmer.clone(), r.count)).collect()
}

fn pairs(v: &[(&str, u64)]) -> Vec<(String, u64)> {
    v.iter().map(|(k, n)| (k.to_string(), *n)).collect()
}

#[test]
fn backgrounds_merged_and_per_sequence() {
    let seqs = ["ACGT", "AACC"];
    let opts = BkgOptions {
        k: vec![2, 1, 1],
        ..Default::default()
    };
    let rows = get_bkg(&seqs, &[], &opts).unwrap();
    assert_eq!(rows.len(), 4 + 16, "k is sorted and deduplicated");
    assert_eq!(
        bkg_counts(&rows[..4]),
        pairs(&[("A", 3), ("C", 3), ("G", 1), ("T", 1)])
    );
    assert!(
        rows.iter()
            .all(|r| r.sequence.is_none() && r.start.is_none())
    );
    let p: f64 = rows[4..].iter().map(|r| r.probability).sum();
    assert!(approx(p, 1.0));
    assert!(
        approx(rows[4 + 1].probability, 2.0 / 6.0),
        "AC: {:?}",
        rows[5]
    );

    let opts = BkgOptions {
        k: vec![1],
        merge: false,
        pseudocount: 1,
        ..Default::default()
    };
    let rows = get_bkg(&seqs, &[], &opts).unwrap();
    // Ordered by k-mer, then sequence.
    let got: Vec<(usize, &str, u64)> = rows
        .iter()
        .map(|r| (r.sequence.unwrap(), r.kmer.as_str(), r.count))
        .collect();
    assert_eq!(
        got,
        [
            (0, "A", 2),
            (1, "A", 3),
            (0, "C", 2),
            (1, "C", 3),
            (0, "G", 2),
            (1, "G", 1),
            (0, "T", 2),
            (1, "T", 1)
        ]
    );
    assert!(approx(rows[1].probability, 3.0 / 8.0));
}

#[test]
fn backgrounds_with_a_given_alphabet() {
    let opts = BkgOptions {
        k: vec![1],
        alphabet: Some(vec!['T', 'A', 'A']),
        ..Default::default()
    };
    let rows = get_bkg(&["AATNC"], &[], &opts).unwrap();
    assert_eq!(bkg_counts(&rows), pairs(&[("A", 2), ("T", 1)]));
    // Nothing counted gives zero probabilities rather than NaN.
    let rows = get_bkg(&["GGGG"], &[], &opts).unwrap();
    assert!(rows.iter().all(|r| r.count == 0 && r.probability == 0.0));
}

#[test]
fn backgrounds_with_reverse_complements() {
    let opts = BkgOptions {
        k: vec![2],
        rc: true,
        ..Default::default()
    };
    let rows = get_bkg(&["AACGTT", "GATC"], &[], &opts).unwrap();
    let count = |k: &str| rows.iter().find(|r| r.kmer == k).unwrap().count;
    // Each k-mer is counted as often as its reverse complement.
    for kmer in kmers(&['A', 'C', 'G', 'T'], 2) {
        let rc: String = kmer
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        assert_eq!(count(&kmer), count(&rc), "{kmer}");
    }
    assert_eq!(count("AA"), 2);
    let rna = get_bkg(&["AACGUU"], &[], &BkgOptions { k: vec![1], ..opts }).unwrap();
    assert_eq!(
        bkg_counts(&rna),
        pairs(&[("A", 4), ("C", 2), ("G", 2), ("U", 4)])
    );
}

#[test]
fn backgrounds_with_reverse_complements_count_letters_of_either_strand() {
    let opts = BkgOptions {
        k: vec![1],
        rc: true,
        ..Default::default()
    };
    // AAC and its reverse complement GTT: the G and T should be counted too.
    let rows = get_bkg(&["AAC"], &[], &opts).unwrap();
    assert_eq!(
        bkg_counts(&rows),
        pairs(&[("A", 2), ("C", 1), ("G", 1), ("T", 2)])
    );
}

#[test]
fn windowed_backgrounds() {
    let seqs = ["AAAACCCC", "GGGGTTTT"];
    let opts = BkgOptions {
        k: vec![1],
        window: Some(WindowOptions {
            size: 4.0,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    let rows = get_bkg(&seqs, &[], &opts).unwrap();
    // Merged windows: ordered by k-mer, then window.
    let got: Vec<(&str, usize, u64)> = rows
        .iter()
        .map(|r| (r.kmer.as_str(), r.start.unwrap(), r.count))
        .collect();
    assert_eq!(
        got,
        [
            ("A", 0, 4),
            ("A", 4, 0),
            ("C", 0, 0),
            ("C", 4, 4),
            ("G", 0, 4),
            ("G", 4, 0),
            ("T", 0, 0),
            ("T", 4, 4)
        ]
    );
    assert!(rows.iter().all(|r| r.end == Some(r.start.unwrap() + 4)));
    assert!(approx(rows[0].probability, 0.5));

    // Per sequence, ordered by name, k-mer and window.
    let names = vec!["zeta".to_string(), "alpha".to_string()];
    let per = BkgOptions {
        merge: false,
        window: Some(WindowOptions {
            size: 0.5,
            overlap: 0.0,
        }),
        ..opts.clone()
    };
    let rows = get_bkg(&seqs, &names, &per).unwrap();
    assert_eq!(rows.len(), 2 * 4 * 2);
    assert!(rows[..8].iter().all(|r| r.sequence == Some(1)));
    assert_eq!((rows[0].kmer.as_str(), rows[0].count), ("A", 0));
    let g: Vec<u64> = rows
        .iter()
        .filter(|r| r.sequence == Some(1) && r.kmer == "G")
        .map(|r| r.count)
        .collect();
    assert_eq!(g, [4, 0]);

    // Fractional overlaps and multiple k.
    let frac = BkgOptions {
        k: vec![1, 2],
        merge: false,
        window: Some(WindowOptions {
            size: 4.0,
            overlap: 0.25,
        }),
        ..opts.clone()
    };
    let rows = get_bkg(&seqs[..1], &[], &frac).unwrap();
    let starts: std::collections::BTreeSet<usize> = rows.iter().map(|r| r.start.unwrap()).collect();
    assert_eq!(starts.into_iter().collect::<Vec<_>>(), [0, 2, 4]);
}

#[test]
fn background_errors() {
    let none: [&str; 0] = [];
    let opts = BkgOptions::default();
    assert_eq!(arg_msg(get_bkg(&none, &[], &opts)), "no sequences given");
    for k in [vec![], vec![0, 1]] {
        assert_eq!(
            arg_msg(get_bkg(
                &["ACGT"],
                &[],
                &BkgOptions {
                    k,
                    ..Default::default()
                }
            )),
            "k must be at least 1"
        );
    }
    assert_eq!(
        arg_msg(get_bkg(&["AC"], &[], &opts)),
        "k must not be greater than the length of the shortest sequence"
    );
    let win = |size, overlap, merge| BkgOptions {
        k: vec![2],
        merge,
        window: Some(WindowOptions { size, overlap }),
        ..Default::default()
    };
    let seqs = ["ACGTACGT", "ACGTAC"];
    assert_eq!(
        arg_msg(get_bkg(&seqs, &[], &win(0.1, 0.0, false))),
        "the window size is too small for these sequences (resolved to 0)"
    );
    assert_eq!(
        arg_msg(get_bkg(&seqs, &[], &win(3.0, 3.0, false))),
        "the window overlap must be smaller than the window size"
    );
    assert_eq!(
        arg_msg(get_bkg(&seqs, &[], &win(1.0, 0.0, false))),
        "the window size must be at least k"
    );
    assert_eq!(
        arg_msg(get_bkg(&seqs, &[], &win(7.0, 0.0, false))),
        "incorrect window and/or overlap size"
    );
    assert_eq!(
        arg_msg(get_bkg(&seqs, &[], &win(3.0, 0.0, true))),
        "merged windowed counts need all sequences to be of equal length"
    );
    assert_eq!(
        arg_msg(get_bkg(&["ACGTACG"], &[], &win(3.0, 0.0, false))),
        "k must not be greater than the smallest window"
    );
}

#[test]
fn backgrounds_of_non_ascii_windows() {
    let opts = BkgOptions {
        k: vec![1],
        merge: false,
        window: Some(WindowOptions {
            size: 2.0,
            overlap: 0.0,
        }),
        ..Default::default()
    };
    let rows = get_bkg(&["αβαα"], &[], &opts).unwrap();
    let got: Vec<(&str, usize, u64)> = rows
        .iter()
        .map(|r| (r.kmer.as_str(), r.start.unwrap(), r.count))
        .collect();
    assert_eq!(got, [("α", 0, 1), ("α", 2, 2), ("β", 0, 1), ("β", 2, 0)]);
}

// Encoded sequences.

#[test]
fn encoding_keeps_unknown_letters() {
    let seqs = ["acgtNNNNACGTRYacgtn", "", "TTTT-ACGT"];
    let e = EncodedSequences::new(&seqs, &Alphabet::Dna, 0).unwrap();
    assert_eq!(e.alphabet(), &Alphabet::Dna);
    assert_eq!((e.len(), e.is_empty(), e.total_len()), (3, false, 28));
    assert_eq!(e.decode(0), "ACGTNNNNACGTRYACGTn");
    assert_eq!(e.decode(1), "");
    assert_eq!(e.decode(2), "TTTT-ACGT");
    // Pieces of a sequence, including ones crossing unknown runs.
    let mut out = Vec::new();
    e.decode_into(0, 2, 6, &mut out);
    e.decode_into(0, 11, 14, &mut out);
    e.decode_into(0, 18, 19, &mut out);
    assert_eq!(out, b"GTNNTRYn");
    assert!(e.heap_bytes() > 0);
    // Every alphabet: RNA treats T as unknown, amino acids match either case,
    // and custom letters match exactly.
    let rna = EncodedSequences::new(&["acguTUx"], &Alphabet::Rna, 1).unwrap();
    assert_eq!(rna.decode(0), "ACGUTUx");
    let aa = EncodedSequences::new(&["mkwy*xBz"], &Alphabet::Aa, 1).unwrap();
    assert_eq!(aa.decode(0), "MKWY*xBz");
    assert!(aa.heap_bytes() >= 8);
    let custom = Alphabet::custom("01xY".chars()).unwrap();
    let c = EncodedSequences::new(&["0110xyYX2"], &custom, 1).unwrap();
    assert_eq!(c.decode(0), "0110xyYX2");
    let mut out = Vec::new();
    c.decode_into(0, 4, 8, &mut out);
    assert_eq!(out, b"xyYX");
}

#[test]
fn encoding_by_records_matches_whole_sequences() {
    for alphabet in [
        Alphabet::Dna,
        Alphabet::Aa,
        Alphabet::custom("ACGT".chars()).unwrap(),
    ] {
        let seqs = [
            random_dna(37, 1) + "NNxx" + &random_dna(9, 2),
            random_dna(3, 3),
        ];
        let whole = EncodedSequences::new(&seqs, &alphabet, 2).unwrap();
        let mut pieces = EncodedSequences::empty(&alphabet).unwrap();
        assert!(pieces.is_empty());
        for s in &seqs {
            let mut r = pieces.record();
            assert!(r.is_empty());
            r.reserve(4);
            for p in s.as_bytes().chunks(5) {
                r.push(p);
            }
            assert_eq!(r.len(), s.len());
            pieces.push(r);
        }
        pieces.shrink_to_fit();
        assert_eq!(pieces, whole, "{alphabet}");
        assert_eq!(pieces.seq_len(0), 50);
        assert_eq!(pieces.decode(0), whole.decode(0));
    }
}

#[test]
fn encoding_rejects_unusable_alphabets() {
    let wide = Alphabet::custom((0..300u32).map(|i| char::from_u32(0x100 + i).unwrap())).unwrap();
    assert_eq!(
        unsupported_msg(EncodedSequences::empty(&wide)),
        "alphabets of more than 255 letters cannot be scanned"
    );
    let accented = Alphabet::custom("aé".chars()).unwrap();
    assert_eq!(
        unsupported_msg(EncodedSequences::new(&["a"], &accented, 1)),
        "alphabets with letters outside ASCII cannot be scanned"
    );
    // All 95 printable ASCII letters are fine.
    let printable = Alphabet::custom((32u8..127).map(char::from)).unwrap();
    let e = EncodedSequences::new(&["Hello, world!"], &printable, 1).unwrap();
    assert_eq!(e.decode(0), "Hello, world!");
}

#[test]
#[should_panic(expected = "a record encoded for another alphabet")]
fn records_cannot_cross_alphabets() {
    let dna = EncodedSequences::empty(&Alphabet::Dna).unwrap();
    let mut aa = EncodedSequences::empty(&Alphabet::Aa).unwrap();
    aa.push(dna.record());
}

// Masking.

#[test]
fn masking_patterns() {
    let seqs = ["acgTTTacgt", "GGGGG", "AC"];
    assert_eq!(
        mask_seqs(&seqs, "ACG", false, 'N', &Alphabet::Dna).unwrap(),
        ["NNNTTTNNNt", "GGGGG", "AC"]
    );
    // Overlapping matches, and reverse complements.
    assert_eq!(
        mask_seqs(&seqs, "GGG", false, '-', &Alphabet::Dna).unwrap()[1],
        "-----"
    );
    assert_eq!(
        mask_seqs(&seqs, "AAA", true, 'x', &Alphabet::Dna).unwrap()[0],
        "acgxxxacgt"
    );
    assert_eq!(
        mask_seqs(&["CCAUUGG"], "AAU", true, 'N', &Alphabet::Rna).unwrap(),
        ["CCNNNGG"]
    );
    // Amino acids and custom alphabets mask the pattern only.
    assert_eq!(
        mask_seqs(&["MKVKV"], "KV", false, 'X', &Alphabet::Aa).unwrap(),
        ["MXXXX"]
    );
    let custom = Alphabet::custom("xyz".chars()).unwrap();
    assert_eq!(
        mask_seqs(&["xyzzy"], "zz", false, '.', &custom).unwrap(),
        ["xy..y"]
    );
}

#[test]
fn masking_errors() {
    assert_eq!(
        arg_msg(mask_seqs(&["ACGT"], "", false, 'N', &Alphabet::Dna)),
        "the pattern must not be empty"
    );
    assert_eq!(
        arg_msg(mask_seqs(&["ACGT"], "AC", false, 'ñ', &Alphabet::Dna)),
        "the masking letter must be a single ASCII character"
    );
    let custom = Alphabet::custom("xy".chars()).unwrap();
    assert_eq!(
        arg_msg(mask_seqs(&["xy"], "x", true, 'N', &custom)),
        "rc is only valid for DNA and RNA sequences"
    );
    assert_eq!(
        arg_msg(mask_ranges("ACGT", &[(-1, 2)], 'N')),
        "the range -1-2 is outside the sequence (length 4)"
    );
    assert_eq!(
        arg_msg(mask_ranges("ACGT", &[(0, 1)], 'é')),
        "the masking letter must be a single ASCII character"
    );
}

#[test]
fn masking_ranges() {
    // Empty and reversed ranges do nothing, and overlaps are fine.
    assert_eq!(
        mask_ranges("ACGTACGT", &[(2, 2), (5, 3), (0, 3), (2, 4), (8, 8)], 'n').unwrap(),
        "nnnnACGT"
    );
    assert_eq!(mask_ranges("ACGT", &[(0, 4)], '*').unwrap(), "****");
    assert_eq!(mask_ranges("", &[], 'N').unwrap(), "");
}

// Removing overlapping hits.

fn hit(
    sequence: usize,
    motif: usize,
    strand: Strand,
    start: i64,
    end: i64,
    priority: f64,
) -> Interval {
    Interval {
        sequence,
        motif,
        strand,
        start,
        end,
        priority,
    }
}

#[test]
fn dedup_keeps_the_best_of_each_cluster() {
    use Strand::*;
    let hits = [
        hit(0, 0, Plus, 0, 10, 1.0),
        hit(0, 0, Plus, 5, 15, 3.0),
        hit(0, 0, Plus, 12, 20, 2.0),
        // Touching is not overlapping.
        hit(0, 0, Plus, 20, 25, 0.5),
        // Other sequences, motifs and strands form their own groups.
        hit(1, 0, Plus, 0, 10, 0.1),
        hit(0, 1, Plus, 0, 10, 0.1),
        hit(0, 0, Minus, 0, 10, 0.1),
    ];
    assert_eq!(
        dedup(&hits, false, false),
        [false, true, false, true, true, true, true]
    );
    // Ignoring strands and motifs merges those groups.
    assert_eq!(
        dedup(&hits, true, true),
        [false, true, false, true, true, false, false]
    );
    assert_eq!(
        dedup(&hits, true, false),
        [false, true, false, true, true, true, false]
    );
    assert_eq!(
        dedup(&hits, false, true),
        [false, true, false, true, true, false, true]
    );
}

#[test]
fn dedup_chains_and_ties() {
    use Strand::*;
    // A chain where the ends do not overlap each other: if the middle is
    // worse, both ends are kept.
    let chain = [
        hit(0, 0, Both, 0, 10, 5.0),
        hit(0, 0, Both, 8, 18, 1.0),
        hit(0, 0, Both, 16, 26, 5.0),
    ];
    assert_eq!(dedup(&chain, false, false), [true, false, true]);
    // Ties go to the earlier hit in the input, whatever its position.
    let tie = [hit(0, 0, Plus, 5, 15, 2.0), hit(0, 0, Plus, 0, 10, 2.0)];
    assert_eq!(dedup(&tie, false, false), [true, false]);
    assert!(dedup(&[], false, false).is_empty());
}

// FASTA.

#[test]
fn fasta_text_round_trips() {
    let records = vec![
        ("chr1 first record".to_string(), random_dna(23, 1)),
        ("empty".to_string(), String::new()),
        ("αβ".to_string(), "αβγδ".to_string()),
    ];
    // Widths that do not split the non-ASCII record; see
    // fasta_wraps_non_ascii_sequences_by_letter for the others.
    for width in [0, 2, 80] {
        let text = fasta::write(&records, width);
        assert_eq!(fasta::read(&text).unwrap(), records, "width {width}");
    }
    assert_eq!(
        fasta::write(&[("a", "ACGTACG"), ("b", "")], 3),
        ">a\nACG\nTAC\nG\n>b\n\n"
    );
    assert_eq!(fasta::write(&[("a", "ACGTACG")], 0), ">a\nACGTACG\n");
    let messy = ";comment\n\n>  one \r\nAC GT\r\n;inner comment\n\tac\n>two\n";
    assert_eq!(
        fasta::read(messy).unwrap(),
        [
            ("one".to_string(), "ACGTac".to_string()),
            ("two".to_string(), String::new())
        ]
    );
    assert!(fasta::read("").unwrap().is_empty());
    assert_eq!(
        arg_msg(fasta::read("ACGT\n>a\n")),
        "FASTA text must start with a '>' header line"
    );
}

#[test]
fn fasta_files_plain_and_gzipped() {
    let text = ">s1 desc\nACGTN\nacgt\n>s2\nMKV\n";
    let plain = temp_path("plain.fa");
    let gz = temp_path("packed.fa.gz");
    std::fs::write(&plain, text).unwrap();
    std::fs::write(&gz, gzip(text.as_bytes())).unwrap();
    let want = fasta::read(text).unwrap();
    assert_eq!(fasta::read_path(&plain, 0).unwrap(), want);
    assert_eq!(fasta::read_path(&gz, 0).unwrap(), want);
    // Skipping the first three lines starts at the second record.
    assert_eq!(fasta::read_path(&gz, 3).unwrap(), want[1..]);
    assert_eq!(
        arg_msg(fasta::read_path(&gz, 1)),
        "FASTA text must start with a '>' header line"
    );

    let (names, seqs) = fasta::read_encoded(&gz, &Alphabet::Aa, 0).unwrap();
    assert_eq!(names, ["s1 desc", "s2"]);
    assert_eq!(seqs.alphabet(), &Alphabet::Aa);
    assert_eq!(seqs.decode(0), "ACGTNACGT");
    assert_eq!(seqs.decode(1), "MKV");
    let (names, seqs) = fasta::read_encoded(&plain, &Alphabet::Dna, 3).unwrap();
    assert_eq!((names.len(), seqs.decode(0)), (1, "MKV".to_string()));
    // A file without records gives nothing.
    std::fs::write(&plain, "; only a comment\n").unwrap();
    let (names, seqs) = fasta::read_encoded(&plain, &Alphabet::Dna, 0).unwrap();
    assert!(names.is_empty() && seqs.is_empty());
    for p in [&plain, &gz] {
        std::fs::remove_file(p).unwrap();
    }

    let missing = temp_path("missing.fa");
    let msg = arg_msg(fasta::read_path(&missing, 0));
    assert!(msg.starts_with("could not open "), "{msg}");
    assert!(fasta::read_encoded(&missing, &Alphabet::Dna, 0).is_err());
    let wide = Alphabet::custom("aé".chars()).unwrap();
    assert!(matches!(
        fasta::read_encoded(&missing, &wide, 0),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn fasta_rejects_invalid_text() {
    let bad: &[u8] = b">a\nAC\xffGT\n";
    let msg = arg_msg(fasta::read_reader(bad, 0));
    assert!(msg.starts_with("could not read the FASTA file"), "{msg}");
}

#[test]
fn fasta_batches_pass_on_errors() {
    let text = ">a\nACGT\n>b\nACGT\n>c\nACGT\n";
    let mut calls = 0;
    let r = fasta::read_encoded_batches(text.as_bytes(), &Alphabet::Dna, 0, 1, |names, _| {
        calls += 1;
        if names == ["b"] {
            Err(Error::Argument("stop at b".into()))
        } else {
            Ok(())
        }
    });
    assert_eq!(arg_msg(r), "stop at b");
    assert_eq!(calls, 2, "reading stops at the first error");
    let r = fasta::read_encoded_batches(text.as_bytes(), &Alphabet::Dna, 0, 100, |_, _| {
        Err(Error::Argument("last batch".into()))
    });
    assert_eq!(arg_msg(r), "last batch");
    let r = fasta::read_encoded_batches("ACGT\n".as_bytes(), &Alphabet::Dna, 0, 1, |_, _| Ok(()));
    assert_eq!(arg_msg(r), "FASTA text must start with a '>' header line");
}

#[test]
fn fasta_wraps_non_ascii_sequences_by_letter() {
    assert_eq!(fasta::write(&[("g", "αβγ")], 2), ">g\nαβ\nγ\n");
}

// MEME.

fn dna_example(name: &str) -> Motif {
    Motif::builder(
        Matrix::from_columns(&[
            [0.5, 0.25, 0.125, 0.125],
            [0.0, 0.0, 1.0, 0.0],
            [0.1, 0.2, 0.3, 0.4],
        ])
        .unwrap(),
        Alphabet::Dna,
    )
    .name(name)
    .motif_type(Some(MotifType::Ppm))
    .build()
    .unwrap()
}

#[test]
fn meme_round_trips_every_alphabet() {
    let bkg = vec![0.3, 0.2, 0.2, 0.3];
    let dna = Motif::builder(dna_example("d").matrix.clone(), Alphabet::Dna)
        .name("d1")
        .altname(Some("alt".into()))
        .motif_type(Some(MotifType::Ppm))
        .nsites(Some(17.0))
        .eval(Some(2.5e-8))
        .bkg(Some(bkg.clone()))
        .strand(Strand::Plus)
        .build()
        .unwrap();
    let text = meme::write(std::slice::from_ref(&dna), &meme::WriteOptions::default()).unwrap();
    assert!(text.starts_with("MEME version 5\n"));
    assert!(text.contains("strands: +\n"));
    assert!(text.contains("A 0.300000 C 0.200000 G 0.200000 T 0.300000"));
    assert!(text.contains("letter-probability matrix: alength= 4 w= 3 nsites= 17 E= 2.5e-08"));
    let back = meme::read(&text, &meme::ReadOptions::default()).unwrap();
    assert_eq!(back.len(), 1);
    let b = &back[0];
    assert_eq!((b.name.as_str(), b.altname.as_deref()), ("d1", Some("alt")));
    assert_eq!(
        (b.nsites, b.eval, b.strand),
        (Some(17.0), Some(2.5e-8), Strand::Plus)
    );
    assert_eq!(b.bkg, bkg);
    assert_eq!(
        b.extrainfo,
        [("eval.string".to_string(), "2.5e-08".to_string())]
    );
    assert_matrix_close(&b.matrix, &dna.matrix, 1e-6, "DNA");

    let rna = ppm(Alphabet::Rna, &[&[0.1, 0.2, 0.3, 0.4]], "r");
    let aa_col: Vec<f64> = (1..=20).map(|i| i as f64 / 210.0).collect();
    let aa = ppm(Alphabet::Aa, &[&aa_col, &aa_col], "a");
    let custom = ppm(
        Alphabet::custom("1Ab".chars()).unwrap(),
        &[&[0.5, 0.3, 0.2], &[0.0, 0.1, 0.9]],
        "c",
    );
    for m in [rna, aa, custom] {
        let text = meme::write(std::slice::from_ref(&m), &meme::WriteOptions::default()).unwrap();
        let back = meme::read(&text, &meme::ReadOptions::default()).unwrap();
        assert_eq!(back[0].alphabet, m.alphabet, "{}", m.alphabet);
        assert_eq!(back[0].strand, Strand::Both);
        assert_matrix_close(&back[0].matrix, &m.matrix, 1e-6, &m.alphabet.name());
    }
}

#[test]
fn meme_write_options() {
    let a = dna_example("a b=c ");
    let b = Motif::builder(a.matrix.clone(), Alphabet::Dna)
        .name("other")
        .bkg(Some(vec![0.4, 0.1, 0.1, 0.4]))
        .strand(Strand::Minus)
        .build()
        .unwrap();
    let text = meme::write(&[a.clone(), b.clone()], &meme::WriteOptions::default()).unwrap();
    // Names lose spaces and equals signs, differing backgrounds become
    // uniform, and mixed strands become both, as does a missing nsites.
    assert!(text.contains("MOTIF a-b.c\n"), "{text}");
    assert!(text.contains("A 0.250000 C 0.250000 G 0.250000 T 0.250000"));
    assert!(text.contains("strands: + -\n"));
    assert!(text.contains("nsites= 100\n"));
    let opts = meme::WriteOptions {
        version: 4,
        bkg: Some(vec![0.1, 0.4, 0.4, 0.1]),
        strand: Some("+".into()),
        ..Default::default()
    };
    let text = meme::write(std::slice::from_ref(&a), &opts).unwrap();
    assert!(text.starts_with("MEME version 4\n"));
    let back = meme::read(&text, &meme::ReadOptions::default()).unwrap();
    assert_eq!(back[0].bkg, [0.1, 0.4, 0.4, 0.1]);
    assert_eq!(back[0].strand, Strand::Plus);

    // Without a header the motifs can be appended to another file.
    let body = meme::write(
        std::slice::from_ref(&b),
        &meme::WriteOptions {
            header: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(body.starts_with("MOTIF other\n"));
    let both = meme::read(&(text + &body), &meme::ReadOptions::default()).unwrap();
    assert_eq!(both.len(), 2);
    assert_eq!(both[1].name, "other");
    assert_eq!(both[1].bkg, [0.1, 0.4, 0.4, 0.1], "the file's background");

    assert_eq!(
        arg_msg(meme::write(
            std::slice::from_ref(&a),
            &meme::WriteOptions {
                bkg: Some(vec![0.5, 0.5]),
                ..Default::default()
            }
        )),
        "bkg must have 4 values"
    );
}

#[test]
fn meme_contribution_matrices() {
    let cwm = motif(
        Alphabet::Dna,
        &[&[0.5, -0.25, 0.0, 0.1], &[-0.1, 0.0, 0.75, 0.0]],
        MotifType::Cwm,
        "cwm",
    );
    let opts = meme::WriteOptions {
        cwm: true,
        ..Default::default()
    };
    let text = meme::write(std::slice::from_ref(&cwm), &opts).unwrap();
    assert!(
        text.contains(" 0.500000  -0.250000  0.000000  0.100000"),
        "{text}"
    );
    let back = meme::read(&text, &meme::ReadOptions { cwm: true }).unwrap();
    assert_eq!(back[0].motif_type, MotifType::Cwm);
    assert_eq!(back[0].matrix, cwm.matrix);
    assert_eq!(
        arg_msg(meme::write(&[dna_example("p")], &opts)),
        "cwm = true requires every motif to be of type CWM"
    );
}

#[test]
fn meme_write_errors() {
    assert_eq!(
        arg_msg(meme::write(&[], &meme::WriteOptions::default())),
        "no motifs to write"
    );
    let rna = ppm(Alphabet::Rna, &[&[0.25; 4]], "r");
    assert_eq!(
        arg_msg(meme::write(
            &[dna_example("d"), rna],
            &meme::WriteOptions::default()
        )),
        "all motifs must use the same alphabet"
    );
    let odd = ppm(Alphabet::custom("@#".chars()).unwrap(), &[&[0.5, 0.5]], "o");
    assert!(
        arg_msg(meme::write(&[odd], &meme::WriteOptions::default()))
            .starts_with("MEME format motifs only allow the following alphabet symbols:")
    );
}

const MEME_HEAD: &str = "MEME version 5\n\nALPHABET= ACGT\n\n";

fn meme_err(text: &str) -> String {
    arg_msg(meme::read(text, &meme::ReadOptions::default()))
}

#[test]
fn meme_reader_errors() {
    assert_eq!(meme_err(MEME_HEAD), "no motifs found");
    assert_eq!(
        meme_err(&format!("{MEME_HEAD}MOTIF m\n")),
        "no motifs found"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}letter-probability matrix: w= 1\n0.25 0.25 0.25 0.25\nMOTIF m\n"
        )),
        "found a letter-probability matrix without a MOTIF line"
    );
    for w in ["", "w= 0", "w= 1.5", "w= x"] {
        assert_eq!(
            meme_err(&format!(
                "{MEME_HEAD}MOTIF m\nletter-probability matrix: {w}\n0.25 0.25 0.25 0.25\n"
            )),
            "motif 'm' has no valid width (w=)",
            "{w}"
        );
    }
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}MOTIF m\nletter-probability matrix: w= 2\n0.25 0.25 0.25 0.25\n"
        )),
        "motif 'm' is truncated"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}MOTIF m\nletter-probability matrix: w=1\n0.5 0.5 0\n"
        )),
        "motif 'm': alphabet length does not match motif length"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}MOTIF m\nletter-probability matrix: w=1\n0.5 0.5 0 x\n"
        )),
        "could not read a probability 'x' as a number"
    );
    assert_eq!(
        meme_err("ALPHABET \"xy\"\nx\ny\nMOTIF m\n"),
        "could not find the end of the custom alphabet"
    );
    assert_eq!(
        meme_err("ALPHABET= ACGTT\n"),
        "could not parse the MEME alphabet"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}Background letter frequencies\nA 0.25 C 0.25\n"
        )),
        "could not parse background frequencies, check that they match the alphabet"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}Background letter frequencies\nA 0.25 C\n"
        )),
        "could not parse background frequencies"
    );
    assert_eq!(
        meme_err(&format!(
            "{MEME_HEAD}Background letter frequencies\nA 0.25 C x G 1 T 1\n"
        )),
        "could not read a background frequency 'x' as a number"
    );
    assert_eq!(
        meme_err(&format!("{MEME_HEAD}strands: x\n")),
        "incorrect strand 'x' (must be one of +, -, +-)"
    );
}

#[test]
fn meme_reader_header_variants() {
    let body =
        "MOTIF m alt\nletter-probability matrix: alength= 4 w= 1 nsites=5\n0.1 0.2 0.3 0.4\n";
    let read =
        |head: &str| meme::read(&format!("{head}{body}"), &meme::ReadOptions::default()).unwrap();
    // No alphabet means DNA, no strands line means +, no background means
    // uniform.
    let m = &read("MEME version 4\n")[0];
    assert_eq!(
        (m.alphabet.clone(), m.strand),
        (Alphabet::Dna, Strand::Plus)
    );
    assert_eq!(m.bkg, [0.25; 4]);
    assert_eq!(
        (m.altname.as_deref(), m.nsites, m.eval),
        (Some("alt"), Some(5.0), None)
    );
    assert!(m.extrainfo.is_empty());
    assert_eq!(read("strands: -\n")[0].strand, Strand::Minus);
    assert_eq!(read("strands: + -\n")[0].strand, Strand::Both);
    // Backgrounds over several lines, with letters the alphabet does not
    // know, are taken in order.
    let m = &read("Background letter frequencies\nW 0.1 X 0.2\nY 0.3 Z 0.4\n")[0];
    assert_eq!(m.bkg, [0.1, 0.2, 0.3, 0.4]);
    let m = &read("Background letter frequencies\nT 0.1 G 0.2 C 0.3 A 0.4\n")[0];
    assert_eq!(m.bkg, [0.4, 0.3, 0.2, 0.1]);
    // RNA and protein alphabets, and a custom one whose letters are listed
    // in MEME's order.
    let rna = meme::read(
        "ALPHABET= ACGU\nMOTIF r\nletter-probability matrix: w= 1\n1 0 0 0\n",
        &meme::ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(rna[0].alphabet, Alphabet::Rna);
    let custom = meme::read(
        "ALPHABET \"mine\"\nb ~ 1\nZ\nEND ALPHABET\nMOTIF c\nletter-probability matrix: w= 1\n0.2 0.3 0.5\n",
        &meme::ReadOptions::default(),
    )
    .unwrap();
    assert_eq!(custom[0].alphabet, Alphabet::custom("1Zb".chars()).unwrap());
    // The file lists Z, b, 1; the alphabet sorts them as 1, Z, b.
    assert_matrix_close(
        &custom[0].matrix,
        &Matrix::from_columns(&[[0.5, 0.2, 0.3]]).unwrap(),
        1e-12,
        "custom",
    );
}

const MEME_OUTPUT: &str = "MEME version 5

ALPHABET= ACGT

MOTIF AAC MEME-1 width = 3 sites = 2
letter-probability matrix: alength= 4 w= 3 nsites= 2 E= 1e-3
 1 0 0 0
 1 0 0 0
 0 1 0 0

--------------------------------------------------------------------------------
\tMotif AAC MEME-1 sites sorted by position p-value
--------------------------------------------------------------------------------
Sequence name            Strand  Start   P-value                 Site
-------------            ------  ----- ---------            ---------
s1                           +      5  1.00e-04 CCCC AAC GGGG
s2                           -      1  2.00e-04 AAC TTT
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
\tMotif AAC MEME-1 block diagrams
--------------------------------------------------------------------------------
\tMotif AAC MEME-1 in BLOCKS format
--------------------------------------------------------------------------------
BL   MOTIF AAC width=3 seqs=2
s1                       (    5) AAC  1
s2                       (    1) AAC  1
//

********************************************************************************
Combined block diagrams: non-overlapping sites with p-value < 0.0001
--------------------------------------------------------------------------------
SEQUENCE NAME            COMBINED P-VALUE  MOTIF DIAGRAM
-------------            ----------------  -------------
s1                               1.2e-05  4_[+1(1.00e-04)]_\\
                                          4
s2                               3.4e-04  [-1(2.00e-04)]_3
--------------------------------------------------------------------------------
";

#[test]
fn meme_output_sections() {
    let sites = meme::read_sites(MEME_OUTPUT).unwrap();
    assert_eq!(sites.len(), 1);
    assert_eq!(sites[0].motif, "AAC");
    let got: Vec<(&str, &str)> = sites[0]
        .sites
        .iter()
        .map(|s| (s.sequence.as_str(), s.site.as_str()))
        .collect();
    assert_eq!(got, [("s1", "AAC"), ("s2", "AAC")]);

    let pv = meme::read_site_pvalues(MEME_OUTPUT).unwrap();
    assert_eq!(pv[0].motif, "AAC");
    let s = &pv[0].sites;
    assert_eq!(s.len(), 2);
    assert_eq!(
        (s[0].sequence.as_str(), s[0].strand, s[0].position),
        ("s1", Some('+'), 4.0)
    );
    assert_eq!((s[0].pvalue, s[0].site.as_str()), (1e-4, "AAC"));
    // Without a left flank the site is found by the motif's width.
    assert_eq!(
        (s[1].strand, s[1].position, s[1].site.as_str()),
        (Some('-'), 0.0, "AAC")
    );

    let combined = meme::read_combined(MEME_OUTPUT).unwrap();
    assert_eq!(combined.len(), 2);
    assert_eq!(combined[0].sequence, "s1");
    assert_eq!(combined[0].pvalue, 1.2e-5);
    assert_eq!(combined[0].diagram, "4_[+1(1.00e-04)]_4");
    assert_eq!(combined[1].diagram, "[-1(2.00e-04)]_3");
    assert!(meme::read_combined(MEME_HEAD).unwrap().is_empty());

    let motifs = meme::read(MEME_OUTPUT, &meme::ReadOptions::default()).unwrap();
    assert_eq!(
        (motifs[0].name.as_str(), motifs[0].altname.as_deref()),
        ("AAC", Some("MEME-1"))
    );
    assert_eq!(motifs[0].eval, Some(1e-3));
}

#[test]
fn meme_site_tables_without_strands() {
    let text = "MOTIF m\n\tMotif m sites sorted by position p-value\n\
                Sequence name  Start  P-value  Site\n\
                s1  3  0.01  AC GTA CG\n";
    let pv = meme::read_site_pvalues(text).unwrap();
    assert_eq!(pv[0].sites[0].strand, None);
    assert_eq!(
        (pv[0].sites[0].position, pv[0].sites[0].site.as_str()),
        (2.0, "GTA")
    );
    assert_eq!(
        arg_msg(meme::read_site_pvalues(&format!("{text}s2 4\n"))),
        "could not read site line 's2 4'"
    );
    assert_eq!(
        arg_msg(meme::read_site_pvalues(&format!("{text}s2 x 0.1 A C G\n"))),
        "could not read a site position 'x' as a number"
    );
}

#[test]
fn meme_blocks_errors() {
    assert_eq!(
        arg_msg(meme::read_sites("Motif m in BLOCKS format\n")),
        "BLOCKS section is truncated"
    );
    assert_eq!(
        arg_msg(meme::read_sites(
            "Motif m in BLOCKS format\nBL MOTIF m width=3\n"
        )),
        "BLOCKS section for 'm' has no seqs="
    );
    assert_eq!(
        arg_msg(meme::read_sites(
            "Motif m in BLOCKS format\nBL MOTIF m width=3 seqs=2\ns1 ( 1) AAA 1\n"
        )),
        "BLOCKS section for 'm' is truncated"
    );
    assert_eq!(
        arg_msg(meme::read_sites(
            "Motif m in BLOCKS format\nBL MOTIF m width=3 seqs=1\ns1 1 AAA\n"
        )),
        "could not read BLOCKS line 's1 1 AAA'"
    );
    assert!(meme::read_sites(MEME_HEAD).unwrap().is_empty());
}

#[test]
fn meme_alphabet_definitions() {
    let def = meme::AlphabetDefinition {
        core: "ACGT".into(),
        complements: Some("TGCA".into()),
        ambiguity: vec![("N".into(), "ACGT".into()), ("R".into(), "AG".into())],
        like: Some("DNA".into()),
        name: Some("my dna".into()),
        ..Default::default()
    };
    assert_eq!(
        meme::write_alphabet(&def).unwrap(),
        "ALPHABET \"my dna\" DNA-LIKE\n\n# Core symbols\nA ~ T\nC ~ G\n\n\
         # Ambiguous symbols\nN = ACGT\nR = AG\n"
    );
    let named = meme::AlphabetDefinition {
        core: "XY".into(),
        letter_names: vec![('X', "ex".into()), ('Y', "why".into())],
        colours: vec![('X', "FF0000".into()), ('Y', "0000FF".into())],
        ..Default::default()
    };
    assert_eq!(
        meme::write_alphabet(&named).unwrap(),
        "ALPHABET\n\n# Core symbols\nX \"ex\" FF0000\nY \"why\" 0000FF\n"
    );
    let errors = [
        (
            meme::AlphabetDefinition {
                like: Some("XNA".into()),
                ..def.clone()
            },
            "like must be one of DNA, RNA and PROTEIN",
        ),
        (
            meme::AlphabetDefinition {
                complements: Some("TGC".into()),
                ..def.clone()
            },
            "core and complements do not have the same number of letters",
        ),
        (
            meme::AlphabetDefinition {
                letter_names: vec![('X', "ex".into())],
                ..named.clone()
            },
            "the letter 'Y' has no name in letter_names",
        ),
        (
            meme::AlphabetDefinition {
                colours: vec![('Y', "00FF00".into())],
                ..named.clone()
            },
            "the letter 'X' has no colour",
        ),
    ];
    for (d, msg) in errors {
        assert_eq!(arg_msg(meme::write_alphabet(&d)), msg);
    }
}

// HOMER.

#[test]
fn homer_round_trips() {
    let m = dna_example("hm");
    let text = homer::write(
        std::slice::from_ref(&m),
        &homer::Threshold::LogoddsAbs(vec![3.5]),
    )
    .unwrap();
    let consensus = m.consensus.clone().unwrap();
    assert!(
        text.starts_with(&format!(">{consensus}\thm\t3.5\n")),
        "{text}"
    );
    assert!(text.contains("0.500\t0.250\t0.125\t0.125\n"));
    let back = homer::read(&text).unwrap();
    assert_eq!(back[0].name, "hm");
    assert_eq!(back[0].alphabet, Alphabet::Dna);
    assert_eq!(
        back[0].extrainfo,
        [("logodds.threshold".to_string(), "3.5".to_string())]
    );
    assert_matrix_close(&back[0].matrix, &m.matrix, 5e-4, "HOMER");
    assert_eq!(
        (back[0].nsites, back[0].bkgsites, back[0].pval),
        (None, None, None)
    );
    // RNA motifs are written as their probabilities and read back as DNA.
    let rna = ppm(Alphabet::Rna, &[&[0.1, 0.2, 0.3, 0.4]], "r");
    let back = homer::read(&homer::write(&[rna], &homer::Threshold::default()).unwrap()).unwrap();
    assert_eq!(back[0].alphabet, Alphabet::Dna);
}

#[test]
fn homer_thresholds() {
    let a = one_hot("ACG");
    let b = one_hot("TTTT");
    let header_value = |text: &str, i: usize| -> f64 {
        let header = text.lines().filter(|l| l.starts_with('>')).nth(i).unwrap();
        header.split('\t').nth(2).unwrap().parse().unwrap()
    };
    // A fraction of the best score, ln(1 / 0.25) per position.
    let text = homer::write(
        &[a.clone(), b.clone()],
        &homer::Threshold::Logodds(vec![0.5]),
    )
    .unwrap();
    assert!((header_value(&text, 0) - 0.5 * 3.0 * 4f64.ln()).abs() < 1e-12);
    assert!((header_value(&text, 1) - 0.5 * 4.0 * 4f64.ln()).abs() < 1e-12);
    // One value per motif.
    let text = homer::write(
        &[a.clone(), b.clone()],
        &homer::Threshold::LogoddsAbs(vec![1.0, 2.0]),
    )
    .unwrap();
    assert_eq!((header_value(&text, 0), header_value(&text, 1)), (1.0, 2.0));
    // Only the best of the 64 sites reaches a P-value of 1/64, so the cutoff
    // sits between the second best score and a little above the best.
    let text = homer::write(
        std::slice::from_ref(&a),
        &homer::Threshold::Pvalue(vec![1.0 / 64.0]),
    )
    .unwrap();
    let best = 3.0 * 4f64.ln();
    let second = 2.0 * 4f64.ln() + (0.001f64 / 0.25).ln();
    let t = header_value(&text, 0);
    assert!(t > second && t <= best + 0.003, "{t}");
    // Non-finite thresholds are written as R writes them.
    let text = homer::write(
        std::slice::from_ref(&a),
        &homer::Threshold::LogoddsAbs(vec![f64::INFINITY]),
    )
    .unwrap();
    assert!(text.starts_with(">ACG\tACG\tInf\n"), "{text}");
    let text = homer::write(
        std::slice::from_ref(&a),
        &homer::Threshold::LogoddsAbs(vec![f64::NAN]),
    )
    .unwrap();
    assert!(text.starts_with(">ACG\tACG\tNaN\n"), "{text}");

    assert_eq!(
        arg_msg(homer::write(
            &[a.clone(), b],
            &homer::Threshold::Logodds(vec![0.1, 0.2, 0.3])
        )),
        "a single threshold or one per motif should be provided"
    );
    let aa = ppm(Alphabet::Aa, &[&[0.05; 20]], "aa");
    assert_eq!(
        unsupported_msg(homer::write(&[aa], &homer::Threshold::default())),
        "HOMER format only supports DNA and RNA motifs"
    );
}

#[test]
fn homer_reads_full_headers() {
    let text = ">ACGT\tRUNX(Runt)/HPC7-Runx1-ChIP-Seq/Homer\t6.5\t-1000.0\t0\t\
                T:500.0(25.00%),B:200.0(5.00%),P:1e-50\n\
                0.97 0.01 0.01 0.01\n0.01 0.97 0.01 0.01\n\
                >GG  plain/name  4\n0.1 0.2 0.3 0.4\n";
    let m = homer::read(text).unwrap();
    assert_eq!(m[0].name, "RUNX");
    assert_eq!(m[0].family.as_deref(), Some("Runt"));
    assert_eq!(
        (m[0].nsites, m[0].bkgsites, m[0].pval),
        (Some(500.0), Some(200.0), Some(1e-50))
    );
    assert_eq!(m[0].extrainfo[0].1, "6.5");
    assert_eq!(m[0].width(), 2);
    // Without tabs, fields are split at spaces; a slash without a family
    // keeps the first part.
    assert_eq!(m[1].name, "plain");
    assert_eq!(m[1].family, None);
}

#[test]
fn homer_reader_errors() {
    assert_eq!(
        arg_msg(homer::read("0.25 0.25 0.25 0.25\n")),
        "no motifs found"
    );
    assert_eq!(
        arg_msg(homer::read(">A\tm\t1\n0.5 0.5 0\n")),
        "motif 'm': expected 4 values per row, found 3"
    );
    assert_eq!(
        arg_msg(homer::read(">A\tm\t1\n>C\tn\t1\n0.25 0.25 0.25 0.25\n")),
        "empty motif matrix for motif 'm'"
    );
    assert_eq!(
        arg_msg(homer::read(">A\tm\t1\n0.5 0.5 0 x\n")),
        "could not read a probability 'x' as a number"
    );
    // NA is not a number in motif files.
    assert_eq!(
        arg_msg(homer::read(">A\tm\t1\n0.5 0.5 0 NA\n")),
        "could not read a probability 'NA' as a number"
    );
}

// JASPAR.

#[test]
fn jaspar_round_trips() {
    let pcm = Motif::builder(
        Matrix::from_columns(&[[5.0, 5.0, 5.0, 5.0], [20.0, 0.0, 0.0, 0.0]]).unwrap(),
        Alphabet::Dna,
    )
    .name("MA0001.1")
    .altname(Some("AGL3".into()))
    .motif_type(Some(MotifType::Pcm))
    .build()
    .unwrap();
    let text = jaspar::write(std::slice::from_ref(&pcm)).unwrap();
    assert_eq!(
        text,
        ">MA0001.1 AGL3\nA [  5 20 ]\nC [  5  0 ]\nG [  5  0 ]\nT [  5  0 ]\n\n"
    );
    let back = jaspar::read(&text).unwrap();
    assert_eq!(back[0].matrix, pcm.matrix);
    assert_eq!(
        (back[0].name.as_str(), back[0].altname.as_deref()),
        ("MA0001.1", Some("AGL3"))
    );
    assert_eq!(
        (back[0].motif_type, back[0].nsites),
        (MotifType::Pcm, Some(20.0))
    );

    // Column widths grow with the site count.
    for (nsites, width) in [(500.0, 3), (5000.0, 4), (50000.0, 5), (500000.0, 6)] {
        let big = Motif::builder(
            Matrix::from_columns(&[[nsites, 0.0, 0.0, 0.0]]).unwrap(),
            Alphabet::Dna,
        )
        .motif_type(Some(MotifType::Pcm))
        .build()
        .unwrap();
        let text = jaspar::write(&[big]).unwrap();
        let c_row = text.lines().nth(2).unwrap();
        assert_eq!(c_row, format!("C [ {:>width$} ]", 0), "{nsites}");
    }

    // RNA, protein and probability matrices.
    let rna = motif(Alphabet::Rna, &[&[1.0, 2.0, 3.0, 4.0]], MotifType::Pcm, "r");
    let back = jaspar::read(&jaspar::write(std::slice::from_ref(&rna)).unwrap()).unwrap();
    assert_eq!(
        (back[0].alphabet.clone(), back[0].matrix.clone()),
        (Alphabet::Rna, rna.matrix.clone())
    );
    let col: Vec<f64> = (0..20).map(|i| (i % 3) as f64).collect();
    let aa = motif(Alphabet::Aa, &[&col], MotifType::Pcm, "aa");
    let back = jaspar::read(&jaspar::write(std::slice::from_ref(&aa)).unwrap()).unwrap();
    assert_eq!(
        (back[0].alphabet.clone(), back[0].matrix.clone()),
        (Alphabet::Aa, aa.matrix.clone())
    );
    let p = Motif::builder(dna_example("p").matrix.clone(), Alphabet::Dna)
        .name("p")
        .nsites(Some(40.0))
        .pseudocount(0.0)
        .build()
        .unwrap();
    let back = jaspar::read(&jaspar::write(&[p]).unwrap()).unwrap();
    assert_eq!(back[0].matrix.col(0), [20.0, 10.0, 5.0, 5.0]);
    assert_eq!(back[0].altname, None);
}

#[test]
fn jaspar_errors() {
    // A bracket-only row in a labelled matrix is skipped.
    let m = jaspar::read(">m\nA [ 1 ]\nC [ 1 ]\n[ ]\nG [ 1 ]\nT [ 1 ]\n").unwrap();
    assert_eq!(m[0].matrix.col(0), [1.0; 4]);
    let custom = ppm(Alphabet::custom("xy".chars()).unwrap(), &[&[0.5, 0.5]], "c");
    assert_eq!(
        unsupported_msg(jaspar::write(&[custom])),
        "JASPAR format only supports DNA, RNA and AA motifs"
    );
    assert_eq!(arg_msg(jaspar::read("A [ 1 2 ]\n")), "no motifs found");
    assert_eq!(
        arg_msg(jaspar::read(">m\n1 2\n3 4\n5 6\n7 8\n9 10\n")),
        "motif 'm' has 5 unlabelled rows; label the rows to give its alphabet"
    );
    assert_eq!(
        arg_msg(jaspar::read(
            ">m\nA [ 1 2 ]\nC [ 3 ]\nG [ 1 1 ]\nT [ 1 1 ]\n"
        )),
        "rows of motif 'm' have different lengths"
    );
    assert_eq!(
        arg_msg(jaspar::read(">m\nA [ 1 ]\nCG [ 3 ]\n")),
        "unexpected row label 'CG' in motif 'm'"
    );
    assert_eq!(
        arg_msg(jaspar::read(">m\n>n\n1\n1\n1\n1\n")),
        "empty motif matrix for motif 'm'"
    );
    assert_eq!(
        arg_msg(jaspar::read(">m\nA [ ]\nC [ ]\n")),
        "empty motif matrix for motif 'm'"
    );
    assert_eq!(
        arg_msg(jaspar::read(">m\nA [ 1 x ]\n")),
        "could not read a count 'x' as a number"
    );
}

// CIS-BP.

#[test]
fn cisbp_round_trips_exactly() {
    let m = Motif::builder(
        Matrix::from_columns(&[
            [0.1, 0.2, 0.3, 0.4],
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        ])
        .unwrap(),
        Alphabet::Dna,
    )
    .name("M00001_2.00")
    .altname(Some("Gene X".into()))
    .family(Some("Homeodomain".into()))
    .organism(Some("Mus musculus".into()))
    .motif_type(Some(MotifType::Ppm))
    .build()
    .unwrap();
    let rna = ppm(Alphabet::Rna, &[&[0.7, 0.1, 0.1, 0.1]], "r");
    let text = cisbp::write(&[m.clone(), rna.clone()]).unwrap();
    assert!(text.starts_with(
        "TF\tM00001_2.00\nTF Name\tGene X\nFamily\tHomeodomain\nSpecies\tMus musculus\n\
         Pos\tA\tC\tG\tT\n1\t0.1\t0.2\t0.3\t0.4\n2\t0.3333333333333333\t"
    ));
    let back = cisbp::read(&text).unwrap();
    assert_eq!(back.len(), 2);
    assert_eq!(back[0].matrix, m.matrix);
    assert_eq!(
        (
            back[0].name.as_str(),
            back[0].altname.as_deref(),
            back[0].family.as_deref(),
            back[0].organism.as_deref()
        ),
        (
            "M00001_2.00",
            Some("Gene X"),
            Some("Homeodomain"),
            Some("Mus musculus")
        )
    );
    // Metadata belongs to the next motif only.
    assert_eq!(
        (back[1].name.as_str(), back[1].altname.as_deref()),
        ("r", None)
    );
    assert_eq!(
        (back[1].alphabet.clone(), back[1].matrix.clone()),
        (Alphabet::Rna, rna.matrix.clone())
    );
    // Columns can come in any order.
    let swapped = cisbp::read("Pos T G C A\n1 0.4 0.3 0.2 0.1\n").unwrap();
    assert_eq!(swapped[0].matrix.col(0), [0.1, 0.2, 0.3, 0.4]);
    // A header without rows gives nothing.
    assert!(cisbp::read("TF x\nPos A C G T\n").unwrap().is_empty());
    // Rows before any header and keys that only start like known ones are
    // ignored.
    let m = cisbp::read("1 0.5 0.5\nTFX ignored\nTF real\nPos A C G T\n1 0.25 0.25 0.25 0.25\n")
        .unwrap();
    assert_eq!((m.len(), m[0].name.as_str()), (1, "real"));
}

#[test]
fn cisbp_errors() {
    let aa = ppm(Alphabet::Aa, &[&[0.05; 20]], "aa");
    assert_eq!(
        unsupported_msg(cisbp::write(&[aa])),
        "CIS-BP format only supports DNA and RNA motifs"
    );
    assert_eq!(
        arg_msg(cisbp::read("TF x\n1 0.25 0.25 0.25 0.25\n")),
        "no motifs found"
    );
    assert_eq!(
        arg_msg(cisbp::read("Pos A CG T\n")),
        "unexpected letter 'CG' in CIS-BP header"
    );
    assert_eq!(
        arg_msg(cisbp::read("Pos\tA\tC\tG\tT\n1\t0.5\t0.5\n")),
        "CIS-BP row '1\t0.5\t0.5' does not match its header"
    );
    assert_eq!(
        arg_msg(cisbp::read("Pos A C G X\n1 0.25 0.25 0.25 0.25\n")),
        "CIS-BP motifs must use the letters ACGT or ACGU"
    );
    assert_eq!(
        arg_msg(cisbp::read("Pos A C G T\n1 0.25 0.25 0.25 y\n")),
        "could not read a probability 'y' as a number"
    );
}

// TRANSFAC.

#[test]
fn transfac_round_trips() {
    let pcm = Motif::builder(
        Matrix::from_columns(&[[3.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 4.0]]).unwrap(),
        Alphabet::Dna,
    )
    .name("M1")
    .altname(Some("factor one".into()))
    .family(Some("bZIP".into()))
    .organism(Some("Homo sapiens".into()))
    .motif_type(Some(MotifType::Pcm))
    .build()
    .unwrap();
    let text = transfac::write(
        std::slice::from_ref(&pcm),
        &transfac::WriteOptions::default(),
    )
    .unwrap();
    assert_eq!(
        text,
        "ID M1\nNA factor one\nHC bZIP\nOS Homo sapiens\nP0\tA\tC\tG\tT\n\
         01\t3\t0\t1\t0\tA\n02\t0\t0\t0\t4\tT\nXX\n//\n"
    );
    let back = transfac::read(&text).unwrap();
    assert_eq!(back[0].matrix, pcm.matrix);
    assert_eq!(back[0].motif_type, MotifType::Pcm);
    assert_eq!(
        (
            back[0].name.as_str(),
            back[0].altname.as_deref(),
            back[0].family.as_deref(),
            back[0].organism.as_deref()
        ),
        ("M1", Some("factor one"), Some("bZIP"), Some("Homo sapiens"))
    );
    // With the accession as the name field, the identifier is the alternate
    // name.
    let opts = transfac::WriteOptions {
        name_tag: "AC".into(),
        altname_tag: "ID".into(),
    };
    let back =
        transfac::read(&transfac::write(std::slice::from_ref(&pcm), &opts).unwrap()).unwrap();
    assert_eq!(
        (back[0].name.as_str(), back[0].altname.as_deref()),
        ("M1", Some("factor one"))
    );
    // A lone NA field is the name, and probabilities are written as counts.
    let back = transfac::read("NA only\nP0 A C G T\n01 1 1 1 1 N\n//\n").unwrap();
    assert_eq!(
        (back[0].name.as_str(), back[0].altname.as_deref()),
        ("only", None)
    );
    let ppm_text =
        transfac::write(&[dna_example("p")], &transfac::WriteOptions::default()).unwrap();
    assert!(ppm_text.contains("03\t10\t20\t30\t40\tN"), "{ppm_text}");
    // RNA motifs go out with a DNA header and come back as DNA.
    let rna = motif(Alphabet::Rna, &[&[1.0, 2.0, 3.0, 4.0]], MotifType::Pcm, "r");
    let back = transfac::read(
        &transfac::write(
            std::slice::from_ref(&rna),
            &transfac::WriteOptions::default(),
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(
        (back[0].alphabet.clone(), back[0].matrix.clone()),
        (Alphabet::Dna, rna.matrix.clone())
    );
}

#[test]
fn transfac_reader_details() {
    // A release header without a matrix is skipped, columns with fewer counts
    // are scaled up to the largest, and rows that are not numbered are
    // ignored.
    let text = "VV  TRANSFAC MATRIX TABLE\nXX\n//\nID m\nP0 A C G T\n01 10 0 0 0 A\n\
                02 0 5 0 0 C\nBF something\n//\n";
    let m = transfac::read(text).unwrap();
    assert_eq!(m.len(), 1);
    // Words in the P0 line other than single letters are ignored.
    let named = transfac::read("ID m\nP0 T G C A Consensus\n01 4 3 2 1 N\n//\n").unwrap();
    assert_eq!(named[0].matrix.col(0), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(m[0].nsites, Some(10.0));
    assert_eq!(m[0].matrix.col(1), [0.0, 10.0, 0.0, 0.0]);
    assert_eq!(
        arg_msg(transfac::read("ID m\nP0 A C G T\n01 1 1 1 1\n")),
        "no motifs found (no '//' separators)"
    );
    assert_eq!(
        arg_msg(transfac::read("VV x\n//\nID y\n//\n")),
        "no motifs found"
    );
    assert_eq!(
        arg_msg(transfac::read("ID m\nP0 A C G T\n01 1 2 3\n//\n")),
        "TRANSFAC matrix row '01 1 2 3' has fewer than four counts"
    );
    assert_eq!(
        arg_msg(transfac::read("ID m\nP0 A C G T\n//\n")),
        "empty motif matrix for motif 'm'"
    );
    assert_eq!(
        arg_msg(transfac::read("ID m\nP0 A C G T\n01 1 2 3 z\n//\n")),
        "could not read a count 'z' as a number"
    );
    let aa = ppm(Alphabet::Aa, &[&[0.05; 20]], "aa");
    assert_eq!(
        unsupported_msg(transfac::write(&[aa], &transfac::WriteOptions::default())),
        "TRANSFAC format only supports DNA and RNA motifs"
    );
}

// UniPROBE.

#[test]
fn uniprobe_round_trips_and_full_files() {
    let m = dna_example("u1");
    let other = ppm(
        Alphabet::Dna,
        &[&[1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0]],
        "u2",
    );
    let text = uniprobe::write(&[m.clone(), other.clone()]).unwrap();
    assert!(
        text.starts_with("u1\nA:\t0.5\t0\t0.1\nC:\t0.25\t0\t0.2\n"),
        "{text}"
    );
    let back = uniprobe::read(&text).unwrap();
    assert_eq!(back.len(), 2);
    assert_eq!((back[0].name.as_str(), &back[0].matrix), ("u1", &m.matrix));
    assert_eq!(
        (back[1].name.as_str(), &back[1].matrix),
        ("u2", &other.matrix)
    );

    let full = "# Motif: GTAAT\tEnrichment Score: 0.49\n1\tGTAAT\t0.49\n\n\
                Enrichment score matrix\n\nA:\t-1\t1\nC:\t0\t0\nG:\t1\t-1\nT:\t0\t0\n\n\
                Probability matrix\n\nA:\t0.1\t0.7\nC:\t0.2\t0.1\nG:\t0.3\t0.1\nT:\t0.4\t0.1\n";
    let m = uniprobe::read(full).unwrap();
    assert_eq!(m.len(), 1);
    assert_eq!(m[0].name, "GTAAT");
    assert_eq!(
        m[0].extrainfo,
        [("enrichment.score".to_string(), "0.49".to_string())]
    );
    assert_eq!(m[0].matrix.col(1), [0.7, 0.1, 0.1, 0.1]);
}

#[test]
fn uniprobe_errors() {
    let rna = ppm(Alphabet::Rna, &[&[0.25; 4]], "r");
    assert_eq!(
        unsupported_msg(uniprobe::write(&[rna])),
        "UniPROBE format only supports DNA motifs"
    );
    assert_eq!(
        arg_msg(uniprobe::read("name\n0.1 0.2\n")),
        "no motifs found"
    );
    assert_eq!(
        arg_msg(uniprobe::read("# Motif: m\n1 m 0.5\n")),
        "every UniPROBE motif needs a probability matrix"
    );
    assert_eq!(
        arg_msg(uniprobe::read("m\nA: 0.5\nG: 0.5\nC: 0\nT: 0\n")),
        "motif 'm': expected a row starting with 'C:'"
    );
    assert_eq!(
        arg_msg(uniprobe::read("m\nA: 0.5\nC: 0.5 0\nG: 0\nT: 0\n")),
        "motif 'm' has rows of different lengths"
    );
    assert_eq!(
        arg_msg(uniprobe::read("m\nA: 0.5\nC: 0.5\nG: 0\n")),
        "motif 'm': expected a row starting with 'T:'"
    );
    assert_eq!(
        arg_msg(uniprobe::read("m\nA: 0.5\nC: q\nG: 0\nT: 0\n")),
        "could not read a probability 'q' as a number"
    );
}

// Plain matrices.

#[test]
fn matrix_round_trips_every_layout() {
    let custom = Alphabet::custom("xyz".chars()).unwrap();
    let motifs = vec![
        motif(
            custom.clone(),
            &[&[1.0, 2.0, 3.0], &[4.0, 1.0, 1.0]],
            MotifType::Pcm,
            "first",
        ),
        motif(
            custom.clone(),
            &[&[0.0, 0.0, 6.0]],
            MotifType::Pcm,
            "second",
        ),
    ];
    let layouts = [
        (matrix::Headers::FirstLine, Some(String::new())),
        (matrix::Headers::FirstLine, Some("//".to_string())),
        (matrix::Headers::None, Some(String::new())),
        (matrix::Headers::Prefix(">".into()), None),
        (
            matrix::Headers::Prefix("motif: ".into()),
            Some(String::new()),
        ),
    ];
    for positions in [matrix::Positions::Columns, matrix::Positions::Rows] {
        for rownames in [false, true] {
            for (headers, sep) in &layouts {
                let w = matrix::WriteOptions {
                    positions,
                    rownames,
                    headers: headers.clone(),
                    sep: sep.clone(),
                    digits: 2,
                    ..Default::default()
                };
                let r = matrix::ReadOptions {
                    positions,
                    rownames,
                    headers: headers.clone(),
                    sep: sep.clone(),
                    alphabet: custom.clone(),
                    motif_type: Some(MotifType::Pcm),
                    ..Default::default()
                };
                let ctx = format!("{positions:?} {rownames} {headers:?} {sep:?}");
                let text = matrix::write(&motifs, &w).unwrap();
                let back = matrix::read(&text, &r).unwrap();
                assert_eq!(back.len(), 2, "{ctx}\n{text}");
                for (a, b) in back.iter().zip(&motifs) {
                    assert_eq!(a.matrix, b.matrix, "{ctx}");
                    if *headers != matrix::Headers::None {
                        assert_eq!(a.name, b.name, "{ctx}");
                    }
                }
            }
        }
    }
    let text = matrix::write(
        &motifs[1..],
        &matrix::WriteOptions {
            rownames: true,
            digits: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(text, "second\nx\t0.0\ny\t0.0\nz\t6.0\n\n");
    let text = matrix::write(
        &motifs[..1],
        &matrix::WriteOptions {
            positions: matrix::Positions::Rows,
            rownames: true,
            headers: matrix::Headers::None,
            sep: None,
            digits: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(text, "x\ty\tz\n1\t2\t3\n4\t1\t1\n");
}

#[test]
fn matrix_types_and_infinities() {
    let p = ppm(Alphabet::Dna, &[&[0.5, 0.5, 0.0, 0.0]], "p");
    // A zero probability with no pseudocount is a weight of minus infinity.
    let text = matrix::write(
        std::slice::from_ref(&p),
        &matrix::WriteOptions {
            motif_type: Some(MotifType::Pwm),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(text, "p\n1.000000\n1.000000\n-Inf\n-Inf\n\n");
    let back = matrix::read(
        &text,
        &matrix::ReadOptions {
            motif_type: Some(MotifType::Pwm),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        back[0].matrix.col(0),
        [1.0, 1.0, f64::NEG_INFINITY, f64::NEG_INFINITY]
    );
    // Without a type it is detected from the values.
    let back = matrix::read("m\n0.1\n0.2\n0.3\n0.4\n", &matrix::ReadOptions::default()).unwrap();
    assert_eq!(back[0].motif_type, MotifType::Ppm);
    let back = matrix::read("m\n1\n2\n3\n4\n", &matrix::ReadOptions::default()).unwrap();
    assert_eq!(back[0].motif_type, MotifType::Pcm);
    // Comments are skipped, and with no separator only exactly empty lines
    // end a motif.
    let opts = matrix::ReadOptions {
        comment: Some("#".into()),
        sep: None,
        ..Default::default()
    };
    let back = matrix::read(
        "# note\na\n1 1\n# more\n1 1\n1 1\n1 1\n\nb\n3\n3\n3\n3\n",
        &opts,
    )
    .unwrap();
    assert_eq!(
        back.iter().map(|m| m.name.as_str()).collect::<Vec<_>>(),
        ["a", "b"]
    );
    assert_eq!(back[0].width(), 2);
}

#[test]
fn matrix_reader_errors() {
    let opts = matrix::ReadOptions::default();
    assert_eq!(arg_msg(matrix::read("\n\n", &opts)), "no motifs found");
    assert_eq!(
        arg_msg(matrix::read(
            "x\n",
            &matrix::ReadOptions {
                headers: matrix::Headers::Prefix(">".into()),
                ..Default::default()
            }
        )),
        "no motifs found"
    );
    assert_eq!(
        arg_msg(matrix::read("m\n1 2\n3\n4 4\n5 5\n", &opts)),
        "motif 'm' has an empty or ragged matrix"
    );
    assert_eq!(
        arg_msg(matrix::read("m\n\n", &opts)),
        "motif 'm' has an empty or ragged matrix"
    );
    assert_eq!(
        arg_msg(matrix::read("m\n", &opts)),
        "motif 'm' has an empty or ragged matrix"
    );
    assert_eq!(
        arg_msg(matrix::read(
            "1\n2\n3\n4\n",
            &matrix::ReadOptions {
                headers: matrix::Headers::None,
                rownames: true,
                ..Default::default()
            }
        )),
        "motif '' has an empty or ragged matrix"
    );
    assert_eq!(
        arg_msg(matrix::read("m\n1\n2\nthree\n4\n", &opts)),
        "could not read a matrix value 'three' as a number"
    );
    // The wrong number of rows for the alphabet is caught when building.
    assert!(matrix::read("m\n1\n2\n3\n", &opts).is_err());
}
