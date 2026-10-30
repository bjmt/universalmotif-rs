//! Behaviour that differs from the R package on purpose, checked directly
//! rather than against reference values.

mod common;

use common::*;
use universalmotif::compare::{CompareOptions, compare_matrix};
use universalmotif::convert::{convert_type, ppm_to_pcm};
use universalmotif::coocc::{CooccHit, CooccOptions, motif_coocc};
use universalmotif::create::{CreateOptions, from_matrix, from_strings};
use universalmotif::io::meme;
use universalmotif::ops::{normalize, round_motif};
use universalmotif::{Alphabet, Matrix, MotifType};

fn opts() -> CreateOptions {
    CreateOptions::default()
}

#[test]
fn counts_are_spread_evenly() {
    assert_eq!(ppm_to_pcm(&[0.25; 4], 50.0), vec![13.0, 13.0, 12.0, 12.0]);
    // A shortfall still goes to the largest count, keeping pure columns pure.
    let mut aa = vec![0.001; 20];
    aa[0] = 0.981;
    let counts = ppm_to_pcm(&aa, 100.0);
    assert_eq!(counts[0], 100.0);
    assert!(counts[1..].iter().all(|&c| c == 0.0));
    // The R package makes a negative count here and fails.
    let o = CreateOptions {
        nsites: Some(50.0),
        ..opts()
    };
    let m = from_strings(&["BZJX"], Some(Alphabet::Aa), &o).unwrap();
    assert!(m.matrix.as_slice().iter().all(|&v| v >= 0.0));
}

#[test]
fn derived_values_are_recomputed() {
    let m = from_strings(&["TTGACC", "TTGACT", "ATGACC", "TTGAAC"], None, &opts()).unwrap();
    let r = round_motif(&convert_type(&m, MotifType::Ppm, None, false).unwrap(), 0.3).unwrap();
    assert!(close(r.icscore, own_icscore(&r)));
    let n = normalize(&m.clone().into_parts().pipe_zero_pseudo()).unwrap();
    assert!(close(n.icscore, own_icscore(&n)));
    let p = convert_type(&m, MotifType::Pwm, Some(5.0), false).unwrap();
    assert_eq!(p.pseudocount, 5.0);
    assert!(p.icscore < m.icscore, "a larger pseudocount lowers the IC");
}

trait ZeroPseudo {
    fn pipe_zero_pseudo(self) -> universalmotif::Motif;
}

impl ZeroPseudo for universalmotif::MotifParts {
    fn pipe_zero_pseudo(mut self) -> universalmotif::Motif {
        self.pseudocount = 0.0;
        universalmotif::Motif::from_parts(self).unwrap()
    }
}

#[test]
fn cwms_have_no_total_ic() {
    let cols = [
        vec![0.5, -0.2, 0.1, 0.0],
        vec![0.0, 0.0, 0.3, -0.1],
        vec![0.1, 0.6, 0.0, 0.0],
    ];
    let o = CreateOptions {
        motif_type: MotifType::Cwm,
        ..opts()
    };
    let m = from_matrix(
        Matrix::from_columns(&cols).unwrap(),
        None,
        Some(Alphabet::Dna),
        &o,
    )
    .unwrap();
    assert!(m.icscore.is_nan());
    assert_eq!(m, m.clone(), "a CWM equals itself");
    assert!(!m.to_string().contains("Total IC"));
    let ppm = convert_type(&m, MotifType::Ppm, None, false).unwrap();
    assert!(close(ppm.icscore, own_icscore(&ppm)));
}

#[test]
fn rounded_probabilities_are_guessed_as_ppm() {
    let cols = [vec![0.33, 0.33, 0.33, 0.0], vec![0.7, 0.1, 0.1, 0.1]];
    let m = from_matrix(
        Matrix::from_columns(&cols).unwrap(),
        None,
        Some(Alphabet::Dna),
        &opts(),
    )
    .unwrap();
    assert_eq!(m.motif_type, MotifType::Ppm);
    assert!(m.matrix.col_sums().iter().all(|s| (s - 1.0).abs() < 1e-12));
}

#[test]
fn comparison_scores_are_symmetric() {
    let seqs = [
        "ACGTNN", "NNACGT", "ACGT", "NACGTN", "TTACGTAA", "CANNTG", "ACNNNGT",
    ];
    let ms: Vec<_> = seqs
        .iter()
        .map(|s| from_strings(&[s], None, &opts()).unwrap())
        .collect();
    let o = CompareOptions {
        min_overlap: 1,
        rc: false,
        ..Default::default()
    };
    let s = compare_matrix(&ms, &o, false).unwrap().scores;
    for i in 0..ms.len() {
        for j in 0..ms.len() {
            assert!(close(s[i][j], s[j][i]), "{} vs {}", seqs[i], seqs[j]);
        }
    }
}

#[test]
fn self_pairs_count_repeated_sites() {
    // Motif 0 has two sites in sequence 0, one in sequence 1, and a site on
    // both strands at one position in sequence 2.
    let hit = |motif, sequence, start| CooccHit {
        motif,
        sequence,
        start,
    };
    let hits = [
        hit(0, 0, 10),
        hit(0, 0, 40),
        hit(0, 1, 5),
        hit(0, 2, 7),
        hit(0, 2, 7),
        hit(1, 0, 20),
    ];
    let o = CooccOptions {
        self_pairs: true,
        max_distance: Some(100),
        ..Default::default()
    };
    let pairs = motif_coocc(&hits, 2, 4, &o).unwrap();
    let own = pairs.iter().find(|p| p.a == 0 && p.b == 0).unwrap();
    assert_eq!((own.both, own.a_only, own.neither), (1, 2, 1));
    assert!(own.pvalue.is_none() && own.odds_ratio.is_none());
    assert_eq!(own.median_distance, Some(30.0));
}

#[test]
fn meme_files_are_tidy() {
    let m = from_strings(&["TTGACC", "TTGACT"], None, &opts()).unwrap();
    let text = meme::write(&[m], &meme::WriteOptions::default()).unwrap();
    assert!(text.lines().all(|l| l == l.trim_end()), "{text}");
    assert!(!text.contains("E="), "{text}");
    let back = meme::read(&text, &meme::ReadOptions::default()).unwrap();
    assert_eq!(back[0].eval, None);
    let def = meme::AlphabetDefinition {
        core: "ACGT".into(),
        complements: Some("TGCA".into()),
        ..Default::default()
    };
    let alph = meme::write_alphabet(&def).unwrap();
    assert_eq!(alph.matches(" ~ ").count(), 2, "{alph}");
}

#[test]
fn custom_alphabet_columns_are_numbered() {
    let m = from_strings(
        &["XYZ", "XYY"],
        Some("XYZ".parse::<Alphabet>().unwrap()),
        &opts(),
    )
    .unwrap();
    let text = m.to_string();
    assert!(text.contains("  1") && !text.contains("[,"), "{text}");
}

#[test]
fn log_pvalues_survive_underflow() {
    use universalmotif::stats::{fisher_greater, ln_fisher_greater, ln_p_adjust_bh, p_adjust_bh};
    for &(a, b, c, d) in &[(10, 2, 5, 20), (40, 3, 10, 60), (1, 1, 1, 1)] {
        let p = fisher_greater(a, b, c, d);
        assert!((ln_fisher_greater(a, b, c, d) - p.ln()).abs() < 1e-9 * p.ln().abs().max(1.0));
    }
    // 1 / choose(6000, 3000), about 1e-1804.
    let ln_p = ln_fisher_greater(3000, 0, 0, 3000);
    assert!(fisher_greater(3000, 0, 0, 3000) == 0.0);
    assert!(
        (ln_p / std::f64::consts::LN_10 + 1804.0).abs() < 1.0,
        "{ln_p}"
    );
    let p = [0.01, 0.04, 0.03, 0.2];
    let q = p_adjust_bh(&p);
    let ln_q = ln_p_adjust_bh(&p.map(f64::ln), 4.0);
    for (q, l) in q.iter().zip(ln_q) {
        assert!((q.ln() - l).abs() < 1e-12);
    }
}

#[test]
fn complete_scores_cover_the_whole_query() {
    use universalmotif::compare::compare_motifs;
    let o = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let make = |s: &str| from_strings(&[s], None, &o).unwrap();
    // A long query, a copy of it, and a short motif matching only its middle.
    let mut motifs = vec![
        make("TTGACGTCAGCCATG"),
        make("TTGACGTCAGCCATG"),
        make("TGACGTCA"),
    ];
    for s in [
        "CACGTGAA", "GATAAGCT", "CCAATCAG", "AAACCCTA", "TATAAAAG", "GGGGCGGG",
    ] {
        motifs.push(make(s));
    }
    let partial = compare_motifs(&motifs, &[0], &CompareOptions::default()).unwrap();
    let complete = compare_motifs(
        &motifs,
        &[0],
        &CompareOptions {
            complete: true,
            ..Default::default()
        },
    )
    .unwrap();
    let p = |rows: &[universalmotif::compare::Comparison], t: usize| {
        rows.iter().find(|c| c.target == t).unwrap().pvalue
    };
    // Scoring the whole query sets the full copy much further apart from the
    // short motif that matches only part of it than scoring the overlap does.
    let gap = |rows: &[universalmotif::compare::Comparison]| p(rows, 2) / p(rows, 1);
    assert!(gap(&complete) > 5.0 * gap(&partial));
    assert!(p(&complete, 1) * 1e4 < p(&complete, 2));
    // The overlap and its mean score are reported either way.
    let short = complete.iter().find(|c| c.target == 2).unwrap();
    assert_eq!(short.overlap, 8);
    assert!(short.score > 0.9);
}

#[test]
fn score_match_agrees_with_scanning() {
    use universalmotif::scan::{ScanOptions, scan_sequences};
    use universalmotif::scoring::score_match;
    // The R package adds the rounded entries as floats and truncates the
    // total again, which can be 0.001 off the score scanning reports.
    let o = CreateOptions {
        pseudocount: 1.0,
        ..opts()
    };
    let motifs = [
        from_strings(&["TGACGTCA", "TGACGTCT", "TTACGTCA", "TGACGACA"], None, &o).unwrap(),
        from_strings(&["TATAAAT", "TATATAT", "TATAAAA"], None, &o).unwrap(),
        from_strings(&["GGCCAATC", "AGCCAATG", "GGCCAATC"], None, &o).unwrap(),
    ];
    let seqs = [
        "TAACGTGGACGTCATTTATATAAGGGCCAATCGAT",
        "GGGTGACGTAAACCATATAAAT",
    ];
    let opts = ScanOptions {
        pvalue: 0.5,
        rc: false,
        nthreads: 1,
        ..Default::default()
    };
    let mut checked = 0;
    for m in &motifs {
        for h in scan_sequences(std::slice::from_ref(m), &seqs, &opts).unwrap() {
            let word = &seqs[h.sequence][h.start..h.start + m.width()];
            assert_eq!(
                score_match(m, &[word], false).unwrap()[0],
                h.score,
                "{word}"
            );
            checked += 1;
        }
    }
    assert!(checked > 10);
}

#[test]
fn bad_matrices_and_frequencies_are_rejected() {
    use universalmotif::shuffle::create_sequences;
    // Frequencies that are all zero would otherwise give only the first letter.
    let e = create_sequences(&['A', 'C', 'G', 'T'], 2, 10, Some(&[0.0; 4]), 1, 1).unwrap_err();
    assert!(e.to_string().contains("above zero"), "{e}");
    // A count matrix with an empty column is not a motif.
    let zero = Matrix::from_vec(4, 2, vec![0.0; 8]).unwrap();
    let e = from_matrix(zero, None, Some(Alphabet::Dna), &opts()).unwrap_err();
    assert!(e.to_string().contains("at least one count"), "{e}");
    // Whole numbers are read as counts, and uneven columns say so.
    let ints = Matrix::from_vec(4, 2, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]).unwrap();
    let e = from_matrix(ints, None, Some(Alphabet::Dna), &opts()).unwrap_err();
    assert!(e.to_string().contains("read as counts"), "{e}");
}

fn dna(site: &str, strand: universalmotif::Strand) -> universalmotif::Motif {
    let o = CreateOptions {
        strand,
        nsites: Some(20.0),
        ..opts()
    };
    from_strings(&[site], Some(Alphabet::Dna), &o).unwrap()
}

#[test]
fn scans_honour_the_motif_strand() {
    use universalmotif::Strand;
    use universalmotif::scan::{ScanOptions, scan_table};
    // The site only appears as its reverse complement.
    let seq = format!("{}GGGTTTTCC{}", "A".repeat(20), "C".repeat(20));
    let strands = |strand: Strand, rc: bool| {
        let o = ScanOptions {
            pvalue: 1e-3,
            rc,
            ..Default::default()
        };
        let t = scan_table(&[dna("GGAAAACCC", strand)], &[seq.as_str()], &o).unwrap();
        t.minus
    };
    assert_eq!(strands(Strand::Both, true), [true]);
    assert!(strands(Strand::Both, false).is_empty());
    // A one-stranded motif ignores rc.
    for rc in [true, false] {
        assert!(strands(Strand::Plus, rc).is_empty(), "rc {rc}");
        assert_eq!(strands(Strand::Minus, rc), [true], "rc {rc}");
    }
    // The same holds for gapped motifs.
    let gapped = |strand: Strand| {
        let m = universalmotif::ops::add_gap(&dna("GGAAAACCC", strand), Some(&[4]), &[0], &[1])
            .unwrap();
        let o = ScanOptions {
            pvalue: 1e-3,
            ..Default::default()
        };
        scan_table(&[m], &[seq.as_str()], &o).unwrap().minus
    };
    assert!(gapped(Strand::Plus).is_empty());
    assert!(gapped(Strand::Minus).iter().all(|&m| m));
    assert!(!gapped(Strand::Minus).is_empty());
}

#[test]
fn enrichment_counts_positions_on_the_strands_scanned() {
    use universalmotif::Strand;
    use universalmotif::enrich::{EnrichOptions, EnrichTest, enrich_motifs};
    let targets = vec![format!("{}GGAAAACCC{}", "A".repeat(20), "C".repeat(20)); 5];
    let bkg = vec!["ACGT".repeat(12); 5];
    let o = EnrichOptions {
        test: EnrichTest::Sites,
        pvalue: 1e-3,
        ..Default::default()
    };
    let both = enrich_motifs(&[dna("GGAAAACCC", Strand::Both)], &targets, Some(&bkg), &o).unwrap();
    let plus = enrich_motifs(&[dna("GGAAAACCC", Strand::Plus)], &targets, Some(&bkg), &o).unwrap();
    assert_eq!(both[0].target_site_hits, plus[0].target_site_hits);
    // The same hits out of half the positions: 41 and 40 windows in each
    // target and background, with a pseudocount of 1 in every cell.
    use universalmotif::stats::fisher_greater;
    let want = |strands: u64| fisher_greater(6, 1, 41 * 5 * strands - 5 + 1, 40 * 5 * strands + 1);
    assert!((both[0].pvalue - want(2)).abs() < 1e-12);
    assert!((plus[0].pvalue - want(1)).abs() < 1e-12);
}

#[test]
fn palindromes_count_once_without_overlaps() {
    use universalmotif::Strand;
    use universalmotif::scan::{NoOverlaps, ScanOptions, scan_table};
    let seq = format!("{}CCACGTGG{}", "A".repeat(20), "T".repeat(20));
    let m = dna("CCACGTGG", Strand::Both);
    let o = ScanOptions {
        pvalue: 1e-3,
        ..Default::default()
    };
    assert_eq!(
        scan_table(std::slice::from_ref(&m), &[seq.as_str()], &o)
            .unwrap()
            .len(),
        2
    );
    let once = ScanOptions {
        no_overlaps: Some(NoOverlaps::default()),
        ..o.clone()
    };
    assert!(NoOverlaps::default().ignore_strand);
    assert_eq!(
        scan_table(std::slice::from_ref(&m), &[seq.as_str()], &once)
            .unwrap()
            .len(),
        1
    );
    let by_strand = ScanOptions {
        no_overlaps: Some(NoOverlaps {
            ignore_strand: false,
            ..Default::default()
        }),
        ..o
    };
    assert_eq!(
        scan_table(&[m], &[seq.as_str()], &by_strand).unwrap().len(),
        2
    );
}

#[test]
fn gapped_pvalues_are_corrected_by_default() {
    use universalmotif::Strand;
    use universalmotif::scan::{ScanOptions, scan_table};
    assert!(ScanOptions::default().correct_gaps);
    let m = universalmotif::ops::add_gap(&dna("TTGACCAAT", Strand::Both), Some(&[3]), &[0], &[3])
        .unwrap();
    let seq = format!(
        "{}TTGACCAAT{}TTGAACCAAT{}",
        "A".repeat(30),
        "C".repeat(30),
        "G".repeat(30)
    );
    let base = ScanOptions {
        pvalue: 1e-2,
        rc: false,
        ..Default::default()
    };
    let plain = ScanOptions {
        correct_gaps: false,
        pvalue: 1e-2 / 4.0,
        ..base.clone()
    };
    let corrected = scan_table(std::slice::from_ref(&m), &[seq.as_str()], &base).unwrap();
    let raw = scan_table(std::slice::from_ref(&m), &[seq.as_str()], &plain).unwrap();
    // Four gap sizes, so four times the P-value at a quarter of the cutoff.
    assert!(!corrected.is_empty());
    assert_eq!(corrected.start, raw.start);
    for (c, r) in corrected.pvalue.iter().zip(&raw.pvalue) {
        assert!((c / r - 4.0).abs() < 1e-9 || *c == 1.0, "{c} {r}");
    }
    // Motifs without gaps are untouched.
    let ungapped = dna("TTGACCAAT", Strand::Both);
    let a = scan_table(std::slice::from_ref(&ungapped), &[seq.as_str()], &base).unwrap();
    let b = scan_table(
        &[ungapped],
        &[seq.as_str()],
        &ScanOptions {
            correct_gaps: false,
            ..base
        },
    )
    .unwrap();
    assert_eq!(a, b);
}

#[test]
fn complete_linkage_keeps_chains_apart() {
    use universalmotif::create::random;
    use universalmotif::merge::{MergeLinkage, MergeSimilarOptions, merge_similar};
    let sites = [
        "GATTACAGCATGCA",
        "CAGCATGCATTGAC",
        "ATGCATTGACCGTA",
        "TTGACCGTAGGTCA",
    ];
    let mut motifs: Vec<_> = sites
        .iter()
        .map(|s| dna(s, universalmotif::Strand::Both))
        .collect();
    motifs.extend((0..8).map(|i| random(10, Alphabet::Dna, &opts(), i).unwrap()));
    let run = |linkage| {
        let o = MergeSimilarOptions {
            linkage,
            ..Default::default()
        };
        merge_similar(&motifs, &o).unwrap()
    };
    let single = run(MergeLinkage::Single);
    assert_eq!(&single.clusters[..5], [0, 0, 0, 0, 1]);
    assert_eq!(single.warnings.len(), 1);
    assert!(
        single.warnings[0].contains("cluster 0 holds 4 of the 12 motifs"),
        "{}",
        single.warnings[0]
    );
    let complete = run(MergeLinkage::Complete);
    assert_eq!(&complete.clusters[..5], [0, 1, 1, 1, 2]);
    assert!(complete.warnings.is_empty());
    assert_eq!(complete.motifs.len(), 10);
    // Every pair within a complete-linkage cluster is linked, so the first
    // and last motifs of the chain are not.
    assert_ne!(complete.clusters[0], complete.clusters[3]);
    assert_eq!(
        "complete".parse::<MergeLinkage>().unwrap(),
        MergeLinkage::Complete
    );
    assert!("average".parse::<MergeLinkage>().is_err());
}

#[test]
fn fasta_names_are_trimmed_to_their_first_word() {
    use universalmotif::io::fasta::{first_word, trim_names};
    assert_eq!(first_word("chr1 some description"), "chr1");
    assert_eq!(first_word("chr1\tx"), "chr1");
    assert_eq!(first_word(""), "");
    let names = vec!["a x".to_string(), "b".into()];
    assert_eq!(trim_names(names).unwrap(), ["a", "b"]);
    let e = trim_names(vec!["chr1 a".into(), "chr1 b".into()]).unwrap_err();
    assert!(
        e.to_string()
            .contains("the sequence name 'chr1' appears more than once"),
        "{e}"
    );
}

#[test]
fn shuffles_keep_dinucleotides_by_default() {
    use universalmotif::shuffle::{ShuffleOptions, shuffle_sequences};
    let seq = "ACGCGTTAGCATCGATCGGGCATTACGCGATATCGCGATCGATCGCGCTAGCTAGG";
    let out = shuffle_sequences(&[seq], &ShuffleOptions::default()).unwrap();
    let pairs = |s: &str| {
        let mut v: Vec<String> = (0..s.len() - 1).map(|i| s[i..i + 2].to_string()).collect();
        v.sort_unstable();
        v
    };
    assert_ne!(out[0], seq);
    assert_eq!(pairs(&out[0]), pairs(seq));
}

#[test]
fn motifs_found_on_one_strand_are_for_that_strand() {
    use universalmotif::Strand;
    use universalmotif::finder::{FinderOptions, motif_finder};
    let site = "TTGACCAT";
    let seqs: Vec<String> = (0..20)
        .map(|i| {
            let bkg =
                universalmotif::shuffle::create_sequences(&['A', 'C', 'G', 'T'], 1, 92, None, i, 1)
                    .unwrap()
                    .remove(0);
            format!("{}{site}{}", &bkg[..46], &bkg[46..])
        })
        .collect();
    for (rc, strand) in [(true, Strand::Both), (false, Strand::Plus)] {
        let o = FinderOptions {
            min_width: 6,
            max_width: 8,
            rc,
            ..Default::default()
        };
        let found = motif_finder(&seqs, None, &o).unwrap();
        assert!(!found.is_empty(), "rc {rc}");
        assert!(found.iter().all(|f| f.motif.strand == strand), "rc {rc}");
    }
}
