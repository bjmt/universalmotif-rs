//! Behaviour that differs from the R package on purpose, checked directly
//! rather than against reference values (see DIVERGENCES.md).

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
        Some(Alphabet::parse("XYZ").unwrap()),
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
