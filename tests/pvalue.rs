//! P-values and score thresholds against `motif_pvalue(method = "dynamic")`.

mod common;

use common::*;
use universalmotif::pvalue::{PvalueOptions, motif_pvalue, pvalue_to_score};

#[test]
fn pvalues_and_scores_match_r() {
    for (i, case) in load("pvalue.json").iter().enumerate() {
        let m = motif_from_record(&case["input"]);
        let a = args(&case["args"]);
        let ctx = format!("case {i} ({}, {:?})", m.name, a.get("bkg"));
        let opts = PvalueOptions {
            bkg: a.get("bkg").map(|b| vec![nums(b)]).unwrap_or_default(),
            ..Default::default()
        };
        let motifs = [m];

        let q: Vec<(usize, f64)> = nums(&case["scores"]).into_iter().map(|s| (0, s)).collect();
        let got = motif_pvalue(&motifs, &q, &opts).unwrap();
        assert_close_vec(&got, &nums(&case["pvalues"]), &format!("{ctx}: P-values"));

        let q: Vec<(usize, f64)> = nums(&case["pvalue_in"])
            .into_iter()
            .map(|p| (0, p))
            .collect();
        let got = pvalue_to_score(&motifs, &q, &opts).unwrap();
        assert_close_vec(&got, &nums(&case["score_out"]), &format!("{ctx}: scores"));
    }
}

/// Values hard-coded in universalmotif's own tests
/// (tests/testthat/test_threshold_contract.R and test_background_contract.R
/// check these against brute-force enumeration).
#[test]
fn scores_are_thresholds_on_the_integer_grid() {
    use universalmotif::{Alphabet, Matrix, Motif, MotifType};
    // Two positions, PWM values on the grid: brute force gives the exact tail.
    let pwm = Matrix::from_columns(&[[1.0, 0.0, -1.0, -2.0], [0.5, 0.5, -0.5, -1.0]]).unwrap();
    let m = Motif::builder(pwm, Alphabet::Dna)
        .motif_type(Some(MotifType::Pwm))
        .build()
        .unwrap();
    let opts = PvalueOptions::default();
    let all: Vec<f64> = [1.0, 0.0, -1.0, -2.0]
        .iter()
        .flat_map(|a| [0.5, 0.5, -0.5, -1.0].map(|b| a + b))
        .collect();
    for s in [-3.0, -1.5, -0.5, 0.0, 0.5, 1.0, 1.5, 1.501] {
        let want = all.iter().filter(|&&x| x >= s - 1e-9).count() as f64 / 16.0;
        let got = motif_pvalue(std::slice::from_ref(&m), &[(0, s)], &opts).unwrap()[0];
        assert!(
            (got - want).abs() < 1e-12,
            "score {s}: got {got}, want {want}"
        );
    }
    // Two sites score 1.5 and the next best score is 0.5, so the lowest grid
    // score reached with probability 2/16 is just above 0.5.
    let s = pvalue_to_score(std::slice::from_ref(&m), &[(0, 2.0 / 16.0)], &opts).unwrap()[0];
    assert_eq!(s, 0.501);
    let s = pvalue_to_score(std::slice::from_ref(&m), &[(0, 1e-9)], &opts).unwrap()[0];
    assert_eq!(s, f64::INFINITY);
}
