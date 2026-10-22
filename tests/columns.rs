//! Column conversion functions against universalmotif's C++ helpers.

mod common;

use common::*;
use universalmotif::convert::*;

#[test]
fn column_functions_match_r() {
    let cases = load("columns.json");
    let mut counts = std::collections::BTreeMap::new();
    for (i, case) in cases.iter().enumerate() {
        let f = case["fn"].as_str().unwrap();
        let a = &case["args"];
        let ctx = format!("case {i} ({f}, args {a})");
        *counts.entry(f).or_insert(0) += 1;
        match f {
            "pcm_to_ppm" => assert_close_vec(
                &pcm_to_ppm(&nums(&a["pos"]), num(&a["pseudocount"])),
                &nums(&case["out"]),
                &ctx,
            ),
            "ppm_to_pcm" => {
                let (pos, n) = (nums(&a["pos"]), num(&a["nsites"]));
                let got = ppm_to_pcm(&pos, n);
                if r_rounding_differs(&pos, n) {
                    // R adds the whole rounding error to the largest count.
                    assert_eq!(got.iter().sum::<f64>(), n, "{ctx}: sum");
                    for (g, p) in got.iter().zip(&pos) {
                        assert!((g - p * n).abs() < 1.0, "{ctx}: {g} for {p}");
                    }
                } else {
                    assert_close_vec(&got, &nums(&case["out"]), &ctx);
                }
            }
            "ppm_to_pwm" => assert_close_vec(
                &ppm_to_pwm(
                    &nums(&a["pos"]),
                    &nums(&a["bkg"]),
                    num(&a["pseudocount"]),
                    num(&a["nsites"]),
                ),
                &nums(&case["out"]),
                &ctx,
            ),
            "pwm_to_ppm" => assert_close_vec(
                &pwm_to_ppm(&nums(&a["pos"]), &nums(&a["bkg"])),
                &nums(&case["out"]),
                &ctx,
            ),
            "ppm_to_icm" => {
                let pos = nums(&a["pos"]);
                // The R package gives NaN for relative entropy where a
                // probability is 0; the port gives 0 (see DIVERGENCES.md).
                let want: Vec<f64> = nums(&case["out"])
                    .iter()
                    .zip(&pos)
                    .map(|(&w, &p)| if w.is_nan() && p == 0.0 { 0.0 } else { w })
                    .collect();
                assert_close_vec(
                    &ppm_to_icm(
                        &pos,
                        &nums(&a["bkg"]),
                        a["relative_entropy"].as_bool().unwrap(),
                    ),
                    &want,
                    &ctx,
                )
            }
            "icm_to_ppm" => {
                assert_close_vec(&icm_to_ppm(&nums(&a["pos"])), &nums(&case["out"]), &ctx)
            }
            "position_icscore" => {
                let got = position_icscore(
                    &nums(&a["pos"]),
                    &nums(&a["bkg"]),
                    motif_type(a["type"].as_str().unwrap()),
                    num(&a["pseudocount"]),
                    num(&a["nsites"]),
                    a["relative_entropy"].as_bool().unwrap(),
                );
                let want = num(&case["out"]);
                if want.is_nan() && a["relative_entropy"].as_bool().unwrap() {
                    // The R package's NaN for a zero probability (see
                    // DIVERGENCES.md); the port's value is finite.
                    assert!(got.is_finite(), "{ctx}: got {got}");
                } else {
                    assert!(close(got, want), "{ctx}: got {got}, want {want}");
                }
            }
            "get_consensus" => {
                let got = consensus_letter(
                    &nums(&a["pos"]),
                    a["alphabet"] == "RNA",
                    motif_type(a["type"].as_str().unwrap()),
                    num(&a["pseudocount"]),
                );
                assert_eq!(got.to_string(), case["out"].as_str().unwrap(), "{ctx}");
            }
            "get_consensus_aa" => {
                let got = consensus_letter_aa(
                    &nums(&a["pos"]),
                    motif_type(a["type"].as_str().unwrap()),
                    num(&a["pseudocount"]),
                );
                assert_eq!(got.to_string(), case["out"].as_str().unwrap(), "{ctx}");
            }
            "consensus_to_ppm" => {
                let l = a["letter"].as_str().unwrap().chars().next().unwrap();
                assert_close_vec(&consensus_to_ppm(l).unwrap(), &nums(&case["out"]), &ctx);
            }
            "consensus_to_ppm_aa" => {
                let l = a["letter"].as_str().unwrap().chars().next().unwrap();
                assert_close_vec(&consensus_to_ppm_aa(l).unwrap(), &nums(&case["out"]), &ctx);
            }
            _ => panic!("unknown function {f}"),
        }
    }
    assert!(
        counts.len() == 11,
        "expected all column functions, got {counts:?}"
    );
}

/// Values hard-coded in universalmotif's own tests (tests/testthat/test_utils.R).
#[test]
fn values_from_universalmotif_tests() {
    use universalmotif::MotifType;
    let rnd = |v: Vec<f64>, d: i32| {
        let f = 10f64.powi(d);
        v.iter().map(|x| (x * f).round() / f).collect::<Vec<_>>()
    };
    assert_eq!(pcm_to_ppm(&[1.0, 1.0], 0.0), vec![0.5, 0.5]);
    assert_eq!(pcm_to_ppm(&[1.0, 0.0], 1.0), vec![0.75, 0.25]);
    assert_eq!(ppm_to_pcm(&[0.25, 0.75], 10.0), vec![3.0, 7.0]);
    assert_eq!(
        rnd(ppm_to_pwm(&[0.25, 0.75], &[], 0.0, 100.0), 3),
        vec![-1.0, 0.585]
    );
    assert_eq!(
        rnd(ppm_to_pwm(&[0.25, 0.75], &[0.4, 0.6], 0.0, 100.0), 3),
        vec![-0.678, 0.322]
    );
    assert_eq!(rnd(pwm_to_ppm(&[-1.0, 0.5], &[]), 4), vec![0.2612, 0.7388]);
    assert_eq!(
        rnd(ppm_to_icm(&[0.25, 0.75], &[], false), 4),
        vec![0.0472, 0.1415]
    );
    assert_eq!(
        rnd(ppm_to_icm(&[0.25, 0.75], &[], true), 4),
        vec![0.0, 0.4387]
    );
    let ic = position_icscore(&[0.25, 0.75], &[], MotifType::Ppm, 1.0, 100.0, false);
    assert_eq!(rnd(vec![ic], 3), vec![0.185]);
    assert_eq!(
        rnd(icm_to_ppm(&[0.5, 0.0, 0.3]), 12),
        vec![0.625, 0.0, 0.375]
    );
}
