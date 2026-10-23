//! Scoring helpers, column metrics, gaps, alphabets, masking and small
//! helpers against the R package.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Motif;
use universalmotif::io::meme::{AlphabetDefinition, write_alphabet};
use universalmotif::mask::{mask_ranges, mask_seqs};
use universalmotif::metrics::{Average, Column, ColumnMetric, average_ic, compare_columns};
use universalmotif::ops::{add_gap, switch_alph, ungap};
use universalmotif::pvalue::log_pvalue_text;
use universalmotif::scoring::{
    ThresholdType, get_matches, get_scores, motif_range, motif_score, prob_match, prob_match_bkg,
    score_match,
};

fn reference() -> Value {
    let path = format!("{}/tests/reference/utils.json", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

fn usizes(v: &Value) -> Vec<usize> {
    nums(v).into_iter().map(|x| x as usize).collect()
}

fn motif(r: &Value, name: &str) -> Motif {
    motif_from_record(&r["motifs"][name])
}

#[test]
fn scoring_helpers_match_r() {
    let r = reference();
    for case in r["scoring"].as_array().unwrap() {
        let f = case["fn"].as_str().unwrap();
        let m = motif(&r, case["motif"].as_str().unwrap());
        let a = &case["args"];
        let af = a["allow_nonfinite"].as_bool().unwrap_or(false);
        let ctx = format!("{f} {} {a}", case["motif"]);
        if case.get("error").is_some() {
            let got = match f {
                "get_matches" => get_matches(&m, num(&a["score"]), af).map(|_| ()),
                _ => panic!("{ctx}: unexpected error case"),
            };
            assert!(got.is_err(), "{ctx}: expected an error");
            continue;
        }
        let want = &case["value"];
        match f {
            "motif_range" => {
                let (lo, hi) = motif_range(&m, af).unwrap();
                assert_close_vec(&[lo, hi], &nums(want), &ctx);
            }
            "motif_score" => {
                let t = a["threshold_type"]
                    .as_str()
                    .unwrap()
                    .parse::<ThresholdType>()
                    .unwrap();
                let got = motif_score(&m, &nums(&a["threshold"]), t, af).unwrap();
                assert_close_vec(&got, &nums(want), &ctx);
            }
            "score_match" => {
                // R adds the rounded entries as floats and truncates the total
                // again, so it can be 0.001 off the integer grid sum that
                // scanning uses and the port returns.
                let got = score_match(&m, &texts(&a["matches"]), af).unwrap();
                let want = nums(want);
                assert_eq!(got.len(), want.len(), "{ctx}");
                for (g, w) in got.iter().zip(&want) {
                    assert!(
                        g == w || (g - w).abs() <= 0.001 + 1e-9,
                        "{ctx}: got {g}, want {w}"
                    );
                }
            }
            "get_scores" => {
                assert_close_vec(&get_scores(&m, af).unwrap(), &nums(want), &ctx);
            }
            "get_matches" => {
                let got = get_matches(&m, num(&a["score"]), af).unwrap();
                assert_eq!(got, texts(want), "{ctx}");
            }
            "prob_match" => {
                let az = a["allow_zero"].as_bool().unwrap();
                let got = prob_match(&m, &texts(&a["matches"]), az).unwrap();
                assert_close_vec(&got, &nums(want), &ctx);
            }
            _ => panic!("unknown function {f}"),
        }
    }
    let p = &r["prob_match_bkg"];
    let bkg: Vec<(char, f64)> = texts(&p["letters"])
        .iter()
        .map(|l| l.chars().next().unwrap())
        .zip(nums(&p["bkg"]))
        .collect();
    let got = prob_match_bkg(&bkg, &texts(&p["matches"])).unwrap();
    assert_close_vec(&got, &nums(&p["value"]), "prob_match_bkg");
    assert!(prob_match_bkg(&bkg, &["ACGN"]).is_err());

    // The R package gives NA here; the lowest score is minus infinity.
    let zero = motif(&r, "zero");
    let (lo, hi) = motif_range(&zero, true).unwrap();
    assert_eq!(lo, f64::NEG_INFINITY);
    assert_eq!(hi, motif_range(&zero, false).unwrap().1.max(hi));
    assert!(score_match(&zero, &["TATAT"], false).is_err());
    assert!(score_match(&zero, &["TATANAT"], false).is_err());
    // Lower case letters match too.
    assert_eq!(
        score_match(&zero, &["tatatat"], true).unwrap(),
        score_match(&zero, &["TATATAT"], true).unwrap()
    );
    assert!(get_matches(&zero, f64::INFINITY, true).unwrap().is_empty());
}

#[test]
fn column_metrics_and_averages_match_r() {
    let r = reference();
    for case in r["compare_columns"].as_array().unwrap() {
        let metric = case["method"]
            .as_str()
            .unwrap()
            .parse::<ColumnMetric>()
            .unwrap();
        let (x, y) = (nums(&case["x"]), nums(&case["y"]));
        let (bkg1, bkg2) = (nums(&case["bkg1"]), nums(&case["bkg2"]));
        let got = compare_columns(
            &Column::new(&x)
                .with_bkg(&bkg1)
                .with_nsites(num(&case["nsites1"])),
            &Column::new(&y)
                .with_bkg(&bkg2)
                .with_nsites(num(&case["nsites2"])),
            metric,
        );
        let ctx = format!("{} {} {}", case["method"], case["x"], case["y"]);
        match case.get("value") {
            // The R package does not normalize the WPCC weights, so the port
            // is checked against a weighted correlation computed here instead.
            Some(_) if metric == ColumnMetric::Wpcc => {
                let want = weighted_pcc(
                    &nums(&case["x"]),
                    &nums(&case["y"]),
                    &nums(&case["bkg1"]),
                    &nums(&case["bkg2"]),
                );
                assert!(close(got.unwrap(), want), "{ctx}");
            }
            Some(v) => assert!(close(got.unwrap(), num(v)), "{ctx}"),
            None => assert!(got.is_err(), "{ctx}"),
        }
    }
    let x = [0.5, 0.1, 0.1, 0.3];
    let y = [0.7, 0.1, 0.1, 0.1];
    let u = [0.25; 4];
    let (cx, cy) = (Column::new(&x).with_bkg(&u), Column::new(&y).with_bkg(&u));
    let pcc = compare_columns(&cx, &cy, ColumnMetric::Pcc).unwrap();
    let wpcc = compare_columns(&cx, &cy, ColumnMetric::Wpcc).unwrap();
    assert!(close(pcc, wpcc), "WPCC with a uniform background is PCC");
    let half = [0.5, 0.5, 0.0, 0.0];
    let kl = compare_columns(&Column::new(&half), &Column::new(&half), ColumnMetric::Kl).unwrap();
    assert_eq!(kl, 0.0, "zeros in both columns add nothing to KL");
    assert!(
        compare_columns(
            &Column::new(&[0.5, 0.5]),
            &Column::new(&[1.0]),
            ColumnMetric::Pcc
        )
        .is_err()
    );
    let avg = &r["average_ic"];
    let motifs: Vec<Motif> = avg["motifs"]
        .as_array()
        .unwrap()
        .iter()
        .map(motif_from_record)
        .collect();
    for case in avg["cases"].as_array().unwrap() {
        let how = case["average"]
            .as_str()
            .unwrap()
            .parse::<Average>()
            .unwrap();
        let got = average_ic(&motifs, how).unwrap();
        assert_close_vec(&got, &nums(&case["value"]), &format!("average_ic {how:?}"));
    }
}

#[test]
fn gaps_and_alphabets_match_r() {
    let r = reference();
    let random = motif(&r, "random");
    for case in r["add_gap"].as_array().unwrap() {
        let a = args(&case["args"]);
        let gaploc = a.get("gaploc").map(usizes);
        let mingap = a.get("mingap").map_or(vec![1], usizes);
        let maxgap = a.get("maxgap").map_or(vec![5], usizes);
        let got = add_gap(&random, gaploc.as_deref(), &mingap, &maxgap);
        let ctx = format!("add_gap {a:?}");
        match case.get("motif") {
            Some(want) => assert_motif_eq(&got.unwrap(), want, &ctx),
            None => assert!(got.is_err(), "{ctx}"),
        }
    }
    let gapped = add_gap(&random, Some(&[2]), &[1], &[4]).unwrap();
    for case in r["ungap"].as_array().unwrap() {
        let delete = case["delete"].as_bool().unwrap();
        assert_motif_eq(&ungap(&gapped, delete), &case["motif"], "ungap");
    }
    for case in r["switch_alph"].as_array().unwrap() {
        let input = motif_from_record(&case["input"]);
        match case.get("output") {
            Some(want) => assert_motif_eq(&switch_alph(&input).unwrap(), want, "switch_alph"),
            None => assert!(switch_alph(&input).is_err()),
        }
    }
}

#[test]
fn masking_matches_r() {
    let r = reference();
    for case in r["mask_seqs"].as_array().unwrap() {
        let alphabet = alphabet(case["alphabet"].as_str().unwrap());
        let letter = case["letter"].as_str().unwrap().chars().next().unwrap();
        let got = mask_seqs(
            &texts(&case["seqs"]),
            case["pattern"].as_str().unwrap(),
            case["rc"].as_bool().unwrap(),
            letter,
            &alphabet,
        );
        let ctx = format!("mask_seqs {case}");
        match case.get("value") {
            Some(v) => assert_eq!(got.unwrap(), texts(v), "{ctx}"),
            None => assert!(got.is_err(), "{ctx}"),
        }
    }
    let mr = &r["mask_ranges"];
    let (seqs, names) = (texts(&mr["seqs"]), texts(&mr["names"]));
    // The R package fails on the first and third cases (when every sequence
    // has exactly one range, or one has none), so those
    // were worked out by hand.
    let fixed: [(usize, [&str; 2]); 2] = [
        (0, ["------GTACGTACGT", "GGGG----AAAATTTT"]),
        (2, ["NNNNNCGTACGTANNN", "GGGGCCCCAAAATTTT"]),
    ];
    for (i, case) in mr["cases"].as_array().unwrap().iter().enumerate() {
        let chroms = texts(&case["chrom"]);
        let (starts, ends) = (nums(&case["start"]), nums(&case["end"]));
        let letter = case["letter"].as_str().unwrap().chars().next().unwrap();
        let got: Result<Vec<String>, String> = (|| {
            if let Some(c) = chroms.iter().find(|c| !names.contains(c)) {
                return Err(format!("{c} is not a sequence"));
            }
            names
                .iter()
                .zip(&seqs)
                .map(|(n, s)| {
                    let ranges: Vec<(i64, i64)> = (0..chroms.len())
                        .filter(|&i| chroms[i] == *n)
                        .map(|i| (starts[i] as i64 - 1, ends[i] as i64))
                        .collect();
                    mask_ranges(s, &ranges, letter).map_err(|e| e.to_string())
                })
                .collect()
        })();
        if let Some((_, want)) = fixed.iter().find(|(k, _)| *k == i) {
            assert_eq!(got.unwrap(), want.to_vec(), "mask_ranges {case}");
            continue;
        }
        match case.get("value") {
            Some(v) => assert_eq!(got.unwrap(), texts(v), "mask_ranges {case}"),
            None => assert!(got.is_err(), "mask_ranges {case}"),
        }
    }
}

#[test]
fn small_helpers_match_r() {
    let r = reference();
    for case in r["meme_alph"].as_array().unwrap() {
        let a = &case["args"];
        let pairs = |v: &Value| -> Vec<(char, String)> {
            v.as_object().map_or(vec![], |o| {
                o.iter()
                    .map(|(k, v)| (k.chars().next().unwrap(), v.as_str().unwrap().into()))
                    .collect()
            })
        };
        let def = AlphabetDefinition {
            core: a["core"].as_str().unwrap().into(),
            complements: opt_str(&a["complements"]),
            ambiguity: a["ambiguity"].as_object().map_or(vec![], |o| {
                o.iter()
                    .map(|(k, v)| (k.clone(), v.as_str().unwrap().into()))
                    .collect()
            }),
            like: opt_str(&a["like"]),
            name: opt_str(&a["alph.name"]),
            letter_names: pairs(&a["letter.names"]),
            colours: pairs(&a["colours"]),
        };
        let got = write_alphabet(&def);
        match case.get("value") {
            Some(v) => {
                // The R package declares each complementary pair twice, once
                // from each letter; the port declares it once.
                let mut seen: Vec<(String, String)> = Vec::new();
                let want: Vec<&str> = v
                    .as_str()
                    .unwrap()
                    .lines()
                    .filter(|l| {
                        let Some((x, y)) = l.split_once(" ~ ") else {
                            return true;
                        };
                        let (x, y) = (x[..1].to_string(), y[..1].to_string());
                        let mirrored = seen.contains(&(y.clone(), x.clone()));
                        seen.push((x, y));
                        !mirrored
                    })
                    .collect();
                assert_eq!(got.unwrap().trim_end(), want.join("\n"), "{a}");
            }
            None => assert!(got.is_err(), "{a}"),
        }
    }
    for case in r["log_string_pval"].as_array().unwrap() {
        let got = log_pvalue_text(case["pval"].as_str().unwrap()).unwrap();
        assert!(close(got, num(&case["value"])), "{case}");
    }
    // Far below the smallest double, where the R package fails on arm64.
    let tiny = log_pvalue_text("2.5e-400").unwrap();
    assert!((tiny - (2.5f64.ln() - 400.0 * std::f64::consts::LN_10)).abs() < 1e-9);
    assert_eq!(log_pvalue_text("0e-5").unwrap(), f64::NEG_INFINITY);
    assert!(log_pvalue_text("-1").is_err() && log_pvalue_text("abc").is_err());
}

/// Pearson correlation of `x` and `y`, with each letter weighted by its share
/// of the two backgrounds.
fn weighted_pcc(x: &[f64], y: &[f64], bkg1: &[f64], bkg2: &[f64]) -> f64 {
    let w: Vec<f64> = bkg1.iter().zip(bkg2).map(|(a, b)| a + b).collect();
    let ws: f64 = w.iter().sum();
    let mean = |v: &[f64]| v.iter().zip(&w).map(|(a, w)| a * w).sum::<f64>() / ws;
    let (mx, my) = (mean(x), mean(y));
    let cov = |a: &[f64], b: &[f64], ma: f64, mb: f64| {
        a.iter()
            .zip(b)
            .zip(&w)
            .map(|((a, b), w)| w * (a - ma) * (b - mb))
            .sum::<f64>()
    };
    let bot = (cov(x, x, mx, mx) * cov(y, y, my, my)).sqrt();
    // Constant columns count as perfectly correlated, as with PCC.
    if bot == 0.0 {
        1.0
    } else {
        cov(x, y, mx, my) / bot
    }
}
