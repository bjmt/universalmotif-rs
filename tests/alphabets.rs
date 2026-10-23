//! Scanning and comparison of amino acid and custom-alphabet motifs against
//! the R package's full `scan_sequences()` and `compare_motifs()`.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Motif;
use universalmotif::compare::{CompareOptions, NullModel, compare_motifs};
use universalmotif::merge::{MergeOptions, merge_motifs};
use universalmotif::scan::{ScanOptions, scan_table};

fn reference() -> Value {
    let path = format!(
        "{}/tests/reference/alphabets.json",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

fn motifs(r: &Value, key: &str) -> Vec<Motif> {
    r[key]
        .as_array()
        .unwrap()
        .iter()
        .map(motif_from_record)
        .collect()
}

fn check_scans(motifs: &[Motif], seqs: &[String], cases: &Value, label: &str) {
    for case in cases.as_array().unwrap() {
        let threshold = num(&case["threshold"]);
        let n = case["table"]["start"].as_array().map_or(0, |a| a.len());
        let ms = &motifs[..motifs.len().min(if label == "aa" && threshold > 0.005 {
            3
        } else {
            usize::MAX
        })];
        let ctx = format!("{label} scan at {threshold}");
        let opts = ScanOptions {
            pvalue: threshold,
            rc: true,
            nthreads: 2,
            ..Default::default()
        };
        let t = scan_table(ms, seqs, &opts).unwrap();
        assert_eq!(t.len(), n, "{ctx}: number of hits");
        let w = &case["table"];
        let f = |v: &[i64]| v.iter().map(|&x| x as f64).collect::<Vec<_>>();
        let plus1 = |v: &[i64]| v.iter().map(|&x| x as f64 + 1.0).collect::<Vec<_>>();
        assert_eq!(plus1(&t.motif), nums(&w["motif.i"]), "{ctx}: motif");
        assert_eq!(
            plus1(&t.sequence),
            nums(&w["sequence.i"]),
            "{ctx}: sequence"
        );
        // R counts positions from 1.
        assert_eq!(plus1(&t.start), nums(&w["start"]), "{ctx}: start");
        assert_eq!(f(&t.end), nums(&w["stop"]), "{ctx}: end");
        assert!(t.minus.iter().all(|m| !m), "{ctx}: only the forward strand");
        assert_close_vec(&t.score, &nums(&w["score"]), &format!("{ctx}: score"));
        let matched: Vec<String> = (0..t.len()).map(|i| t.matched(i, seqs)).collect();
        assert_eq!(matched, texts(&w["match"]), "{ctx}: match");
        // The R package's P-values come from the same distribution but are
        // printed through a different path, so allow for rounding.
        for (g, e) in t.pvalue.iter().zip(nums(&w["pvalue"])) {
            assert!(
                (g - e).abs() <= 1e-9 * e.max(1e-300) + 1e-15,
                "{ctx}: pvalue {g} vs {e}"
            );
        }
    }
}

#[test]
fn protein_and_custom_scans_match_r() {
    let r = reference();
    // R's rounding of counts is reproduced here: see `r_style_pwm`.
    let aa: Vec<Motif> = motifs(&r, "aa_motifs").iter().map(r_style_pwm).collect();
    check_scans(&aa, &texts(&r["aa_sequences"]), &r["aa_scan"], "aa");
    check_scans(
        &motifs(&r, "custom_motifs"),
        &texts(&r["custom_sequences"]),
        &r["custom_scan"],
        "custom",
    );
}

fn check_compare(motifs: &[Motif], cases: &Value, subsets: &[&[usize]], label: &str) {
    for (case, subset) in cases.as_array().unwrap().iter().zip(subsets) {
        let ms: Vec<Motif> = subset.iter().map(|&i| motifs[i].clone()).collect();
        let n = ms.len();
        let opts = CompareOptions {
            min_overlap: num(&case["min_overlap"]) as usize,
            nthreads: 2,
            ..Default::default()
        };
        let rows = compare_motifs(&ms, &(0..n).collect::<Vec<_>>(), &opts).unwrap();
        let want: Vec<Vec<f64>> = case["sums"].as_array().unwrap().iter().map(nums).collect();
        for c in &rows {
            let sum = c.score * c.overlap as f64;
            assert!(
                close(sum, want[c.query][c.target]),
                "{label} {}x{}: {sum} vs {}",
                c.query,
                c.target,
                want[c.query][c.target]
            );
            assert!(
                (0.0..=1.0).contains(&c.pvalue),
                "{label}: pvalue {}",
                c.pvalue
            );
            assert_eq!(c.query_consensus.len(), c.overlap);
        }
    }
}

#[test]
fn protein_and_custom_comparisons_match_r() {
    let r = reference();
    let aa = motifs(&r, "aa_motifs");
    check_compare(
        &aa,
        &r["aa_compare"],
        &[&[0, 1, 2, 3, 4], &[0, 2, 3, 4]],
        "aa",
    );
    let cu = motifs(&r, "custom_motifs");
    check_compare(&cu, &r["custom_compare"], &[&[0, 1, 2]], "custom");

    // A motif aligned to itself is a perfect match on the forward strand,
    // and its own comparison is the most significant.
    let rows = compare_motifs(&aa, &[0], &CompareOptions::default()).unwrap();
    assert!(
        rows.iter()
            .all(|c| c.strand == universalmotif::Strand::Plus)
    );
    let best = rows
        .iter()
        .min_by(|a, b| a.pvalue.total_cmp(&b.pvalue))
        .unwrap();
    assert_eq!(best.target, 0);
    assert_eq!(rows[4].query_consensus, "MKVLAW");
    // Equal K and R makes an ambiguous consensus letter.
    assert_eq!(rows[4].target_consensus, "MXVLAW");

    // The parametric null is only for DNA and RNA.
    let opts = CompareOptions {
        null: NullModel::Parametric,
        ..Default::default()
    };
    assert!(compare_motifs(&aa, &[0], &opts).is_err());
}

#[test]
fn protein_motifs_merge() {
    let r = reference();
    let aa = motifs(&r, "aa_motifs");
    let merged = merge_motifs(
        &[aa[0].clone(), aa[4].clone()],
        &MergeOptions {
            min_overlap: 4,
            ..Default::default()
        },
    )
    .unwrap();
    // The short motif lies inside the long one, so the frame is the long
    // motif's and every column is a probability column.
    assert_eq!(merged.width(), aa[0].width());
    assert_eq!(merged.alphabet, universalmotif::Alphabet::Aa);
    for c in merged.matrix.columns() {
        assert!((c.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }
    assert!(merged.consensus.as_deref().unwrap().starts_with("MKVLAW"));
}
