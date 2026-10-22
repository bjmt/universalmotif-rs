//! Gapped scanning against the R package's `scan_sequences()`, which only
//! manages the forward strand (see DIVERGENCES.md); both strands are checked
//! against plain scoring in the crate's unit tests.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Motif;
use universalmotif::scan::{NoOverlaps, ScanOptions, scan_table};

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

#[test]
fn gapped_scans_match_r() {
    let path = format!("{}/tests/reference/gapped.json", env!("CARGO_MANIFEST_DIR"));
    let r: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let motifs: Vec<Motif> = r["motifs"]
        .as_array()
        .unwrap()
        .iter()
        .map(motif_from_record)
        .collect();
    assert!(motifs[0].gap.isgapped && motifs[1].gap.isgapped);
    let seqs = texts(&r["sequences"]);
    for case in r["cases"].as_array().unwrap() {
        let opts = ScanOptions {
            pvalue: num(&case["threshold"]),
            rc: false,
            no_overlaps: case["no_overlaps"]
                .as_bool()
                .unwrap()
                .then(NoOverlaps::default),
            nthreads: 2,
            ..Default::default()
        };
        let ctx = format!(
            "gapped scan at {} {}",
            case["threshold"], case["no_overlaps"]
        );
        let t = scan_table(&motifs, &seqs, &opts).unwrap();
        let mut got: Vec<(i64, i64, i64, i64, String, f64, f64)> = (0..t.len())
            .map(|i| {
                (
                    t.motif[i] + 1,
                    t.sequence[i] + 1,
                    t.start[i],
                    t.end[i],
                    t.matched(i, &seqs),
                    t.score[i],
                    t.pvalue[i],
                )
            })
            .collect();
        got.sort_by(|a, b| (a.0, a.1, a.2, &a.4).cmp(&(b.0, b.1, b.2, &b.4)));
        let w = &case["table"];
        let want_match = texts(&w["match"]);
        assert_eq!(got.len(), want_match.len(), "{ctx}: number of hits");
        let col = |k: &str| nums(&w[k]);
        let (m, s, st, sp, sc, pv) = (
            col("motif.i"),
            col("sequence.i"),
            col("start"),
            col("stop"),
            col("score"),
            col("pvalue"),
        );
        for (i, g) in got.iter().enumerate() {
            let want = (m[i] as i64, s[i] as i64, st[i] as i64, sp[i] as i64);
            assert_eq!((g.0, g.1, g.2, g.3), want, "{ctx}: row {i}");
            assert_eq!(g.4, want_match[i], "{ctx}: row {i}");
            assert!(close(g.5, sc[i]), "{ctx}: score {} vs {}", g.5, sc[i]);
            assert!(
                (g.6 - pv[i]).abs() <= 1e-9 * pv[i],
                "{ctx}: pvalue {} vs {}",
                g.6,
                pv[i]
            );
        }
    }
}
