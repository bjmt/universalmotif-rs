//! Sequence scanning against `scan_sequences_lite()`.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Strand;
use universalmotif::io::fasta;
use universalmotif::scan::{NoOverlaps, OverlapPriority, ScanOptions, scan_sequences};

/// R's make.unique(): later duplicates get ".1", ".2" and so on.
fn make_unique(names: &[String]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for n in names {
        let mut name = n.clone();
        let mut k = 0;
        while out.contains(&name) {
            k += 1;
            name = format!("{n}.{k}");
        }
        out.push(name);
    }
    out
}

fn ints(v: &Value) -> Vec<usize> {
    nums(v).iter().map(|&x| x as usize).collect()
}

fn texts(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap().to_string())
        .collect()
}

#[test]
fn scans_match_r() {
    let motifs: Vec<_> = load_key("scan.json", "motifs")
        .iter()
        .map(motif_from_record)
        .collect();
    let names = make_unique(&motifs.iter().map(|m| m.name.clone()).collect::<Vec<_>>());

    for (i, case) in load("scan.json").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("scan case {i} ({} {:?})", case["file"], a);
        let records = fasta::read(&reference_file(case["file"].as_str().unwrap())).unwrap();
        let records: Vec<_> = if case["file"] == "extra.fa" {
            records.into_iter().filter(|(n, _)| n != "short").collect()
        } else {
            records
        };
        let seqs: Vec<&str> = records.iter().map(|(_, s)| s.as_str()).collect();
        let b = |k: &str| a.get(k).and_then(Value::as_bool).unwrap_or(false);
        let opts = ScanOptions {
            pvalue: a.get("pvalue").map(num).unwrap_or(1e-4),
            rc: a.get("RC").and_then(Value::as_bool).unwrap_or(true),
            no_overlaps: b("no.overlaps").then(|| NoOverlaps {
                by: match a.get("no.overlaps.by").and_then(Value::as_str) {
                    Some("score") => OverlapPriority::Score,
                    _ => OverlapPriority::Pvalue,
                },
                ignore_strand: b("no.overlaps.by.strand"),
                ignore_motif: b("no.overlaps.by.motif"),
            }),
            nthreads: 1,
            ..Default::default()
        };
        let hits = scan_sequences(&motifs, &seqs, &opts).unwrap();
        let w = &case["hits"];
        assert_eq!(hits.len(), ints(&w["start"]).len(), "{ctx}: number of hits");

        let col = |f: &dyn Fn(&universalmotif::scan::Hit) -> String| {
            hits.iter().map(f).collect::<Vec<_>>()
        };
        assert_eq!(
            col(&|h| names[h.motif].clone()),
            texts(&w["motif"]),
            "{ctx}: motif"
        );
        assert_eq!(
            hits.iter().map(|h| h.motif + 1).collect::<Vec<_>>(),
            ints(&w["motif_i"]),
            "{ctx}: motif index"
        );
        assert_eq!(
            col(&|h| records[h.sequence].0.clone()),
            texts(&w["sequence"]),
            "{ctx}: sequence"
        );
        assert_eq!(
            hits.iter().map(|h| h.start).collect::<Vec<_>>(),
            ints(&w["start"]),
            "{ctx}: start"
        );
        assert_eq!(
            hits.iter().map(|h| h.end).collect::<Vec<_>>(),
            ints(&w["end"]),
            "{ctx}: end"
        );
        assert_eq!(
            col(&|h| if h.strand == Strand::Minus {
                "-".into()
            } else {
                "+".into()
            }),
            texts(&w["strand"]),
            "{ctx}: strand"
        );
        assert_eq!(
            col(&|h| h.matched.clone()),
            texts(&w["match"]),
            "{ctx}: match"
        );
        let f =
            |g: &dyn Fn(&universalmotif::scan::Hit) -> f64| hits.iter().map(g).collect::<Vec<_>>();
        assert_close_vec(
            &f(&|h| h.score),
            &nums(&w["score"]),
            &format!("{ctx}: score"),
        );
        assert_close_vec(
            &f(&|h| h.score_pct),
            &nums(&w["score_pct"]),
            &format!("{ctx}: score_pct"),
        );
        assert_close_vec(
            &f(&|h| h.pvalue),
            &nums(&w["pvalue"]),
            &format!("{ctx}: pvalue"),
        );
    }
}

#[test]
fn short_sequences_and_lower_case() {
    let m = universalmotif::create::from_strings(
        &["TTGACT"],
        None,
        &universalmotif::create::CreateOptions::default(),
    )
    .unwrap();
    let opts = ScanOptions {
        pvalue: 0.01,
        ..Default::default()
    };
    let hits = scan_sequences(std::slice::from_ref(&m), &["ACG", "aattgactaa"], &opts).unwrap();
    assert!(hits.iter().all(|h| h.sequence == 1));
    let fwd = hits.iter().find(|h| h.strand == Strand::Plus).unwrap();
    assert_eq!((fwd.start, fwd.end, fwd.matched.as_str()), (3, 8, "TTGACT"));
}
