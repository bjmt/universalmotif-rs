//! Comparison, merging and trees against the R package's lite functions.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Motif;
use universalmotif::cluster::{Method, TreeOptions, hclust, motif_tree};
use universalmotif::compare::{CompareOptions, NullModel, compare_matrix, compare_motifs};
use universalmotif::merge::{MergeOptions, MergeSimilarOptions, merge_motifs, merge_similar};

fn set(name: &str) -> Vec<Motif> {
    load_key("compare.json", name)
        .iter()
        .map(motif_from_record)
        .collect()
}

fn compare_opts(a: &serde_json::Map<String, Value>) -> CompareOptions {
    let mut o = CompareOptions::default();
    if let Some(v) = a.get("RC").and_then(Value::as_bool) {
        o.rc = v;
    }
    if let Some(v) = a.get("min.overlap") {
        o.min_overlap = num(v) as usize;
    }
    if a.get("null").and_then(Value::as_str) == Some("parametric") {
        o.null = NullModel::Parametric;
    }
    if let Some(b) = a.get("bkg") {
        let b = nums(b);
        o.bkg = [b[0], b[1], b[2], b[3]];
    }
    o
}

fn ints(v: &Value) -> Vec<i64> {
    nums(v).iter().map(|&x| x as i64).collect()
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

#[test]
fn comparisons_match_r() {
    for (i, case) in load_key("compare.json", "compare").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("compare case {i} ({} {a:?})", case["set"]);
        let motifs = set(case["set"].as_str().unwrap());
        let opts = compare_opts(&a);
        if let Some(rows) = case["matrix"].as_array() {
            let out = a
                .get("matrix.out")
                .and_then(Value::as_str)
                .unwrap_or("score");
            let m = compare_matrix(&motifs, &opts, out != "score").unwrap();
            let got = match out {
                "score" => m.scores,
                "pvalue" => m.pvalues,
                _ => m.qvalues,
            };
            for (r, (g, w)) in got.iter().zip(rows).enumerate() {
                assert_close_vec(g, &nums(w), &format!("{ctx}: row {r}"));
            }
        } else {
            let w = &case["long"];
            let queries: Vec<usize> = nums(&a["compare.to"])
                .iter()
                .map(|&x| x as usize - 1)
                .collect();
            let qmax = a.get("qvalue").map(num).unwrap_or(0.1);
            let mut rows: Vec<_> = compare_motifs(&motifs, &queries, &opts)
                .unwrap()
                .into_iter()
                .filter(|c| c.qvalue <= qmax && c.overlap > 0)
                .collect();
            rows.sort_by(|x, y| {
                x.qvalue
                    .total_cmp(&y.qvalue)
                    .then(x.pvalue.total_cmp(&y.pvalue))
                    .then(x.target.cmp(&y.target))
            });
            let col = |f: &dyn Fn(&universalmotif::compare::Comparison) -> i64| {
                rows.iter().map(f).collect::<Vec<_>>()
            };
            assert_eq!(
                col(&|c| c.query as i64 + 1),
                ints(&w["subject_i"]),
                "{ctx}: subject"
            );
            assert_eq!(
                col(&|c| c.target as i64 + 1),
                ints(&w["target_i"]),
                "{ctx}: target"
            );
            assert_eq!(col(&|c| c.offset), ints(&w["offset"]), "{ctx}: offset");
            assert_eq!(
                col(&|c| c.overlap as i64),
                ints(&w["overlap"]),
                "{ctx}: overlap"
            );
            let strands: Vec<String> = rows.iter().map(|c| c.strand.as_str().to_string()).collect();
            assert_eq!(strands, texts(&w["strand"]), "{ctx}: strand");
            let qc: Vec<String> = rows.iter().map(|c| c.query_consensus.clone()).collect();
            let tc: Vec<String> = rows.iter().map(|c| c.target_consensus.clone()).collect();
            assert_eq!(
                qc,
                texts(&w["subject_consensus"]),
                "{ctx}: subject consensus"
            );
            assert_eq!(tc, texts(&w["target_consensus"]), "{ctx}: target consensus");
            let f = |g: &dyn Fn(&universalmotif::compare::Comparison) -> f64| {
                rows.iter().map(g).collect::<Vec<_>>()
            };
            assert_close_vec(
                &f(&|c| c.score),
                &nums(&w["score"]),
                &format!("{ctx}: score"),
            );
            assert_close_vec(
                &f(&|c| c.pvalue),
                &nums(&w["pvalue"]),
                &format!("{ctx}: pvalue"),
            );
            assert_close_vec(
                &f(&|c| c.qvalue),
                &nums(&w["qvalue"]),
                &format!("{ctx}: qvalue"),
            );
        }
    }
}

#[test]
fn merges_match_r() {
    for (i, case) in load_key("compare.json", "merge").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("merge case {i} ({} {:?} {a:?})", case["set"], case["idx"]);
        let all = set(case["set"].as_str().unwrap());
        let motifs: Vec<Motif> = nums(&case["idx"])
            .iter()
            .map(|&k| all[k as usize - 1].clone())
            .collect();
        let opts = MergeOptions {
            min_overlap: a.get("min.overlap").map(|v| num(v) as usize).unwrap_or(5),
            rc: a.get("RC").and_then(Value::as_bool).unwrap_or(true),
            name: a.get("new.name").and_then(Value::as_str).map(String::from),
            weighted: a.get("weighted").and_then(Value::as_bool).unwrap_or(false),
            nthreads: 1,
        };
        let got = merge_motifs(&motifs, &opts).unwrap();
        assert_motif_eq(&got, &case["motif"], &ctx);
    }
}

#[test]
fn merge_similar_matches_r() {
    let motifs = set("dna");
    for (i, case) in load_key("compare.json", "similar").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("merge_similar case {i} ({a:?})");
        let c = compare_opts(&a);
        let opts = MergeSimilarOptions {
            qvalue: a.get("qvalue").map(num).unwrap_or(0.01),
            null: c.null,
            min_overlap: c.min_overlap,
            rc: c.rc,
            bkg: c.bkg,
            weighted: a.get("weighted").and_then(Value::as_bool).unwrap_or(false),
            nthreads: 1,
        };
        let (clusters, merged) = merge_similar(&motifs, &opts).unwrap();
        let want: Vec<usize> = nums(&case["clusters"])
            .iter()
            .map(|&x| x as usize - 1)
            .collect();
        assert_eq!(clusters, want, "{ctx}: clusters");
        let want = case["motifs"].as_array().unwrap();
        assert_eq!(merged.len(), want.len(), "{ctx}: number of motifs");
        for (k, (g, w)) in merged.iter().zip(want).enumerate() {
            assert_motif_eq(g, w, &format!("{ctx}, motif {k}"));
        }
    }
}

#[test]
fn trees_match_r() {
    let motifs = set("dna");
    for (i, case) in load_key("compare.json", "tree").iter().enumerate() {
        let method = Method::parse(case["method"].as_str().unwrap()).unwrap();
        let ctx = format!("tree case {i} ({} {:?})", case["source"], method);
        let dist: Vec<Vec<f64>> = case["dist"].as_array().unwrap().iter().map(nums).collect();
        let tree = hclust(&dist, method).unwrap();
        let merge: Vec<Vec<i64>> = case["merge"].as_array().unwrap().iter().map(ints).collect();
        assert_eq!(
            tree.merge.iter().map(|m| m.to_vec()).collect::<Vec<_>>(),
            merge,
            "{ctx}: merge"
        );
        assert_close_vec(
            &tree.height,
            &nums(&case["height"]),
            &format!("{ctx}: height"),
        );
        let order: Vec<i64> = tree.order.iter().map(|&o| o as i64 + 1).collect();
        assert_eq!(order, ints(&case["order"]), "{ctx}: order");
        let labels = texts(&case["labels"]);
        assert_eq!(
            tree.newick(&labels),
            case["newick"].as_str().unwrap(),
            "{ctx}: newick"
        );

        if case["source"] == "motifs" {
            let a = args(&case["args"]);
            let opts = TreeOptions {
                min_overlap: num(&a["min.overlap"]) as usize,
                rc: a["RC"].as_bool().unwrap(),
                method,
                nthreads: 1,
            };
            // Heights can differ from R in the last bits, because the
            // correlation sums round slightly differently.
            let own = motif_tree(&motifs, &opts).unwrap();
            assert_eq!(own.merge, tree.merge, "{ctx}: motif_tree merge");
            assert_eq!(own.order, tree.order, "{ctx}: motif_tree order");
            assert_close_vec(
                &own.height,
                &tree.height,
                &format!("{ctx}: motif_tree height"),
            );
        }
    }
}

#[test]
fn linkage_matrix() {
    let dist = vec![
        vec![0.0, 1.0, 4.0, 5.0],
        vec![1.0, 0.0, 3.0, 4.5],
        vec![4.0, 3.0, 0.0, 2.0],
        vec![5.0, 4.5, 2.0, 0.0],
    ];
    let t = hclust(&dist, Method::Complete).unwrap();
    assert_eq!(
        t.linkage(),
        vec![
            [0.0, 1.0, 1.0, 2.0],
            [2.0, 3.0, 2.0, 2.0],
            [4.0, 5.0, 5.0, 4.0]
        ]
    );
}

#[test]
fn strand_ties_go_to_the_forward_strand() {
    use universalmotif::Strand;
    use universalmotif::compare::align;
    use universalmotif::create::{CreateOptions, from_strings};
    // TTGACGTCAA is its own reverse complement, so both strands of the
    // query score the same against it.
    let opts = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let q = from_strings(&["GGGGCCCCAT"], None, &opts).unwrap();
    let t = from_strings(&["TTGACGTCAA"], None, &opts).unwrap();
    let al = align(&[q, t], &[(0, 1), (1, 0)], &CompareOptions::default()).unwrap();
    assert!(al.iter().all(|a| a.strand == Strand::Plus), "{al:?}");
}
