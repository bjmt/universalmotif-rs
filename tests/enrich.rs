//! Statistics, cutree, enrichment, co-occurrence and positional enrichment
//! against the R package.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::Motif;
use universalmotif::cluster::{Linkage, hclust};
use universalmotif::coocc::{CooccHit, CooccOptions, motif_coocc};
use universalmotif::enrich::{EnrichOptions, EnrichTest, enrich_motifs};
use universalmotif::peaks::{PeakHit, PeaksMode, PeaksOptions, motif_peaks};
use universalmotif::scan::{NoOverlaps, OverlapPriority, ScanOptions, scan_sequences};
use universalmotif::stats;

const FILE: &str = "enrich.json";

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

fn reference() -> Value {
    let path = format!("{}/tests/reference/{FILE}", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn opt(x: Option<f64>) -> f64 {
    x.unwrap_or(f64::NAN)
}

fn motifs(r: &Value) -> Vec<Motif> {
    r["motifs"]
        .as_array()
        .unwrap()
        .iter()
        .map(motif_from_record)
        .collect()
}

/// R's make.unique().
fn make_unique(names: &[String]) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    names
        .iter()
        .map(|n| {
            let (mut new, mut k) = (n.clone(), 0);
            while seen.contains(&new) {
                k += 1;
                new = format!("{n}.{k}");
            }
            seen.insert(new.clone());
            new
        })
        .collect()
}

#[test]
fn distributions_match_r() {
    let r = reference();
    let s = &r["stats"];
    for c in s["phyper"].as_array().unwrap() {
        let (x, m, n, k) = (num(&c["x"]), num(&c["m"]), num(&c["n"]), num(&c["k"]));
        let (xi, m, n, k) = (x as i64, m as u64, n as u64, k as u64);
        let ctx = format!("phyper({x}, {m}, {n}, {k})");
        assert!(
            close(stats::phyper(xi, m, n, k, true), num(&c["lower"])),
            "{ctx}: lower"
        );
        assert!(
            close(stats::phyper(xi, m, n, k, false), num(&c["upper"])),
            "{ctx}: upper"
        );
        if xi >= 0 {
            assert!(
                close(stats::dhyper(xi as u64, m, n, k), num(&c["density"])),
                "{ctx}: d"
            );
        }
    }
    for c in s["pbinom"].as_array().unwrap() {
        let (k, n, p) = (num(&c["k"]) as u64, num(&c["n"]) as u64, num(&c["p"]));
        let ctx = format!("pbinom({k}, {n}, {p})");
        let got = stats::pbinom_upper(k, n, p);
        assert!(
            close(got, num(&c["upper"])),
            "{ctx}: {got} vs {}",
            c["upper"]
        );
        assert!(
            close(stats::dbinom(k, n, p), num(&c["density"])),
            "{ctx}: density"
        );
    }
    for c in s["fisher"].as_array().unwrap() {
        let t: Vec<u64> = nums(&c["table"]).iter().map(|&x| x as u64).collect();
        let got = stats::fisher_greater(t[0], t[1], t[2], t[3]);
        assert!(close(got, num(&c["pvalue"])), "fisher {t:?}: {got}");
    }
    let b = &s["bonferroni"];
    assert_close_vec(
        &stats::p_adjust_bonferroni(&nums(&b["p"])),
        &nums(&b["adjusted"]),
        "bonferroni",
    );
}

#[test]
fn cutree_matches_r() {
    let r = reference();
    for c in r["cutree"].as_array().unwrap() {
        let method = c["method"].as_str().unwrap().parse::<Linkage>().unwrap();
        let dist: Vec<Vec<f64>> = c["dist"].as_array().unwrap().iter().map(nums).collect();
        let tree = hclust(&dist, method).unwrap();
        let ids = |v: &Vec<usize>| v.iter().map(|&x| x as f64 + 1.0).collect::<Vec<_>>();
        for (k, want) in c["k"].as_array().unwrap().iter().enumerate() {
            let got = tree.cutree(k + 1).unwrap();
            assert_eq!(ids(&got), nums(want), "{method:?} k={}", k + 1);
        }
        for (h, want) in nums(&c["h"]).iter().zip(c["by_h"].as_array().unwrap()) {
            let got = tree.cutree_height(*h).unwrap();
            assert_eq!(ids(&got), nums(want), "{method:?} h={h}");
        }
        assert!(tree.cutree(0).is_err() && tree.cutree(13).is_err());
    }
}

#[test]
fn enrichment_matches_r() {
    let r = reference();
    let ms = motifs(&r);
    let names = make_unique(&ms.iter().map(|m| m.name.clone()).collect::<Vec<_>>());
    let tgt = texts(&r["targets"]);
    let bkg = texts(&r["background"]);
    for case in r["enrich"].as_array().unwrap() {
        let a = args(&case["args"]);
        let ctx = format!("enrich {a:?}");
        let opts = EnrichOptions {
            pvalue: a.get("pvalue").map(num).unwrap_or(1e-4),
            test: a
                .get("test")
                .and_then(Value::as_str)
                .unwrap_or("seqs")
                .parse::<EnrichTest>()
                .unwrap(),
            rc: a.get("RC").and_then(Value::as_bool).unwrap_or(true),
            pseudocount: a.get("pseudocount").map(num).unwrap_or(1.0) as u64,
            ..Default::default()
        };
        let qmax = num(&a["qvalue"]);
        let mut got: Vec<_> = enrich_motifs(&ms, &tgt, Some(&bkg), &opts)
            .unwrap()
            .into_iter()
            .filter(|e| e.qvalue <= qmax)
            .collect();
        got.sort_by(|x, y| x.pvalue.total_cmp(&y.pvalue).then(x.motif.cmp(&y.motif)));
        let t = &case["table"];
        let gn: Vec<String> = got.iter().map(|e| names[e.motif].clone()).collect();
        assert_eq!(gn, texts(&t["motif"]), "{ctx}: motifs");
        let ints = |f: fn(&universalmotif::enrich::Enrichment) -> usize| {
            got.iter().map(|e| f(e) as f64).collect::<Vec<_>>()
        };
        assert_eq!(ints(|e| e.target_seq_n), nums(&t["target.seq.n"]), "{ctx}");
        assert_eq!(
            ints(|e| e.target_seq_hits),
            nums(&t["target.seq.hits"]),
            "{ctx}"
        );
        assert_eq!(
            ints(|e| e.target_site_hits),
            nums(&t["target.site.hits"]),
            "{ctx}"
        );
        assert_eq!(ints(|e| e.bkg_seq_n), nums(&t["bkg.seq.n"]), "{ctx}");
        assert_eq!(ints(|e| e.bkg_seq_hits), nums(&t["bkg.seq.hits"]), "{ctx}");
        assert_eq!(
            ints(|e| e.bkg_site_hits),
            nums(&t["bkg.site.hits"]),
            "{ctx}"
        );
        let f = |g: fn(&universalmotif::enrich::Enrichment) -> f64| {
            got.iter().map(g).collect::<Vec<_>>()
        };
        assert_close_vec(
            &f(|e| e.enrichment),
            &nums(&t["enrichment"]),
            &format!("{ctx}: enr"),
        );
        assert_close_vec(
            &f(|e| e.log2_enrichment),
            &nums(&t["log2.enrichment"]),
            &format!("{ctx}: log2"),
        );
        assert_close_vec(
            &f(|e| e.pvalue),
            &nums(&t["pvalue"]),
            &format!("{ctx}: pvalue"),
        );
        assert_close_vec(
            &f(|e| e.qvalue),
            &nums(&t["qvalue"]),
            &format!("{ctx}: qvalue"),
        );
    }
}

fn check_coocc(got: &[universalmotif::coocc::CooccPair], t: &Value, names: &[String], ctx: &str) {
    // The port does not test a motif against itself, so
    // self pairs are checked on their own and left out of the comparison with
    // R, whose q-values are recomputed without them.
    for p in got.iter().filter(|p| p.a == p.b) {
        assert!(
            p.odds_ratio.is_none() && p.pvalue.is_none() && p.qvalue.is_none(),
            "{ctx}"
        );
        assert_eq!(p.b_only, 0, "{ctx}: self pair b_only");
    }
    let got: Vec<_> = got.iter().filter(|p| p.a != p.b).collect();
    let ra = texts(&t["motif_a"]);
    let rb = texts(&t["motif_b"]);
    let keep: Vec<usize> = (0..ra.len()).filter(|&k| ra[k] != rb[k]).collect();
    let sel = |v: Vec<f64>| keep.iter().map(|&k| v[k]).collect::<Vec<f64>>();
    let sel_t = |v: Vec<String>| keep.iter().map(|&k| v[k].clone()).collect::<Vec<String>>();
    let na: Vec<String> = got.iter().map(|p| names[p.a].clone()).collect();
    let nb: Vec<String> = got.iter().map(|p| names[p.b].clone()).collect();
    assert_eq!(na, sel_t(ra.clone()), "{ctx}: motif_a");
    assert_eq!(nb, sel_t(rb.clone()), "{ctx}: motif_b");
    let c = |f: fn(&universalmotif::coocc::CooccPair) -> u64| {
        got.iter().map(|p| f(p) as f64).collect::<Vec<_>>()
    };
    assert_eq!(c(|p| p.a_only), sel(nums(&t["a_only"])), "{ctx}: a_only");
    assert_eq!(c(|p| p.b_only), sel(nums(&t["b_only"])), "{ctx}: b_only");
    assert_eq!(c(|p| p.both), sel(nums(&t["both"])), "{ctx}: both");
    assert_eq!(c(|p| p.neither), sel(nums(&t["neither"])), "{ctx}: neither");
    let o = |f: fn(&universalmotif::coocc::CooccPair) -> Option<f64>| {
        got.iter().map(|p| opt(f(p))).collect::<Vec<_>>()
    };
    assert_close_vec(
        &o(|p| p.odds_ratio),
        &sel(nums(&t["odds_ratio"])),
        &format!("{ctx}: or"),
    );
    let rp = sel(nums(&t["pvalue"]));
    assert_close_vec(&o(|p| p.pvalue), &rp, &format!("{ctx}: pvalue"));
    let tested: Vec<usize> = (0..rp.len()).filter(|&k| !rp[k].is_nan()).collect();
    let q = stats::p_adjust_bh(&tested.iter().map(|&k| rp[k]).collect::<Vec<_>>());
    let mut rq = vec![f64::NAN; rp.len()];
    for (&k, q) in tested.iter().zip(q) {
        rq[k] = q;
    }
    assert_close_vec(&o(|p| p.qvalue), &rq, &format!("{ctx}: qvalue"));
    if !t["both.clustered"].is_null() {
        let bc: Vec<f64> = got
            .iter()
            .map(|p| p.both_clustered.unwrap() as f64)
            .collect();
        assert_eq!(bc, sel(nums(&t["both.clustered"])), "{ctx}: both.clustered");
        assert_close_vec(
            &o(|p| p.median_distance),
            &sel(nums(&t["median.distance"])),
            &format!("{ctx}: median.distance"),
        );
    } else {
        assert!(got.iter().all(|p| p.both_clustered.is_none()));
    }
}

fn coocc_opts(a: &serde_json::Map<String, Value>) -> CooccOptions {
    CooccOptions {
        max_distance: a.get("max.distance").map(|v| num(v) as u64),
        min_coocc: a.get("min.coocc").map(num).unwrap_or(1.0) as u64,
        pseudocount: a.get("pseudocount").map(num).unwrap_or(0.0) as u64,
        self_pairs: a
            .get("self.pairs")
            .and_then(Value::as_bool)
            .unwrap_or(false),
        nthreads: 2,
    }
}

#[test]
fn coocc_matches_r() {
    let r = reference();
    let ms = motifs(&r);
    let names = make_unique(&ms.iter().map(|m| m.name.clone()).collect::<Vec<_>>());
    let tgt = texts(&r["targets"]);
    for case in r["coocc"].as_array().unwrap() {
        let a = args(&case["args"]);
        let ctx = format!("coocc {a:?}");
        let scan = ScanOptions {
            pvalue: num(&a["pvalue"]),
            rc: a.get("RC").and_then(Value::as_bool).unwrap_or(true),
            no_overlaps: (a.get("no.overlaps").and_then(Value::as_bool) == Some(true)).then_some(
                NoOverlaps {
                    by: OverlapPriority::Pvalue,
                    ignore_strand: false,
                    ignore_motif: false,
                },
            ),
            nthreads: 1,
            ..Default::default()
        };
        let hits: Vec<CooccHit> = scan_sequences(&ms, &tgt, &scan)
            .unwrap()
            .iter()
            .map(|h| CooccHit {
                motif: h.motif,
                sequence: h.sequence,
                start: h.start as i64,
            })
            .collect();
        let got = motif_coocc(&hits, ms.len(), tgt.len(), &coocc_opts(&a)).unwrap();
        check_coocc(&got, &case["table"], &names, &ctx);
    }

    let h = &r["coocc_hits"];
    let hits: Vec<CooccHit> = (0..nums(&h["motif.i"]).len())
        .map(|i| CooccHit {
            motif: nums(&h["motif.i"])[i] as usize - 1,
            sequence: nums(&h["sequence.i"])[i] as usize - 1,
            start: nums(&h["start"])[i] as i64,
        })
        .collect();
    for case in r["coocc_hit_cases"].as_array().unwrap() {
        let a = args(&case["args"]);
        let got = motif_coocc(&hits, 4, 10, &coocc_opts(&a)).unwrap();
        check_coocc(
            &got,
            &case["table"],
            &names[..4],
            &format!("coocc hits {a:?}"),
        );
    }
}

#[test]
fn peaks_match_r() {
    let r = reference();
    let h = &r["peak_hits"];
    let motif_names = texts(&h["motif"]);
    let mut sorted = motif_names.clone();
    sorted.sort();
    sorted.dedup();
    let seq_names = texts(&h["sequence"]);
    let (start, end, score) = (nums(&h["start"]), nums(&h["end"]), nums(&h["score"]));
    for case in r["peaks"].as_array().unwrap() {
        let a = args(&case["args"]);
        let ctx = format!("peaks {:?}", a.keys().collect::<Vec<_>>());
        let strand = a.get("seq.strand").and_then(Value::as_object);
        let hits: Vec<PeakHit> = (0..motif_names.len())
            .map(|i| {
                // R counts positions from 1.
                let mut centre = (start[i] - 1.0 + end[i] - 1.0) / 2.0;
                if strand.is_some_and(|s| s[&seq_names[i]] == "-") {
                    centre = 300.0 - 1.0 - centre;
                }
                PeakHit {
                    motif: sorted.iter().position(|m| *m == motif_names[i]).unwrap(),
                    sequence: seq_names[i][1..].parse().unwrap(),
                    centre,
                    score: score[i],
                }
            })
            .collect();
        let opts = PeaksOptions {
            mode: a
                .get("mode")
                .and_then(Value::as_str)
                .unwrap_or("central")
                .parse::<PeaksMode>()
                .unwrap(),
            min_window: a.get("min.window").map(num).unwrap_or(10.0) as usize,
            max_window: a.get("max.window").map(|v| num(v) as usize),
            window_step: a.get("window.step").map(num).unwrap_or(10.0) as usize,
            position_step: a.get("position.step").map(num).unwrap_or(5.0) as usize,
            nthreads: 2,
        };
        let qmax = num(&a["qvalue"]);
        let mut got: Vec<_> = motif_peaks(&hits, 300, &opts)
            .unwrap()
            .into_iter()
            .filter(|p| p.qvalue <= qmax)
            .collect();
        got.sort_by(|x, y| x.pvalue.total_cmp(&y.pvalue).then(x.motif.cmp(&y.motif)));
        let t = &case["table"];
        let names: Vec<String> = got.iter().map(|p| sorted[p.motif].clone()).collect();
        assert_eq!(names, texts(&t["motif"]), "{ctx}: motif");
        let idx: Vec<f64> = got.iter().map(|p| p.motif as f64 + 1.0).collect();
        assert_eq!(idx, nums(&t["motif.i"]), "{ctx}: motif.i");
        let u = |f: fn(&universalmotif::peaks::MotifPeak) -> usize| {
            got.iter().map(|p| f(p) as f64).collect::<Vec<_>>()
        };
        assert_eq!(u(|p| p.nhits), nums(&t["nhits"]), "{ctx}: nhits");
        assert_eq!(u(|p| p.hits_in), nums(&t["hits.in"]), "{ctx}: hits.in");
        assert_eq!(u(|p| p.hits_out), nums(&t["hits.out"]), "{ctx}: hits.out");
        let o = |f: fn(&universalmotif::peaks::MotifPeak) -> Option<f64>| {
            got.iter().map(|p| opt(f(p))).collect::<Vec<_>>()
        };
        assert_close_vec(
            &o(|p| p.best_window.map(|w| w as f64)),
            &nums(&t["best.window"]),
            &format!("{ctx}: best.window"),
        );
        assert_close_vec(
            &o(|p| p.best_centre.map(|c| c + 1.0)),
            &nums(&t["best.center"]),
            &format!("{ctx}: centre"),
        );
        assert_close_vec(
            &o(|p| p.expected_in),
            &nums(&t["expected.in"]),
            &format!("{ctx}: exp"),
        );
        assert_close_vec(
            &o(|p| p.enrichment),
            &nums(&t["enrichment"]),
            &format!("{ctx}: enr"),
        );
        assert_close_vec(
            &o(|p| p.log2_enrichment),
            &nums(&t["log2.enrichment"]),
            &format!("{ctx}: log2"),
        );
        let f =
            |g: fn(&universalmotif::peaks::MotifPeak) -> f64| got.iter().map(g).collect::<Vec<_>>();
        assert_close_vec(
            &f(|p| p.pvalue),
            &nums(&t["pvalue"]),
            &format!("{ctx}: pvalue"),
        );
        assert_close_vec(
            &f(|p| p.qvalue),
            &nums(&t["qvalue"]),
            &format!("{ctx}: qvalue"),
        );
        for (p, w) in got.iter().zip(t["centers"].as_array().unwrap()) {
            let c: Vec<f64> = p.centres.iter().map(|&x| x as f64 + 1.0).collect();
            assert_eq!(c, nums(w), "{ctx}: centres");
        }
    }
}
