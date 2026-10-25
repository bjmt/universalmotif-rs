//! Randomized motifs, sequences and arguments checked against reference
//! values from the R package, exported by scripts/make_reference_random.R.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::compare::{CompareOptions, Comparison, NullModel, compare_matrix};
use universalmotif::complexity::{
    ComplexityMethod, ComplexityOptions, complexity, sequence_complexity,
};
use universalmotif::convert::convert_type;
use universalmotif::enrich::{EnrichOptions, EnrichTest, Enrichment, enrich_motifs};
use universalmotif::kmer::{BkgOptions, WindowOptions, get_bkg};
use universalmotif::merge::{MergeOptions, merge_motifs};
use universalmotif::metrics::{Average, Column, ColumnMetric, average_ic, compare_columns};
use universalmotif::ops::{TrimFrom, reverse_complement, round_motif, trim_motif};
use universalmotif::pvalue::{PvalueOptions, motif_pvalue, pvalue_to_score};
use universalmotif::scan::{NoOverlaps, OverlapPriority, ScanOptions, scan_sequences, scan_table};
use universalmotif::scoring::{
    ThresholdType, get_matches, get_scores, motif_range, motif_score, prob_match, score_match,
};
use universalmotif::{Motif, MotifType, Strand};

fn reference(file: &str) -> Value {
    let path = format!("{}/tests/reference/{file}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn cases<'a>(r: &'a Value, key: &str) -> &'a Vec<Value> {
    r[key].as_array().unwrap()
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

fn idx(v: &Value) -> Vec<usize> {
    nums(v).into_iter().map(|x| x as usize).collect()
}

fn all_motifs() -> Vec<Motif> {
    cases(&reference("random_motifs.json"), "motifs")
        .iter()
        .map(motif_from_record)
        .collect()
}

fn pick(motifs: &[Motif], v: &Value) -> Vec<Motif> {
    idx(v).into_iter().map(|i| motifs[i].clone()).collect()
}

/// The total IC the port gives a motif converted with a new pseudocount:
/// that of the input motif with the new pseudocount.
fn icscore_with_pseudocount(m: &Motif, pseudocount: f64) -> f64 {
    let mut parts = m.clone().into_parts();
    parts.pseudocount = pseudocount;
    own_icscore(&Motif::from_parts(parts).unwrap())
}

#[test]
fn convert_type_matches_r() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "convert").iter().enumerate() {
        let input = &motifs[num(&case["input"]) as usize];
        let a = &case["args"];
        let ctx = format!(
            "convert case {i} ({} {} -> {a})",
            input.name, input.motif_type
        );
        let to = motif_type(a["type"].as_str().unwrap());
        let pseudo = opt_num(&a["pseudocount"]);
        let got = convert_type(
            input,
            to,
            pseudo,
            a["relative_entropy"].as_bool().unwrap_or(false),
        )
        .unwrap_or_else(|e| panic!("{ctx}: {e}"));
        let want = &case["motif"];
        // With a new pseudocount the R package now recomputes the total IC
        // from the new matrix as its new type, while the port gives the total
        // IC of the input with the new pseudocount. The consensus is checked
        // in `consensus_follows_a_new_pseudocount`.
        let new_pseudo = pseudo.filter(|&p| p != input.pseudocount && to != input.motif_type);
        let skip = Skip {
            derived: new_pseudo.is_some(),
            ..Default::default()
        };
        assert_motif_eq_skip(&got, want, &ctx, skip);
        if let Some(p) = new_pseudo {
            let ic = icscore_with_pseudocount(input, p);
            assert!(close(got.icscore, ic), "{ctx}: icscore");
            if to != MotifType::Pwm {
                assert_eq!(
                    got.consensus,
                    opt_str(&want["consensus"]),
                    "{ctx}: consensus"
                );
            }
        } else {
            // Checked here as well, since the shared helper leaves out the
            // total IC of PPMs that R used to round to counts differently.
            assert!(close(got.icscore, num(&want["icscore"])), "{ctx}: icscore");
        }
    }
}

/// A new pseudocount changes the probabilities of a PCM or PPM, so the
/// consensus of the converted motif can change, as it does in the R package.
/// The port keeps the input's consensus. PWMs are left out: their consensus
/// in the R package ignores a non-uniform background.
#[test]
fn consensus_follows_a_new_pseudocount() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "convert").iter().enumerate() {
        let input = &motifs[num(&case["input"]) as usize];
        let a = &case["args"];
        let to = motif_type(a["type"].as_str().unwrap());
        let Some(p) = opt_num(&a["pseudocount"]) else {
            continue;
        };
        if p == input.pseudocount
            || to == input.motif_type
            || to == MotifType::Pwm
            || !matches!(input.motif_type, MotifType::Pcm | MotifType::Ppm)
        {
            continue;
        }
        let got = convert_type(input, to, Some(p), false).unwrap();
        assert_eq!(
            got.consensus,
            opt_str(&case["motif"]["consensus"]),
            "convert case {i} ({} {} -> {a})",
            input.name,
            input.motif_type
        );
    }
}

#[test]
fn motif_operations_match_r() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "ops").iter().enumerate() {
        let op = case["op"].as_str().unwrap();
        let input = &motifs[num(&case["input"]) as usize];
        let a = &case["args"];
        let ctx = format!(
            "ops case {i} ({op} on {} {}, {a})",
            input.name, input.motif_type
        );
        let mut skip = Skip::default();
        let got = match op {
            "motif_rc" => reverse_complement(input, false).map(Some),
            "trim_motifs" => {
                // The R package stores nsites = 100 in trimmed motifs that had
                // none (port-bugs.md, "trim_motifs() invents nsites").
                skip.nsites = input.nsites.is_none();
                let tf = a["trim_from"]
                    .as_str()
                    .unwrap()
                    .parse::<TrimFrom>()
                    .unwrap();
                trim_motif(input, num(&a["min_ic"]), tf)
            }
            "round_motif" => {
                // The R package keeps the old total IC and consensus.
                skip.derived = true;
                round_motif(input, num(&a["pct_tolerance"])).map(Some)
            }
            _ => panic!("unknown op {op}"),
        };
        match (&case["motif"], got) {
            (Value::Null, Ok(None)) => {
                let e = case["error"].as_str().unwrap();
                assert!(e.contains("completely trimmed"), "{ctx}: {e}");
            }
            // A column with every probability below the tolerance has
            // nothing left to renormalize, which both refuse.
            (Value::Null, Err(_)) => {
                assert_eq!(op, "round_motif", "{ctx}");
                assert!(case["error"].as_str().unwrap().contains("NaN"), "{ctx}");
            }
            (Value::Null, Ok(Some(m))) => {
                panic!("{ctx}: R failed with {}, but got\n{m}", case["error"])
            }
            (_, Err(e)) => panic!("{ctx}: {e}"),
            (_, Ok(None)) => panic!("{ctx}: unexpectedly trimmed away"),
            (want, Ok(Some(m))) => {
                assert_motif_eq_skip(&m, want, &ctx, skip);
                if skip.derived {
                    assert!(close(m.icscore, own_icscore(&m)), "{ctx}: icscore");
                } else {
                    assert!(close(m.icscore, num(&want["icscore"])), "{ctx}: icscore");
                }
                if skip.nsites {
                    assert_eq!(m.nsites, None, "{ctx}: nsites");
                }
            }
        }
    }
}

#[test]
fn average_ic_matches_r() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for case in cases(&r, "average_ic") {
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
fn scoring_helpers_match_r() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "scoring").iter().enumerate() {
        let f = case["fn"].as_str().unwrap();
        let m = &motifs[num(&case["motif"]) as usize];
        let a = &case["args"];
        let ctx = format!("scoring case {i} ({f} on {} {}, {a})", m.name, m.motif_type);
        let want = &case["value"];
        match f {
            "motif_range" => {
                let (lo, hi) = motif_range(m, false).unwrap();
                assert_close_vec(&[lo, hi], &nums(want), &ctx);
            }
            "motif_score" => {
                let t = a["threshold_type"]
                    .as_str()
                    .unwrap()
                    .parse::<ThresholdType>()
                    .unwrap();
                let got = motif_score(m, &nums(&a["threshold"]), t, false).unwrap();
                assert_close_vec(&got, &nums(want), &ctx);
            }
            "score_match" => {
                // The R package adds the truncated entries as floats and
                // truncates the total again, so it can be 0.001 off the
                // integer sum the port uses (port-bugs.md, "score_match() can
                // be 0.001 off the scan score").
                let got = score_match(m, &texts(&a["matches"]), false).unwrap();
                let want = nums(want);
                assert_eq!(got.len(), want.len(), "{ctx}");
                for (g, w) in got.iter().zip(&want) {
                    assert!(
                        g == w || (g - w).abs() <= 0.001 + 1e-9,
                        "{ctx}: got {g}, want {w}"
                    );
                }
            }
            "prob_match" => {
                let az = a["allow_zero"].as_bool().unwrap();
                let got = prob_match(m, &texts(&a["matches"]), az).unwrap();
                assert_close_vec(&got, &nums(want), &ctx);
            }
            "get_scores" => {
                assert_close_vec(&get_scores(m, false).unwrap(), &nums(want), &ctx);
            }
            "get_matches" => {
                let s = num(&a["score"]);
                let got = get_matches(m, s, false).unwrap();
                if case["error"].is_string() {
                    // The R package fails for single-column motifs, where
                    // its score matrix collapses to a vector; every letter
                    // scoring at least `s` is checked for instead.
                    assert_eq!(m.width(), 1, "{ctx}: {}", case["error"]);
                    let letters: Vec<String> =
                        m.alphabet.letters().iter().map(|c| c.to_string()).collect();
                    let scores = score_match(m, &letters, false).unwrap();
                    let mut want: Vec<(f64, String)> = scores
                        .into_iter()
                        .zip(letters)
                        .filter(|(sc, _)| *sc >= s - 1e-9)
                        .collect();
                    want.sort_by(|x, y| x.1.cmp(&y.1));
                    let mut got = got.clone();
                    got.sort();
                    let want: Vec<String> = want.into_iter().map(|w| w.1).collect();
                    assert_eq!(got, want, "{ctx}");
                } else {
                    let mut got = got;
                    let mut want = texts(want);
                    got.sort();
                    want.sort();
                    assert_eq!(got, want, "{ctx}");
                }
            }
            _ => panic!("unknown function {f}"),
        }
    }
}

#[test]
fn pvalues_and_scores_match_r() {
    let r = reference("random_motifs.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "pvalue").iter().enumerate() {
        let m = &motifs[num(&case["motif"]) as usize];
        let ctx = format!(
            "pvalue case {i} ({} {}, bkg {})",
            m.name, m.motif_type, case["bkg"]
        );
        let opts = PvalueOptions {
            bkg: if case["bkg"].is_null() {
                Vec::new()
            } else {
                vec![nums(&case["bkg"])]
            },
            ..Default::default()
        };
        let ms = std::slice::from_ref(m);
        let q: Vec<(usize, f64)> = nums(&case["scores"]).into_iter().map(|s| (0, s)).collect();
        let got = motif_pvalue(ms, &q, &opts).unwrap();
        assert_close_vec(&got, &nums(&case["pvalues"]), &format!("{ctx}: P-values"));
        let q: Vec<(usize, f64)> = nums(&case["pvalue_in"])
            .into_iter()
            .map(|p| (0, p))
            .collect();
        let got = pvalue_to_score(ms, &q, &opts).unwrap();
        assert_close_vec(&got, &nums(&case["score_out"]), &format!("{ctx}: scores"));
    }
}

// The R script scans PCMs with zero counts and no pseudocount as PPMs, since
// scan_sequences_lite() reports hits at nearly every position for them; the
// port scans the PCMs as they are.
#[test]
fn scans_match_r() {
    let r = reference("random_scan.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "scan").iter().enumerate() {
        let a = args(&case["args"]);
        let set = case["set"].as_str().unwrap();
        let ctx = format!("scan case {i} ({set} {a:?})");
        let ms = pick(&motifs, &case["motifs"]);
        let seqs = texts(&r[format!("{set}_sequences")]);
        let b = |k: &str| a.get(k).and_then(Value::as_bool).unwrap_or(false);
        let opts = ScanOptions {
            pvalue: num(&a["pvalue"]),
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
        let hits = scan_sequences(&ms, &seqs, &opts).unwrap();
        let w = &case["hits"];
        assert_eq!(hits.len(), idx(&w["start"]).len(), "{ctx}: number of hits");
        let ints = |f: &dyn Fn(&universalmotif::scan::Hit) -> usize| {
            hits.iter().map(f).collect::<Vec<_>>()
        };
        assert_eq!(ints(&|h| h.motif + 1), idx(&w["motif_i"]), "{ctx}: motif");
        assert_eq!(
            ints(&|h| h.sequence + 1),
            idx(&w["sequence_i"]),
            "{ctx}: sequence"
        );
        // R counts positions from 1 and includes the end.
        assert_eq!(ints(&|h| h.start + 1), idx(&w["start"]), "{ctx}: start");
        assert_eq!(ints(&|h| h.end), idx(&w["end"]), "{ctx}: end");
        let strands: Vec<String> = hits
            .iter()
            .map(|h| if h.strand == Strand::Minus { "-" } else { "+" }.to_string())
            .collect();
        assert_eq!(strands, texts(&w["strand"]), "{ctx}: strand");
        let matched: Vec<String> = hits.iter().map(|h| h.matched.clone()).collect();
        assert_eq!(matched, texts(&w["match"]), "{ctx}: match");
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
fn protein_and_custom_scans_match_r() {
    let r = reference("random_scan.json");
    let motifs = all_motifs();
    for case in cases(&r, "full_scan") {
        let set = case["set"].as_str().unwrap();
        let threshold = num(&case["threshold"]);
        let ctx = format!("{set} scan at {threshold}");
        let ms = pick(&motifs, &case["motifs"]);
        let seqs = texts(&r[format!("{set}_sequences")]);
        let opts = ScanOptions {
            pvalue: threshold,
            nthreads: 1,
            ..Default::default()
        };
        let t = scan_table(&ms, &seqs, &opts).unwrap();
        let w = &case["table"];
        assert_eq!(t.len(), nums(&w["start"]).len(), "{ctx}: number of hits");
        let plus1 = |v: &[i64]| v.iter().map(|&x| x as f64 + 1.0).collect::<Vec<_>>();
        assert_eq!(plus1(&t.motif), nums(&w["motif.i"]), "{ctx}: motif");
        assert_eq!(
            plus1(&t.sequence),
            nums(&w["sequence.i"]),
            "{ctx}: sequence"
        );
        assert_eq!(plus1(&t.start), nums(&w["start"]), "{ctx}: start");
        let end: Vec<f64> = t.end.iter().map(|&x| x as f64).collect();
        assert_eq!(end, nums(&w["stop"]), "{ctx}: end");
        assert!(t.minus.iter().all(|m| !m), "{ctx}: only the forward strand");
        assert_close_vec(&t.score, &nums(&w["score"]), &format!("{ctx}: score"));
        let matched: Vec<String> = (0..t.len()).map(|i| t.matched(i, &seqs)).collect();
        assert_eq!(matched, texts(&w["match"]), "{ctx}: match");
        assert_close_vec(&t.pvalue, &nums(&w["pvalue"]), &format!("{ctx}: pvalue"));
    }
}

#[test]
fn enrichment_matches_r() {
    let r = reference("random_scan.json");
    let motifs = pick(&all_motifs(), &r["enrich_motifs"]);
    let tgt = texts(&r["enrich_targets"]);
    let bkg = texts(&r["enrich_background"]);
    for case in cases(&r, "enrich") {
        let a = args(&case["args"]);
        let ctx = format!("enrich {a:?}");
        let opts = EnrichOptions {
            pvalue: num(&a["pvalue"]),
            test: a
                .get("test")
                .and_then(Value::as_str)
                .unwrap_or("seqs")
                .parse::<EnrichTest>()
                .unwrap(),
            rc: a.get("RC").and_then(Value::as_bool).unwrap_or(true),
            pseudocount: a.get("pseudocount").map(num).unwrap_or(1.0) as u64,
            nthreads: 1,
            ..Default::default()
        };
        let mut got = enrich_motifs(&motifs, &tgt, Some(&bkg), &opts).unwrap();
        got.sort_by(|x, y| x.pvalue.total_cmp(&y.pvalue).then(x.motif.cmp(&y.motif)));
        let t = &case["table"];
        let names: Vec<String> = got.iter().map(|e| motifs[e.motif].name.clone()).collect();
        assert_eq!(names, texts(&t["motif"]), "{ctx}: motifs");
        let counts =
            |f: fn(&Enrichment) -> usize| got.iter().map(|e| f(e) as f64).collect::<Vec<_>>();
        assert_eq!(
            counts(|e| e.target_seq_n),
            nums(&t["target.seq.n"]),
            "{ctx}"
        );
        assert_eq!(
            counts(|e| e.target_seq_hits),
            nums(&t["target.seq.hits"]),
            "{ctx}"
        );
        assert_eq!(
            counts(|e| e.target_site_hits),
            nums(&t["target.site.hits"]),
            "{ctx}"
        );
        assert_eq!(counts(|e| e.bkg_seq_n), nums(&t["bkg.seq.n"]), "{ctx}");
        assert_eq!(
            counts(|e| e.bkg_seq_hits),
            nums(&t["bkg.seq.hits"]),
            "{ctx}"
        );
        assert_eq!(
            counts(|e| e.bkg_site_hits),
            nums(&t["bkg.site.hits"]),
            "{ctx}"
        );
        let f = |g: fn(&Enrichment) -> f64| got.iter().map(g).collect::<Vec<_>>();
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

/// The R package gives a P-value of 0 to scores beyond what the parametric
/// null can reach, where the port gives the null's highest score's
/// (port-bugs.md, "The parametric null gives P-values of exactly 0"); such
/// values are checked to be small and then set to R's 0.
fn beyond_null(got: &[f64], want: &[f64], null: NullModel) -> Vec<f64> {
    got.iter()
        .zip(want)
        .map(|(&g, &w)| {
            if null == NullModel::Parametric && w == 0.0 && g != 0.0 {
                assert!(g > 0.0 && g < 1e-6, "{g} for R's 0");
                0.0
            } else {
                g
            }
        })
        .collect()
}

#[test]
fn comparisons_match_r() {
    let r = reference("random_compare.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "compare").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("compare case {i} ({} {a:?})", case["set"]);
        let ms = pick(&motifs, &case["motifs"]);
        let opts = compare_opts(&a);
        if let Some(rows) = case["matrix"].as_array() {
            let out = a
                .get("matrix.out")
                .and_then(Value::as_str)
                .unwrap_or("score");
            let m = compare_matrix(&ms, &opts, out != "score").unwrap();
            let got = match out {
                "score" => m.scores,
                "pvalue" => m.pvalues,
                _ => m.qvalues,
            };
            for (k, (g, w)) in got.iter().zip(rows).enumerate() {
                let g = beyond_null(g, &nums(w), opts.null);
                assert_close_vec(&g, &nums(w), &format!("{ctx}: row {k}"));
            }
        } else {
            let w = &case["long"];
            let queries: Vec<usize> = idx(&a["compare.to"]).iter().map(|&x| x - 1).collect();
            let mut rows: Vec<Comparison> =
                universalmotif::compare::compare_motifs(&ms, &queries, &opts)
                    .unwrap()
                    .into_iter()
                    .filter(|c| c.overlap > 0)
                    .collect();
            rows.sort_by(|x, y| {
                x.qvalue
                    .total_cmp(&y.qvalue)
                    .then(x.pvalue.total_cmp(&y.pvalue))
                    .then(x.target.cmp(&y.target))
            });
            let ints = |f: &dyn Fn(&Comparison) -> i64| rows.iter().map(f).collect::<Vec<_>>();
            let want_ints = |k: &str| nums(&w[k]).iter().map(|&x| x as i64).collect::<Vec<_>>();
            assert_eq!(
                ints(&|c| c.query as i64 + 1),
                want_ints("subject_i"),
                "{ctx}"
            );
            assert_eq!(
                ints(&|c| c.target as i64 + 1),
                want_ints("target_i"),
                "{ctx}"
            );
            assert_eq!(ints(&|c| c.offset), want_ints("offset"), "{ctx}: offset");
            assert_eq!(
                ints(&|c| c.overlap as i64),
                want_ints("overlap"),
                "{ctx}: overlap"
            );
            let strands: Vec<String> = rows.iter().map(|c| c.strand.as_str().to_string()).collect();
            assert_eq!(strands, texts(&w["strand"]), "{ctx}: strand");
            let f = |g: &dyn Fn(&Comparison) -> f64| rows.iter().map(g).collect::<Vec<_>>();
            assert_close_vec(
                &f(&|c| c.score),
                &nums(&w["score"]),
                &format!("{ctx}: score"),
            );
            assert_close_vec(
                &beyond_null(&f(&|c| c.pvalue), &nums(&w["pvalue"]), opts.null),
                &nums(&w["pvalue"]),
                &format!("{ctx}: pvalue"),
            );
            assert_close_vec(
                &beyond_null(&f(&|c| c.qvalue), &nums(&w["qvalue"]), opts.null),
                &nums(&w["qvalue"]),
                &format!("{ctx}: qvalue"),
            );
        }
    }
}

#[test]
fn protein_and_custom_comparisons_match_r() {
    let r = reference("random_compare.json");
    let motifs = all_motifs();
    for case in cases(&r, "full_compare") {
        let ms = pick(&motifs, &case["motifs"]);
        let n = ms.len();
        let opts = CompareOptions {
            min_overlap: num(&case["min_overlap"]) as usize,
            ..Default::default()
        };
        let ctx = format!("{} compare, min_overlap {}", case["set"], opts.min_overlap);
        let rows = universalmotif::compare::compare_motifs(&ms, &(0..n).collect::<Vec<_>>(), &opts)
            .unwrap();
        let want: Vec<Vec<f64>> = cases(case, "sums").iter().map(nums).collect();
        for c in &rows {
            let sum = c.score * c.overlap as f64;
            let w = want[c.query][c.target];
            assert!(
                close(sum, w) || (sum - w).abs() < 1e-12,
                "{ctx} {}x{}: {sum} vs {w}",
                c.query,
                c.target
            );
        }
    }
}

#[test]
fn merges_match_r() {
    let r = reference("random_compare.json");
    let motifs = all_motifs();
    for (i, case) in cases(&r, "merge").iter().enumerate() {
        let a = args(&case["args"]);
        let ctx = format!("merge case {i} ({:?} {a:?})", idx(&case["motifs"]));
        let ms = pick(&motifs, &case["motifs"]);
        let opts = MergeOptions {
            min_overlap: num(&a["min.overlap"]) as usize,
            rc: a["RC"].as_bool().unwrap(),
            weighted: a["weighted"].as_bool().unwrap(),
            nthreads: 1,
            ..Default::default()
        };
        let got = merge_motifs(&ms, &opts).unwrap_or_else(|e| panic!("{ctx}: {e}"));
        assert_motif_eq(&got, &case["motif"], &ctx);
    }
}

/// The background-weighted Pearson correlation with weights that sum to 1,
/// which the R package does not normalize (port-bugs.md, "WPCC weights are
/// not normalized").
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
    if bot == 0.0 {
        1.0
    } else {
        cov(x, y, mx, my) / bot
    }
}

#[test]
fn column_metrics_match_r() {
    let r = reference("random_compare.json");
    for case in cases(&r, "compare_columns") {
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
        )
        .unwrap();
        let ctx = format!("{metric} {x:?} {y:?}");
        let want = num(&case["value"]);
        if metric == ColumnMetric::Wpcc {
            assert!(
                close(got, weighted_pcc(&x, &y, &bkg1, &bkg2)),
                "{ctx}: {got}"
            );
        } else if want.is_nan() {
            // The R package's 0 * log(0) is NaN where the port takes it as 0
            // (port-bugs.md, "KL and ALLR give NaN for zero probabilities").
            assert!(
                matches!(
                    metric,
                    ColumnMetric::Kl | ColumnMetric::Allr | ColumnMetric::AllrLl
                ),
                "{ctx}: R gave NaN"
            );
            assert!(!got.is_nan(), "{ctx}");
        } else {
            assert!(close(got, want), "{ctx}: got {got}, want {want}");
        }
    }
}

#[test]
fn backgrounds_match_r() {
    let r = reference("random_sequences.json");
    let seqs = texts(&r["sequences"]);
    let names = texts(&r["names"]);
    for case in cases(&r, "bkg") {
        let a = args(&case["args"]);
        let ctx = format!("get_bkg {a:?}");
        let mut opts = BkgOptions::default();
        if let Some(k) = a.get("k") {
            opts.k = idx(k);
        }
        if let Some(p) = a.get("pseudocount") {
            opts.pseudocount = num(p) as u64;
        }
        opts.rc = a.get("RC").and_then(Value::as_bool).unwrap_or(false);
        opts.merge = a.get("merge.res").and_then(Value::as_bool).unwrap_or(true);
        if a.get("window").and_then(Value::as_bool) == Some(true) {
            opts.window = Some(WindowOptions {
                size: num(&a["window.size"]),
                overlap: num(&a["window.overlap"]),
            });
        }
        let got = get_bkg(&seqs, &names, &opts);
        if case["error"].is_string() {
            assert!(got.is_err(), "{ctx}: R failed with {}", case["error"]);
            continue;
        }
        let got = got.unwrap();
        let t = &case["table"];
        let kmer: Vec<String> = got.iter().map(|r| r.kmer.clone()).collect();
        assert_eq!(kmer, texts(&t["klet"]), "{ctx}: kmer");
        let count: Vec<f64> = got.iter().map(|r| r.count as f64).collect();
        assert_eq!(count, nums(&t["count"]), "{ctx}: count");
        let prob: Vec<f64> = got.iter().map(|r| r.probability).collect();
        assert_close_vec(&prob, &nums(&t["probability"]), &format!("{ctx}: prob"));
        if !t["sequence"].is_null() {
            let sq: Vec<String> = got
                .iter()
                .map(|r| names[r.sequence.unwrap() % names.len()].clone())
                .collect();
            assert_eq!(sq, texts(&t["sequence"]), "{ctx}: sequence");
        }
        if !t["start"].is_null() {
            let start: Vec<usize> = got.iter().map(|r| r.start.unwrap() + 1).collect();
            let end: Vec<usize> = got.iter().map(|r| r.end.unwrap()).collect();
            assert_eq!(start, idx(&t["start"]), "{ctx}: start");
            assert_eq!(end, idx(&t["stop"]), "{ctx}: stop");
        }
    }
}

#[test]
fn complexity_matches_r() {
    let r = reference("random_sequences.json");
    for case in cases(&r, "complexity") {
        let s = case["string"].as_str().unwrap();
        let method = case["method"]
            .as_str()
            .unwrap()
            .parse::<ComplexityMethod>()
            .unwrap();
        let max_word = num(&case["max_word"]) as usize;
        let got = complexity(s, method, Some("ACGT"), max_word);
        let want = num(&case["value"]);
        assert!(
            close(got, want),
            "{s} {method:?} {max_word}: {got} vs {want}"
        );
    }
    let seqs = texts(&r["sequences"]);
    let names = texts(&r["names"]);
    for case in cases(&r, "sequence_complexity") {
        let a = args(&case["args"]);
        let mut opts = ComplexityOptions::default();
        if let Some(m) = a.get("method") {
            opts.method = m.as_str().unwrap().parse::<ComplexityMethod>().unwrap();
        }
        if let Some(w) = a.get("window.size") {
            opts.window = num(w);
        }
        if let Some(o) = a.get("window.overlap") {
            opts.overlap = Some(num(o));
        }
        if let Some(m) = a.get("trifonov.max.word.size") {
            opts.max_word = num(m) as usize;
        }
        let ctx = format!("sequence_complexity {a:?}");
        let got = sequence_complexity(&seqs, &opts).unwrap();
        let t = &case["table"];
        let seqnames: Vec<String> = got.iter().map(|r| names[r.sequence].clone()).collect();
        assert_eq!(seqnames, texts(&t["sequence"]), "{ctx}: sequence");
        let start: Vec<usize> = got.iter().map(|r| r.start + 1).collect();
        let end: Vec<usize> = got.iter().map(|r| r.end).collect();
        assert_eq!(start, idx(&t["start"]), "{ctx}: start");
        assert_eq!(end, idx(&t["stop"]), "{ctx}: stop");
        let cx: Vec<f64> = got.iter().map(|r| r.complexity).collect();
        assert_close_vec(&cx, &nums(&t["complexity"]), &format!("{ctx}: complexity"));
    }
}
