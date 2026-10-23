//! Motif file readers and writers against universalmotif.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::io::{cisbp, homer, jaspar, matrix, meme, transfac, uniprobe};
use universalmotif::{Motif, Result};

/// Strings from R. `read.table()` turns numeric-looking columns (such as
/// sequence names "1", "2") into numbers, so those are accepted too.
fn strs(v: &Value) -> Vec<String> {
    let one = |s: &Value| match s {
        Value::String(s) => s.clone(),
        other => other.to_string(),
    };
    match v {
        Value::Array(a) => a.iter().map(one).collect(),
        Value::String(s) => vec![s.clone()],
        _ => Vec::new(),
    }
}

fn matrix_headers(a: &serde_json::Map<String, Value>) -> matrix::Headers {
    match a.get("headers") {
        Some(Value::Bool(false)) => matrix::Headers::None,
        Some(Value::String(p)) => matrix::Headers::Prefix(p.clone()),
        _ => matrix::Headers::FirstLine,
    }
}

fn positions(a: &serde_json::Map<String, Value>) -> matrix::Positions {
    match a.get("positions").and_then(Value::as_str) {
        Some("rows") => matrix::Positions::Rows,
        _ => matrix::Positions::Columns,
    }
}

fn read_case(case: &Value) -> Result<Vec<Motif>> {
    let text = reference_file(case["file"].as_str().unwrap());
    let a = args(&case["args"]);
    match case["format"].as_str().unwrap() {
        "jaspar" => jaspar::read(&text),
        "homer" => homer::read(&text),
        "transfac" => transfac::read(&text),
        "cisbp" => cisbp::read(&text),
        "uniprobe" => uniprobe::read(&text),
        "meme" => meme::read(
            &text,
            &meme::ReadOptions {
                cwm: a.get("cwm").and_then(Value::as_bool).unwrap_or(false),
            },
        ),
        "matrix" => matrix::read(
            &text,
            &matrix::ReadOptions {
                positions: positions(&a),
                headers: matrix_headers(&a),
                rownames: a.get("rownames").and_then(Value::as_bool).unwrap_or(false),
                motif_type: a.get("type").and_then(Value::as_str).map(motif_type),
                ..Default::default()
            },
        ),
        f => panic!("unknown format {f}"),
    }
}

#[test]
fn readers_match_r() {
    for (i, case) in load_key("io.json", "read").iter().enumerate() {
        let ctx = format!(
            "read case {i} ({} {} {})",
            case["format"], case["file"], case["args"]
        );
        let got = read_case(case);
        let Some(want) = case["motifs"].as_array() else {
            assert!(
                got.is_err(),
                "{ctx}: R failed with {}, the crate did not",
                case["error"]
            );
            continue;
        };
        let got = got.unwrap_or_else(|e| panic!("{ctx}: {e}"));
        assert_eq!(got.len(), want.len(), "{ctx}: number of motifs");
        for (k, (g, w)) in got.iter().zip(want).enumerate() {
            assert_motif_eq(g, w, &format!("{ctx}, motif {k}"));
        }

        let text = reference_file(case["file"].as_str().unwrap());
        if let Some(sites) = case["sites"].as_array() {
            let got = meme::read_sites(&text).unwrap();
            assert_eq!(got.len(), sites.len(), "{ctx}: sites");
            for (g, w) in got.iter().zip(sites) {
                let names: Vec<&str> = g.sites.iter().map(|s| s.sequence.as_str()).collect();
                let seqs: Vec<&str> = g.sites.iter().map(|s| s.site.as_str()).collect();
                assert_eq!(names, strs(&w["names"]), "{ctx}: site names");
                assert_eq!(seqs, strs(&w["sites"]), "{ctx}: sites");
            }
        }
        if let Some(meta) = case["sites_meta"].as_array() {
            let got = meme::read_site_pvalues(&text).unwrap();
            assert_eq!(got.len(), meta.len(), "{ctx}: site tables");
            for (g, w) in got.iter().zip(meta) {
                let seq: Vec<&str> = g.sites.iter().map(|s| s.sequence.as_str()).collect();
                let site: Vec<&str> = g.sites.iter().map(|s| s.site.as_str()).collect();
                // R counts positions from 1.
                let pos: Vec<f64> = g.sites.iter().map(|s| s.position + 1.0).collect();
                let pv: Vec<f64> = g.sites.iter().map(|s| s.pvalue).collect();
                assert_eq!(seq, strs(&w["sequence"]), "{ctx}: site table sequences");
                assert_eq!(site, strs(&w["site"]), "{ctx}: site table sites");
                assert_close_vec(&pos, &nums(&w["position"]), &format!("{ctx}: positions"));
                assert_close_vec(&pv, &nums(&w["pvalue"]), &format!("{ctx}: P-values"));
            }
        }
        if case["combined"].is_object() {
            let w = &case["combined"];
            let got = meme::read_combined(&text).unwrap();
            let seq: Vec<&str> = got.iter().map(|s| s.sequence.as_str()).collect();
            let pv: Vec<f64> = got.iter().map(|s| s.pvalue).collect();
            assert_eq!(seq, strs(&w["sequence"]), "{ctx}: combined sequences");
            assert_close_vec(
                &pv,
                &nums(&w["pvalue"]),
                &format!("{ctx}: combined P-values"),
            );
            // universalmotif drops the third and later lines of a diagram,
            // so its diagrams can be truncated.
            for (g, d) in got.iter().zip(strs(&w["diagram"])) {
                assert!(
                    g.diagram == d || (g.diagram.starts_with(&d) && d.ends_with('_')),
                    "{ctx}: diagram {} vs {d}",
                    g.diagram
                );
            }
        }
    }
}

#[test]
fn combined_diagrams_span_lines() {
    let text = reference_file("meme_full.txt");
    let got = meme::read_combined(&text).unwrap();
    let row = got.iter().find(|s| s.sequence == "AT1G56240").unwrap();
    assert_eq!(
        row.diagram,
        "38_[+2(1.30e-06)]_26_[+3(8.51e-05)]_18_[-1(5.87e-06)]_323_[+1(1.39e-06)]_\
         234_[+3(2.60e-05)]_172_[+3(6.10e-06)]_97"
    );
}

fn write_case(case: &Value, motifs: &[Motif]) -> Result<String> {
    let a = args(&case["args"]);
    let b = |k: &str| a.get(k).and_then(Value::as_bool).unwrap_or(false);
    match case["format"].as_str().unwrap() {
        "meme" => meme::write(
            motifs,
            &meme::WriteOptions {
                version: a.get("version").map(|v| num(v) as u32).unwrap_or(5),
                bkg: a.get("bkg").map(nums),
                strand: a.get("strand").and_then(Value::as_str).map(String::from),
                cwm: b("cwm"),
                header: true,
            },
        ),
        "jaspar" => jaspar::write(motifs),
        "homer" => {
            let t = a.get("threshold").map(nums).unwrap_or(vec![0.8]);
            let th = match a.get("threshold_type").and_then(Value::as_str) {
                Some("logodds_abs") => homer::Threshold::LogoddsAbs(t),
                _ => homer::Threshold::Logodds(t),
            };
            homer::write(motifs, &th)
        }
        "transfac" => {
            let mut o = transfac::WriteOptions::default();
            if let Some(t) = a.get("name_tag").and_then(Value::as_str) {
                o.name_tag = t.into();
            }
            if let Some(t) = a.get("altname_tag").and_then(Value::as_str) {
                o.altname_tag = t.into();
            }
            transfac::write(motifs, &o)
        }
        "matrix" => matrix::write(
            motifs,
            &matrix::WriteOptions {
                positions: positions(&a),
                rownames: b("rownames"),
                motif_type: a.get("type").and_then(Value::as_str).map(motif_type),
                sep: if b("sep_none") {
                    None
                } else {
                    Some(
                        a.get("sep")
                            .and_then(Value::as_str)
                            .unwrap_or("")
                            .to_string(),
                    )
                },
                headers: matrix_headers(&a),
                digits: a.get("digits").map(|v| num(v) as usize).unwrap_or(6),
            },
        ),
        f => panic!("unknown format {f}"),
    }
}

#[test]
fn writers_match_r() {
    for (i, case) in load_key("io.json", "write").iter().enumerate() {
        let ctx = format!("write case {i} ({} {})", case["format"], case["args"]);
        let motifs: Vec<Motif> = case["motifs"]
            .as_array()
            .unwrap()
            .iter()
            .map(motif_from_record)
            .collect();
        let got = write_case(case, &motifs);
        match case["text"].as_str() {
            None => assert!(got.is_err(), "{ctx}: R failed with {}", case["error"]),
            Some(want) => {
                let got = got.unwrap_or_else(|e| panic!("{ctx}: {e}"));
                if case["format"] == "homer" {
                    // Thresholds follow HOMER's own scoring rather than the
                    // R package's (see homer_thresholds()).
                    let strip = |t: &str| {
                        t.lines()
                            .map(|l| match l.strip_prefix('>') {
                                Some(h) => h.split('\t').take(2).collect::<Vec<_>>().join("\t"),
                                None => l.to_string(),
                            })
                            .collect::<Vec<_>>()
                    };
                    assert_eq!(strip(&got), strip(want), "{ctx}");
                    continue;
                }
                // The R package leaves trailing spaces in MEME files and writes
                // E= 0 for motifs without an E-value; the port does neither.
                let tidied;
                let want = if case["format"] == "meme" {
                    tidied = want
                        .lines()
                        .map(|l| l.trim_end().trim_end_matches(" E= 0"))
                        .collect::<Vec<_>>()
                        .join("\n")
                        + if want.ends_with('\n') { "\n" } else { "" };
                    tidied.as_str()
                } else {
                    want
                };
                if got != want {
                    for (k, (g, w)) in got.lines().zip(want.lines()).enumerate() {
                        assert_eq!(g, w, "{ctx}: line {}", k + 1);
                    }
                    assert_eq!(got, want, "{ctx}");
                }
            }
        }
    }
}

#[test]
fn written_files_read_back() {
    let motifs = jaspar::read(&reference_file("jaspar.txt")).unwrap();
    let again = jaspar::read(&jaspar::write(&motifs).unwrap()).unwrap();
    assert_eq!(again, motifs);

    let motifs = meme::read(&reference_file("meme_full.txt"), &Default::default()).unwrap();
    let text = meme::write(&motifs, &Default::default()).unwrap();
    let again = meme::read(&text, &Default::default()).unwrap();
    for (a, m) in again.iter().zip(&motifs) {
        assert_eq!(a.name, m.name);
        assert_close_vec(a.matrix.as_slice(), m.matrix.as_slice(), "meme round trip");
    }
}

/// HOMER scores a site as the sum of ln(max(p, 0.001) / 0.25) over the
/// probabilities written in the file; check the written thresholds on that
/// scale.
#[test]
fn homer_thresholds() {
    let motifs = jaspar::read(&reference_file("jaspar.txt")).unwrap();
    let read_back = |text: &str| -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        let mut max = 0.0;
        let mut thr = None;
        for l in text.lines() {
            if let Some(h) = l.strip_prefix('>') {
                if let Some(t) = thr.take() {
                    out.push((t, max));
                }
                thr = Some(h.split('\t').nth(2).unwrap().parse::<f64>().unwrap());
                max = 0.0;
            } else {
                max += l
                    .split('\t')
                    .map(|v| (v.parse::<f64>().unwrap().max(0.001) / 0.25).ln())
                    .fold(f64::NEG_INFINITY, f64::max);
            }
        }
        out.push((thr.unwrap(), max));
        out
    };
    let text = homer::write(&motifs, &homer::Threshold::Logodds(vec![0.8])).unwrap();
    for (t, max) in read_back(&text) {
        assert!((t - 0.8 * max).abs() < 1e-12, "{t} vs {}", 0.8 * max);
    }
    let text = homer::write(&motifs, &homer::Threshold::Pvalue(vec![1e-4])).unwrap();
    for (t, max) in read_back(&text) {
        assert!(t > 0.0 && t <= max, "{t} vs max {max}");
    }
    // A 6-mer with no ambiguity: only the consensus site has a P-value of
    // 1/4^6, so the cutoff must sit above the best site with one mismatch.
    let m = universalmotif::create::from_strings(&["ACGTAC"], None, &Default::default()).unwrap();
    let text = homer::write(&[m], &homer::Threshold::Pvalue(vec![0.25f64.powi(6)])).unwrap();
    let (t, max) = read_back(&text)[0];
    let one_mismatch = 5.0 * 4f64.ln() + (0.001f64 / 0.25).ln();
    assert!(t > one_mismatch && t <= max, "{t} vs {one_mismatch}, {max}");
}
