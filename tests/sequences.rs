//! k-mers, backgrounds, complexity and motif discovery against the R
//! package.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::complexity::{
    ComplexityMethod, ComplexityOptions, complexity, sequence_complexity,
};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::kmer::{BkgOptions, WindowOptions, count_kmers, get_bkg, kmers, windows};

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a.iter().map(|s| s.as_str().unwrap().to_string()).collect(),
        Value::String(s) => vec![s.clone()],
        _ => vec![],
    }
}

fn chars(s: &str) -> Vec<char> {
    s.chars().collect()
}

fn usizes(v: &Value) -> Vec<usize> {
    nums(v).iter().map(|&x| x as usize).collect()
}

#[test]
fn kmers_match_r() {
    for case in load_key("sequences.json", "klets") {
        let letters: Vec<char> = texts(&case["letters"])
            .iter()
            .map(|s| s.chars().next().unwrap())
            .collect();
        let k = num(&case["k"]) as usize;
        assert_eq!(kmers(&letters, k), texts(&case["klets"]));
    }
    for case in load_key("sequences.json", "count") {
        let s = case["string"].as_str().unwrap();
        let k = num(&case["k"]) as usize;
        let mut alph = match case["alph"].as_str() {
            Some(a) => chars(a),
            None => chars(s),
        };
        alph.sort_unstable();
        alph.dedup();
        assert_eq!(kmers(&alph, k), texts(&case["klets"]), "{s} k={k}");
        let got: Vec<f64> = count_kmers(s, k, &alph).iter().map(|&c| c as f64).collect();
        assert_eq!(got, nums(&case["counts"]), "{s} k={k}");
    }
}

#[test]
fn windows_match_r() {
    for case in load_key("sequences.json", "windows") {
        let (n, w, o) = (
            num(&case["n"]) as usize,
            num(&case["window"]) as usize,
            num(&case["overlap"]) as usize,
        );
        let incomplete = case["incomplete"].as_bool().unwrap();
        if !incomplete {
            // The R package also drops a complete last window, so these are
            // checked against the full set of windows instead.
            continue;
        }
        let got = windows(n, w, o, incomplete);
        // R counts positions from 1, with the end included.
        let start: Vec<usize> = got.iter().map(|x| x.0 + 1).collect();
        let stop: Vec<usize> = got.iter().map(|x| x.1).collect();
        let ctx = format!("n={n} window={w} overlap={o} incomplete={incomplete}");
        assert_eq!(start, usizes(&case["start"]), "{ctx}: start");
        assert_eq!(stop, usizes(&case["stop"]), "{ctx}: stop");
        let complete: Vec<(usize, usize)> =
            got.iter().copied().filter(|&(a, b)| b - a == w).collect();
        assert_eq!(windows(n, w, o, false), complete, "{ctx}: complete only");
    }
}

#[test]
fn complexity_matches_r() {
    for case in load_key("sequences.json", "complexity") {
        let s = case["string"].as_str().unwrap();
        let method = case["method"]
            .as_str()
            .unwrap()
            .parse::<ComplexityMethod>()
            .unwrap();
        let alph = case["alph"].as_str();
        let max_word = num(&case["max_word"]) as usize;
        let got = complexity(s, method, alph, max_word);
        let want = num(&case["value"]);
        assert!(
            close(got, want),
            "{s} {method:?} {alph:?} {max_word}: {got} vs {want}"
        );
    }
}

#[test]
fn sequence_complexity_matches_r() {
    let seqs = texts(&load_key("sequences.json", "complexity_sequences").into());
    let names = texts(&load_key("sequences.json", "complexity_names").into());
    for case in load_key("sequences.json", "sequence_complexity") {
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
        let stop: Vec<usize> = got.iter().map(|r| r.end).collect();
        assert_eq!(start, usizes(&t["start"]), "{ctx}: start");
        assert_eq!(stop, usizes(&t["stop"]), "{ctx}: stop");
        let cx: Vec<f64> = got.iter().map(|r| r.complexity).collect();
        assert_close_vec(&cx, &nums(&t["complexity"]), &format!("{ctx}: complexity"));
    }
}

#[test]
fn backgrounds_match_r() {
    let seqs = texts(&load_key("sequences.json", "bkg_sequences").into());
    let names = texts(&load_key("sequences.json", "bkg_names").into());
    for case in load_key("sequences.json", "bkg") {
        let a = args(&case["args"]);
        let ctx = format!("get_bkg {a:?}");
        let mut opts = BkgOptions::default();
        if let Some(k) = a.get("k") {
            opts.k = usizes(k);
        }
        if let Some(p) = a.get("pseudocount") {
            opts.pseudocount = num(p) as u64;
        }
        if let Some(al) = a.get("alphabet") {
            opts.alphabet = Some(chars(al.as_str().unwrap()));
        }
        opts.rc = a.get("RC").and_then(Value::as_bool).unwrap_or(false);
        opts.merge = a.get("merge.res").and_then(Value::as_bool).unwrap_or(true);
        if a.get("window").and_then(Value::as_bool) == Some(true) {
            opts.window = Some(WindowOptions {
                size: a.get("window.size").map(num).unwrap_or(0.1),
                overlap: a.get("window.overlap").map(num).unwrap_or(0.0),
            });
        }
        let t = &case["table"];
        let got = get_bkg(&seqs, &names, &opts).unwrap();
        let kmer: Vec<String> = got.iter().map(|r| r.kmer.clone()).collect();
        assert_eq!(kmer, texts(&t["klet"]), "{ctx}: kmer");
        let count: Vec<f64> = got.iter().map(|r| r.count as f64).collect();
        assert_eq!(count, nums(&t["count"]), "{ctx}: count");
        let prob: Vec<f64> = got.iter().map(|r| r.probability).collect();
        assert_close_vec(&prob, &nums(&t["probability"]), &format!("{ctx}: prob"));
        if !t["sequence"].is_null() {
            let n = names.len();
            let sq: Vec<String> = got
                .iter()
                .map(|r| names[r.sequence.unwrap() % n].clone())
                .collect();
            assert_eq!(sq, texts(&t["sequence"]), "{ctx}: sequence");
        }
        if !t["start"].is_null() {
            let start: Vec<usize> = got.iter().map(|r| r.start.unwrap() + 1).collect();
            let stop: Vec<usize> = got.iter().map(|r| r.end.unwrap()).collect();
            assert_eq!(start, usizes(&t["start"]), "{ctx}: start");
            assert_eq!(stop, usizes(&t["stop"]), "{ctx}: stop");
        }
    }
}

#[test]
fn motif_finder_matches_r() {
    let pos = texts(&load_key("sequences.json", "finder_positives").into());
    let neg = texts(&load_key("sequences.json", "finder_negatives").into());
    for case in load_key("sequences.json", "finder") {
        let a = args(&case["args"]);
        let ctx = format!("motif_finder {a:?}");
        let mut o = FinderOptions::default();
        let int = |k: &str, d: usize| a.get(k).map(|v| num(v) as usize).unwrap_or(d);
        o.min_width = int("min.width", o.min_width);
        o.max_width = int("max.width", o.max_width);
        o.nmotifs = int("nmotifs", o.nmotifs);
        o.pseudocount = int("pseudocount", o.pseudocount as usize) as u32;
        o.hit_pvalue = a.get("hit.pvalue").map(num).unwrap_or(o.hit_pvalue);
        o.stop_pvalue = a.get("stop.pvalue").map(num).unwrap_or(o.stop_pvalue);
        // The R package's default cutoff, which the reference values used.
        o.qvalue = a.get("qvalue").map(num).unwrap_or(0.05);
        o.dedup_overlap = a.get("dedup.overlap").map(num).unwrap_or(o.dedup_overlap);
        o.rc = a.get("RC").and_then(Value::as_bool).unwrap_or(true);
        // Both correct P-values for the words searched and count every word
        // in the q-values. Every motif R reports is looked for among all of
        // the port's, whose P-value and q-value must be R's.
        let kept = motif_finder(&pos, Some(&neg), &o).unwrap();
        o.qvalue = 1.0;
        let got = motif_finder(&pos, Some(&neg), &o).unwrap();
        let want = case["motifs"].as_array().unwrap();
        assert!(kept.len() <= want.len(), "{ctx}: number of motifs");
        let s = &case["stats"];
        for (i, w) in want.iter().enumerate() {
            let ctx = format!("{ctx}, motif {i}");
            let wm = motif_from_record(w);
            let g = got
                .iter()
                .find(|g| g.motif.matrix == wm.matrix)
                .unwrap_or_else(|| panic!("{ctx}: not found"));
            assert_eq!(g.motif.nsites, wm.nsites, "{ctx}: nsites");
            assert_eq!(g.motif.bkg, wm.bkg, "{ctx}: bkg");
            let col = |k: &str| nums(&s[k])[i];
            assert_eq!(g.seqs_pos as f64, col("seqs_pos"), "{ctx}: seqs_pos");
            assert_eq!(g.seqs_neg as f64, col("seqs_neg"), "{ctx}: seqs_neg");
            assert_eq!(g.sites_pos as f64, col("sites_pos"), "{ctx}: sites_pos");
            assert_eq!(g.sites_neg as f64, col("sites_neg"), "{ctx}: sites_neg");
            assert_eq!(g.n_pos as f64, col("n_pos"), "{ctx}: n_pos");
            assert_eq!(g.n_neg as f64, col("n_neg"), "{ctx}: n_neg");
            let (want_p, want_q) = (col("pvalue"), col("qvalue"));
            assert!(
                (g.pvalue - want_p).abs() <= 1e-9 * want_p,
                "{ctx}: pvalue {} vs {want_p}",
                g.pvalue
            );
            assert!(
                (g.qvalue - want_q).abs() <= 1e-9 * want_q,
                "{ctx}: qvalue {} vs {want_q}",
                g.qvalue
            );
            assert!(
                (g.log10_pvalue - g.pvalue.log10()).abs() < 1e-9,
                "{ctx}: log10"
            );
        }
        for g in &kept {
            assert!(
                g.qvalue <= a.get("qvalue").map(num).unwrap_or(0.05),
                "{ctx}"
            );
        }
    }
}
