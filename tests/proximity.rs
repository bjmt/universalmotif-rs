//! Proximity to genomic anchors and background matching against the R
//! package.

mod common;

use common::*;
use serde_json::Value;
use universalmotif::background::{
    BinType, MatchOptions, gc_fraction, genome_windows, match_axis, match_bkg,
};
use universalmotif::proximity::{
    ProximityCount, ProximityHit, ProximityMethod, ProximityOptions, motif_proximity,
};
use universalmotif::ranges::{GRange, GStrand, Orientation};

fn reference() -> Value {
    let path = format!(
        "{}/tests/reference/proximity.json",
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

fn opt(x: Option<f64>) -> f64 {
    x.unwrap_or(f64::NAN)
}

/// Ranges from R, moved from 1-based inclusive to 0-based half-open positions.
fn ranges(t: &Value, chroms: &[String]) -> Vec<GRange> {
    let names = texts(&t["seqnames"]);
    let (s, e) = (nums(&t["start"]), nums(&t["end"]));
    let strand = texts(&t["strand"]);
    (0..names.len())
        .map(|i| GRange {
            chrom: chroms.iter().position(|c| *c == names[i]).unwrap(),
            start: s[i] as i64 - 1,
            end: e[i] as i64,
            strand: strand[i].parse::<GStrand>().unwrap(),
        })
        .collect()
}

#[test]
fn proximity_matches_r() {
    let r = reference();
    let chroms = texts(&r["chromosome_names"]);
    let seqs = texts(&r["chromosomes"]);
    let lengths: Vec<i64> = seqs.iter().map(|s| s.len() as i64).collect();
    let hit_ranges = ranges(&r["hits"], &chroms);
    let motif_names = texts(&r["hits"]["motif"]);
    let mut sorted = motif_names.clone();
    sorted.sort();
    sorted.dedup();
    let scores = nums(&r["hits"]["score"]);
    let hits: Vec<ProximityHit> = (0..hit_ranges.len())
        .map(|i| ProximityHit {
            motif: sorted.iter().position(|m| *m == motif_names[i]).unwrap(),
            range: hit_ranges[i],
            score: scores[i],
        })
        .collect();
    let anchors = ranges(&r["anchors"], &chroms);
    for case in r["proximity"].as_array().unwrap() {
        let a = args(&case["args"]);
        let ctx = format!("proximity {a:?} universe={}", case["universe"]);
        let (universe, chrom_lengths): (Vec<GRange>, Vec<Option<i64>>) =
            if case["universe"].as_bool().unwrap() {
                (
                    ranges(&r["universe"], &chroms),
                    lengths.iter().map(|&l| Some(l)).collect(),
                )
            } else {
                (
                    lengths
                        .iter()
                        .enumerate()
                        .map(|(c, &l)| GRange::new(c, 0, l))
                        .collect(),
                    lengths.iter().map(|&l| Some(l)).collect(),
                )
            };
        let s = |k: &str, d: &str| a.get(k).and_then(Value::as_str).unwrap_or(d).to_string();
        let opts = ProximityOptions {
            method: s("method", "binomial").parse::<ProximityMethod>().unwrap(),
            count: s("count", "hits").parse::<ProximityCount>().unwrap(),
            orientation: s("orientation", "undirected")
                .parse::<Orientation>()
                .unwrap(),
            dist_thresholds: a
                .get("dist.thresholds")
                .map(|v| nums(v).iter().map(|&x| x as i64).collect()),
            ignore_strand: a
                .get("ignore.strand")
                .and_then(Value::as_bool)
                .unwrap_or(true),
            ..Default::default()
        };
        let mut got = motif_proximity(&hits, &anchors, &universe, &chrom_lengths, &opts).unwrap();
        let t = &case["table"];
        // The R package breaks ties in P-value by its own motif index, the
        // order in which motifs first appear among the hits.
        let r_index: std::collections::HashMap<String, f64> = texts(&t["motif"])
            .into_iter()
            .zip(nums(&t["motif.i"]))
            .collect();
        let r_order = |m: usize| r_index.get(&sorted[m]).copied().unwrap_or(f64::INFINITY);
        got.sort_by(|x, y| {
            x.pvalue
                .total_cmp(&y.pvalue)
                .then(r_order(x.motif).total_cmp(&r_order(y.motif)))
        });
        let names: Vec<String> = got.iter().map(|p| sorted[p.motif].clone()).collect();
        assert_eq!(names, texts(&t["motif"]), "{ctx}: motif");
        let nh: Vec<f64> = got.iter().map(|p| p.nhits as f64).collect();
        assert_eq!(nh, nums(&t["nhits"]), "{ctx}: nhits");
        type Get = fn(&universalmotif::proximity::Proximity) -> Option<f64>;
        let cols: [(&str, Get); 10] = [
            ("best.dist", |p| p.best_dist.map(|x| x as f64)),
            ("hits.near", |p| p.hits_near.map(|x| x as f64)),
            ("hits.far", |p| p.hits_far.map(|x| x as f64)),
            ("frac.near", |p| p.frac_near),
            ("expected.near", |p| p.expected_near),
            ("enrichment", |p| p.enrichment),
            ("log2.enrichment", |p| p.log2_enrichment),
            ("pvalue", |p| Some(p.pvalue)),
            ("qvalue", |p| Some(p.qvalue)),
            ("ks.D", |p| p.ks_d),
        ];
        for (col, f) in cols {
            let g: Vec<f64> = got.iter().map(|p| opt(f(p))).collect();
            assert_close_vec(&g, &nums(&t[col]), &format!("{ctx}: {col}"));
        }
        for (p, w) in got.iter().zip(t["dist.to.nearest"].as_array().unwrap()) {
            let g: Vec<f64> = p.dist_to_nearest.iter().map(|&x| opt(x)).collect();
            assert_close_vec(&g, &nums(w), &format!("{ctx}: distances"));
        }
    }
}

#[test]
fn permutation_test_behaves() {
    let r = reference();
    let chroms = texts(&r["chromosome_names"]);
    let lengths: Vec<i64> = texts(&r["chromosomes"])
        .iter()
        .map(|s| s.len() as i64)
        .collect();
    let motif_names = texts(&r["hits"]["motif"]);
    let hit_ranges = ranges(&r["hits"], &chroms);
    let scores = nums(&r["hits"]["score"]);
    let hits: Vec<ProximityHit> = (0..hit_ranges.len())
        .map(|i| ProximityHit {
            motif: usize::from(motif_names[i] != "creb"),
            range: hit_ranges[i],
            score: scores[i],
        })
        .collect();
    let anchors = ranges(&r["anchors"], &chroms);
    let universe: Vec<GRange> = lengths
        .iter()
        .enumerate()
        .map(|(c, &l)| GRange::new(c, 0, l))
        .collect();
    let cl: Vec<Option<i64>> = lengths.iter().map(|&l| Some(l)).collect();
    let mut opts = ProximityOptions {
        method: ProximityMethod::Permutation,
        n_perm: 199,
        seed: 3,
        ..Default::default()
    };
    let a = motif_proximity(&hits, &anchors, &universe, &cl, &opts).unwrap();
    // The planted motif beats every permutation; the rest do not.
    assert!((a[0].pvalue - 1.0 / 200.0).abs() < 1e-12, "{}", a[0].pvalue);
    assert!(a[1].pvalue > 0.05, "{}", a[1].pvalue);
    opts.nthreads = 2;
    assert_eq!(
        a,
        motif_proximity(&hits, &anchors, &universe, &cl, &opts).unwrap()
    );
    opts.perm_per_chrom = false;
    opts.count = ProximityCount::Anchors;
    let b = motif_proximity(&hits, &anchors, &universe, &cl, &opts).unwrap();
    assert!(b[0].pvalue < 0.05);
}

#[test]
fn match_axes_match_r() {
    let r = reference();
    let t = texts(&r["match_targets"]);
    let u = texts(&r["match_universe"]);
    let gc_t: Vec<f64> = t.iter().map(|s| gc_fraction(s)).collect();
    let gc_u: Vec<f64> = u.iter().map(|s| gc_fraction(s)).collect();
    assert_close_vec(&gc_t, &nums(&r["gc_targets"]), "gc targets");
    assert_close_vec(&gc_u, &nums(&r["gc_universe"]), "gc universe");
    for (i, c) in r["axes"].as_array().unwrap().iter().enumerate() {
        let range = (!c["range"].is_null()).then(|| {
            let v = nums(&c["range"]);
            (v[0], v[1])
        });
        let bt = c["type"].as_str().unwrap().parse::<BinType>().unwrap();
        let ax = match_axis(
            &nums(&c["targets"]),
            &nums(&c["universe"]),
            num(&c["n_bins"]) as usize,
            bt,
            range,
        );
        let one = |v: &[usize]| v.iter().map(|&x| x as f64 + 1.0).collect::<Vec<_>>();
        assert_eq!(one(&ax.bin_targets), nums(&c["bin.t"]), "axis {i}: targets");
        assert_eq!(
            one(&ax.bin_universe),
            nums(&c["bin.u"]),
            "axis {i}: universe"
        );
        assert_eq!(ax.n_bins as f64, num(&c["n.bins"]), "axis {i}: bins");
    }
}

#[test]
fn matched_backgrounds_follow_composition() {
    let r = reference();
    let t = texts(&r["match_targets"]);
    let u = texts(&r["match_universe"]);
    let opts = MatchOptions {
        n_per_target: 3,
        seed: 5,
        ..Default::default()
    };
    let m = match_bkg(&t, &u, &[], &opts).unwrap();
    assert_eq!(m.universe_index.len(), 36);
    // Draws repeat only for targets whose nearby bins ran out.
    let uniq: std::collections::HashSet<_> = m.universe_index.iter().collect();
    assert!(uniq.len() == 36 || m.fallback_targets > 0);
    let coarse = MatchOptions {
        n_bins_gc: 4,
        n_bins_length: 3,
        ..opts.clone()
    };
    let c = match_bkg(&t, &u, &[], &coarse).unwrap();
    assert_eq!(c.fallback_targets, 0);
    let uniq: std::collections::HashSet<_> = c.universe_index.iter().collect();
    assert_eq!(uniq.len(), 36, "unique draws");
    // Matched GC should track the targets far better than random draws.
    let mean_gap: f64 = m
        .universe_index
        .iter()
        .enumerate()
        .map(|(k, &i)| (m.gc_universe[i] - m.gc_targets[k / 3]).abs())
        .sum::<f64>()
        / 36.0;
    let random_gap: f64 = (0..36)
        .map(|k| (m.gc_universe[(k * 7) % u.len()] - m.gc_targets[k / 3]).abs())
        .sum::<f64>()
        / 36.0;
    assert!(mean_gap < random_gap, "{mean_gap} vs {random_gap}");
    assert_eq!(m, match_bkg(&t, &u, &[], &opts).unwrap());

    // A universe too small for unique draws falls back to repeats.
    let small = match_bkg(&t, &u[..10], &[], &opts).unwrap();
    assert!(small.fallback_targets > 0);
    assert_eq!(small.universe_index.len(), 36);

    let w = genome_windows(
        &[5000, 200, 3000],
        &[150],
        500,
        &[GRange::new(0, 0, 2500)],
        2,
    )
    .unwrap();
    assert_eq!(w.len(), 500);
    assert!(w.iter().all(|g| g.width() == 150));
    assert!(w.iter().all(|g| g.chrom != 0 || g.start >= 2500));
    assert!(w.iter().all(|g| g.end <= [5000, 200, 3000][g.chrom]));
}
