//! Edge cases of enrichment, positional and genomic proximity testing,
//! co-occurrence, background matching, motif discovery, ranges and genomes
//! on disk. Counts come from small hand-made inputs and P-values are checked
//! against binomial and hypergeometric sums computed here from scratch.

use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::DeflateEncoder;
use universalmotif::create::{CreateOptions, from_strings};
use universalmotif::{Alphabet, Error, Motif};

use universalmotif::background::{
    BinType, Covariate, MatchOptions, gc_fraction, genome_windows, match_axis, match_bkg,
};
use universalmotif::coocc::{CooccHit, CooccOptions, motif_coocc};
use universalmotif::enrich::{EnrichOptions, EnrichTest, enrich_motifs};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::genome::Genome;
use universalmotif::peaks::{PeakHit, PeaksMode, PeaksOptions, motif_peaks};
use universalmotif::proximity::{
    ProximityCount, ProximityHit, ProximityMethod, ProximityOptions, motif_proximity,
};
use universalmotif::ranges::{
    GRange, GStrand, Orientation, distance, intersect, overlaps_any, reduce, total_width, zone,
};

// Independent reference calculations.

fn ln_fact(n: u64) -> f64 {
    (2..=n).map(|i| (i as f64).ln()).sum()
}

fn ln_choose(n: u64, k: u64) -> f64 {
    ln_fact(n) - ln_fact(k) - ln_fact(n - k)
}

/// P(X >= a) for the top-left cell of `[[a, b], [c, d]]` with fixed margins.
fn fisher_upper(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let (r1, r2, c1) = (a + b, c + d, a + c);
    let n = r1 + r2;
    (a..=r1.min(c1))
        .filter(|&x| c1 - x <= r2)
        .map(|x| (ln_choose(r1, x) + ln_choose(r2, c1 - x) - ln_choose(n, c1)).exp())
        .sum()
}

/// P(X >= k) for X binomial with `n` trials and success probability `p`.
fn binom_upper(k: u64, n: u64, p: f64) -> f64 {
    (k..=n)
        .map(|i| (ln_choose(n, i) + i as f64 * p.ln() + (n - i) as f64 * (1.0 - p).ln()).exp())
        .sum::<f64>()
        .min(1.0)
}

/// Benjamini-Hochberg q-values, written plainly as the minimum over larger
/// P-values of `m p / rank`.
fn bh(p: &[f64]) -> Vec<f64> {
    let m = p.len() as f64;
    p.iter()
        .map(|&pi| {
            p.iter()
                .filter(|&&pj| pj >= pi)
                .map(|&pj| {
                    let rank = p.iter().filter(|&&x| x <= pj).count() as f64;
                    m * pj / rank
                })
                .fold(1.0f64, f64::min)
        })
        .collect()
}

fn near(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-300)
}

#[track_caller]
fn assert_near(got: f64, want: f64, ctx: &str) {
    assert!(near(got, want, 1e-9), "{ctx}: got {got}, want {want}");
}

#[track_caller]
fn arg_msg<T: std::fmt::Debug>(r: universalmotif::Result<T>) -> String {
    match r {
        Err(Error::Argument(m)) => m,
        other => panic!("expected an argument error, got {other:?}"),
    }
}

fn motif(site: &str) -> Motif {
    let o = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    from_strings(&[site], None, &o).unwrap()
}

fn tempdir() -> PathBuf {
    static COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "universalmotif-genomic-edge-{}-{n}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

// enrich_motifs

/// Targets and background of A and T only, with G and GC runs planted, so
/// that only exact copies of the planted words are hits.
fn enrich_inputs() -> (Vec<Motif>, Vec<&'static str>, Vec<&'static str>) {
    let motifs = vec![motif("GGGGGGGG"), motif("CACGTGCA"), motif("GCGCGCGC")];
    let targets = vec![
        "AAAAGGGGGGGGTTTTTTTT",
        "ATATATGGGGGGGGGATATA",
        "TTGGGGGGGGTTTTAAAAAA",
        "AAAAAAGCGCGCGCAAAAAA",
    ];
    let bkg = vec![
        "AGGGGGGGGAAAAATTTTTT",
        "ATATATATATATATATATAT",
        "TTTTTTTTTTAAAAAAAAAA",
        "AAAATTTTAAAATTTTAAAA",
        "TATATATATATATATAAAAA",
    ];
    (motifs, targets, bkg)
}

#[test]
fn enrich_counts_sequences_with_hits() {
    let (motifs, targets, bkg) = enrich_inputs();
    let opts = EnrichOptions {
        rc: false,
        ..Default::default()
    };
    let res = enrich_motifs(&motifs, &targets, Some(&bkg), &opts).unwrap();
    assert_eq!(res.len(), 3);
    let table: Vec<_> = res
        .iter()
        .map(|e| {
            (
                e.motif,
                e.target_seq_n,
                e.target_seq_hits,
                e.target_site_hits,
                e.bkg_seq_n,
                e.bkg_seq_hits,
                e.bkg_site_hits,
            )
        })
        .collect();
    // The run of nine Gs holds two overlapping hits.
    assert_eq!(
        table,
        vec![
            (0, 4, 3, 4, 5, 1, 1),
            (1, 4, 0, 0, 5, 0, 0),
            (2, 4, 1, 1, 5, 0, 0)
        ]
    );
    // Hit rates of 3/4 against 1/5; with no background hits, the background
    // rate is floored at half a sequence.
    assert_near(res[0].enrichment, 3.75, "enrichment");
    assert_near(res[0].log2_enrichment, 3.75f64.log2(), "log2");
    assert_eq!(res[1].enrichment, 0.0);
    assert_eq!(res[1].log2_enrichment, f64::NEG_INFINITY);
    assert_near(res[2].enrichment, 0.25 / (0.5 / 5.0), "floored enrichment");
    // A pseudocount of 1 goes into every cell.
    let want = [
        fisher_upper(4, 2, 2, 5),
        fisher_upper(1, 1, 5, 6),
        fisher_upper(2, 1, 4, 6),
    ];
    for (e, w) in res.iter().zip(want) {
        assert_near(e.pvalue, w, "pvalue");
        assert_near(e.log10_pvalue, w.log10(), "log10 pvalue");
    }
    let q = bh(&want);
    for (e, w) in res.iter().zip(&q) {
        assert_near(e.qvalue, *w, "qvalue");
        assert_eq!(e.log10_qvalue, e.qvalue.log10());
    }
    // Without the pseudocount, the table is used as it is.
    let raw = enrich_motifs(
        &motifs,
        &targets,
        Some(&bkg),
        &EnrichOptions {
            pseudocount: 0,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_near(raw[0].pvalue, fisher_upper(3, 1, 1, 4), "no pseudocount");
    assert_near(raw[1].pvalue, 1.0, "no hits anywhere");
    assert_eq!(raw[0].enrichment, res[0].enrichment);
}

#[test]
fn enrich_counts_sites_per_position_and_strand() {
    let (motifs, targets, bkg) = enrich_inputs();
    // Twenty letters hold 13 positions for a motif of width 8.
    for (rc, strands) in [(false, 1u64), (true, 2)] {
        let opts = EnrichOptions {
            rc,
            test: EnrichTest::Sites,
            ..Default::default()
        };
        let res = enrich_motifs(&motifs, &targets, Some(&bkg), &opts).unwrap();
        let (ppos, npos) = (4 * 13 * strands, 5 * 13 * strands);
        // GCGCGCGC is its own reverse complement, so with both strands each
        // copy is two hits.
        let gc_hits = strands;
        assert_eq!(res[0].target_site_hits, 4, "rc {rc}");
        assert_eq!(res[0].bkg_site_hits, 1, "rc {rc}");
        assert_eq!(res[2].target_site_hits as u64, gc_hits, "rc {rc}");
        assert_eq!(res[2].target_seq_hits, 1, "rc {rc}");
        assert_near(
            res[0].enrichment,
            (4.0 / ppos as f64) / (1.0 / npos as f64),
            "site enrichment",
        );
        assert_near(
            res[2].enrichment,
            (gc_hits as f64 / ppos as f64) / (0.5 / npos as f64),
            "floored site enrichment",
        );
        assert_eq!(res[1].enrichment, 0.0);
        assert_near(
            res[0].pvalue,
            fisher_upper(4 + 1, 1 + 1, ppos - 4 + 1, npos - 1 + 1),
            "site pvalue",
        );
        assert_near(
            res[2].pvalue,
            fisher_upper(gc_hits + 1, 1, ppos - gc_hits + 1, npos + 1),
            "site pvalue without background hits",
        );
    }
}

#[test]
fn enrich_sites_with_sequences_shorter_than_the_motif() {
    let m = vec![motif("GGGGGGGG")];
    let res = enrich_motifs(
        &m,
        &["GGGG", "AAA"],
        Some(&["GGGGG"][..]),
        &EnrichOptions {
            test: EnrichTest::Sites,
            ..Default::default()
        },
    )
    .unwrap();
    // No position fits, so every count is zero and the table is pure
    // pseudocount.
    assert_eq!(res[0].target_site_hits, 0);
    assert_eq!(res[0].enrichment, 0.0);
    assert_eq!(res[0].log2_enrichment, f64::NEG_INFINITY);
    assert_near(res[0].pvalue, fisher_upper(1, 1, 1, 1), "pvalue");
}

#[test]
fn enrich_protein_motifs_scan_one_strand() {
    let o = CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let m = vec![from_strings(&["WWWWWW"], Some(Alphabet::Aa), &o).unwrap()];
    let targets = ["AAAAWWWWWWAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA"];
    let bkg = ["AAAAAAAAAAAAAAAAAAAA"; 3];
    let res = enrich_motifs(
        &m,
        &targets,
        Some(&bkg[..]),
        &EnrichOptions {
            test: EnrichTest::Sites,
            pvalue: 1e-7,
            rc: true,
            ..Default::default()
        },
    )
    .unwrap();
    // Only nucleotides have a reverse strand, so positions are not doubled.
    let (ppos, npos) = (2 * 15, 3 * 15);
    assert_eq!(res[0].target_site_hits, 1);
    assert_eq!(res[0].bkg_site_hits, 0);
    assert_near(
        res[0].enrichment,
        (1.0 / ppos as f64) / (0.5 / npos as f64),
        "enrichment",
    );
    assert_near(
        res[0].pvalue,
        fisher_upper(2, 1, ppos as u64, npos as u64 + 1),
        "pvalue",
    );
}

#[test]
fn enrich_shuffled_background_is_seeded() {
    let (motifs, targets, _) = enrich_inputs();
    let opts = EnrichOptions {
        seed: 9,
        shuffle_k: 1,
        ..Default::default()
    };
    let a = enrich_motifs(&motifs, &targets, None, &opts).unwrap();
    assert!(a.iter().all(|e| e.bkg_seq_n == targets.len()));
    let b = enrich_motifs(
        &motifs,
        &targets,
        None,
        &EnrichOptions {
            nthreads: 2,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_eq!(a, b);
    // Target counts do not depend on the background.
    let c = enrich_motifs(&motifs, &targets, None, &EnrichOptions { seed: 10, ..opts }).unwrap();
    for (x, y) in a.iter().zip(&c) {
        assert_eq!(x.target_seq_hits, y.target_seq_hits);
        assert_eq!(x.target_site_hits, y.target_site_hits);
    }
}

#[test]
fn enrich_rejects_bad_input() {
    let (motifs, targets, bkg) = enrich_inputs();
    let o = EnrichOptions::default();
    let none: [&str; 0] = [];
    assert_eq!(
        arg_msg(enrich_motifs(&[], &targets, Some(&bkg), &o)),
        "no motifs given"
    );
    let rna = from_strings(
        &["GGGGGGGG"],
        Some(Alphabet::Rna),
        &CreateOptions::default(),
    )
    .unwrap();
    assert_eq!(
        arg_msg(enrich_motifs(
            &[motifs[0].clone(), rna],
            &targets,
            Some(&bkg),
            &o
        )),
        "all motifs must share the same alphabet"
    );
    for p in [0.0, 1.0, -1.0, f64::NAN] {
        assert_eq!(
            arg_msg(enrich_motifs(
                &motifs,
                &targets,
                Some(&bkg),
                &EnrichOptions {
                    pvalue: p,
                    ..o.clone()
                }
            )),
            "pvalue must be in (0, 1)"
        );
    }
    assert_eq!(
        arg_msg(enrich_motifs(&motifs, &none, Some(&bkg), &o)),
        "no sequences given"
    );
    assert_eq!(
        arg_msg(enrich_motifs(&motifs, &targets, Some(&none), &o)),
        "no background sequences"
    );
    // A bad shuffle is only an error when the background is shuffled.
    let bad_k = EnrichOptions {
        shuffle_k: 0,
        ..o.clone()
    };
    assert_eq!(
        arg_msg(enrich_motifs(&motifs, &targets, None, &bad_k)),
        "k must be at least 1"
    );
    assert!(enrich_motifs(&motifs, &targets, Some(&bkg), &bad_k).is_ok());
}

#[test]
fn enum_names_round_trip() {
    for t in [EnrichTest::Seqs, EnrichTest::Sites] {
        assert_eq!(t.to_string().parse::<EnrichTest>().unwrap(), t);
    }
    assert_eq!(
        arg_msg("site".parse::<EnrichTest>()),
        "unknown test 'site'; expected 'seqs' or 'sites'"
    );
    for m in [PeaksMode::Central, PeaksMode::Local] {
        assert_eq!(m.to_string().parse::<PeaksMode>().unwrap(), m);
    }
    assert_eq!(
        arg_msg("middle".parse::<PeaksMode>()),
        "unknown mode 'middle'; expected 'central' or 'local'"
    );
    for m in [
        ProximityMethod::Binomial,
        ProximityMethod::Permutation,
        ProximityMethod::Ks,
    ] {
        assert_eq!(m.to_string().parse::<ProximityMethod>().unwrap(), m);
    }
    assert_eq!(
        arg_msg("fisher".parse::<ProximityMethod>()),
        "unknown method 'fisher'; expected 'binomial', 'permutation' or 'ks'"
    );
    for c in [ProximityCount::Hits, ProximityCount::Anchors] {
        assert_eq!(c.to_string().parse::<ProximityCount>().unwrap(), c);
    }
    assert_eq!(
        arg_msg("sites".parse::<ProximityCount>()),
        "unknown count 'sites'; expected 'hits' or 'anchors'"
    );
    for o in [
        Orientation::Undirected,
        Orientation::Upstream,
        Orientation::Downstream,
    ] {
        assert_eq!(o.to_string().parse::<Orientation>().unwrap(), o);
    }
    assert_eq!(
        arg_msg("both".parse::<Orientation>()),
        "unknown orientation 'both'; expected 'undirected', 'upstream' or 'downstream'"
    );
    for s in [GStrand::Plus, GStrand::Minus, GStrand::Any] {
        assert_eq!(s.to_string().parse::<GStrand>().unwrap(), s);
    }
    assert_eq!("".parse::<GStrand>().unwrap(), GStrand::Any);
    assert_eq!("+-".parse::<GStrand>().unwrap(), GStrand::Any);
    assert_eq!(
        arg_msg(".".parse::<GStrand>()),
        "unknown strand '.'; expected '+', '-' or '*'"
    );
    for b in [BinType::Quantile, BinType::EqualWidth] {
        assert_eq!(b.to_string().parse::<BinType>().unwrap(), b);
    }
    assert_eq!(
        arg_msg("equal".parse::<BinType>()),
        "unknown bin type 'equal'; expected 'quantile' or 'equalwidth'"
    );
}

// motif_peaks

fn peak(motif: usize, sequence: usize, centre: f64, score: f64) -> PeakHit {
    PeakHit {
        motif,
        sequence,
        centre,
        score,
    }
}

#[test]
fn peaks_central_windows_match_the_binomial() {
    // Motif 7 has six of its ten best hits within five positions of the
    // middle, position 49; motif 2 has three hits spread out.
    let mut hits = Vec::new();
    for (s, c) in [47.0, 48.0, 48.5, 49.5, 51.0, 54.0, 5.0, 15.0, 85.0, 95.0]
        .into_iter()
        .enumerate()
    {
        hits.push(peak(7, s, c, 10.0));
    }
    // Weaker hits elsewhere in the same sequences are ignored, as is an
    // equally good one after the first.
    hits.push(peak(7, 0, 90.0, 5.0));
    hits.push(peak(7, 1, 2.0, 10.0));
    for (s, c) in [10.0, 50.0, 90.0].into_iter().enumerate() {
        hits.push(peak(2, s, c, 1.0));
    }
    // A better hit replaces an earlier one.
    hits.push(peak(2, 0, 11.0, 2.0));
    let opts = PeaksOptions {
        min_window: 10,
        max_window: Some(30),
        window_step: 10,
        ..Default::default()
    };
    let res = motif_peaks(&hits, 100, &opts).unwrap();
    assert_eq!(res.iter().map(|p| p.motif).collect::<Vec<_>>(), [2, 7]);
    let p7 = &res[1];
    assert_eq!(p7.nhits, 10);
    assert_eq!(p7.best_window, Some(10));
    assert_eq!(p7.best_centre, Some(49.0));
    assert_eq!((p7.hits_in, p7.hits_out), (6, 4));
    assert_eq!(p7.expected_in, Some(1.0));
    assert_eq!(p7.enrichment, Some(6.0));
    assert_eq!(p7.log2_enrichment, Some(6f64.log2()));
    // Three widths were tried, each with the same six hits inside.
    let p_best = binom_upper(6, 10, 0.1);
    assert!(p_best < binom_upper(6, 10, 0.2));
    assert_near(p7.pvalue, 3.0 * p_best, "pvalue");
    // Centres are rounded half to even as positions counted from 1.
    assert_eq!(p7.centres, [47, 48, 49, 49, 51, 54, 5, 15, 85, 95]);
    let p2 = &res[0];
    assert_eq!(p2.nhits, 3);
    // The better hit in sequence 0 came last, so its centre does too.
    assert_eq!(p2.centres, [50, 90, 11]);
    assert_eq!(p2.best_window, Some(10));
    assert_eq!(p2.hits_in, 1);
    assert_near(p2.pvalue, 3.0 * binom_upper(1, 3, 0.1), "spread pvalue");
    let q = bh(&[p2.pvalue, p7.pvalue]);
    assert_near(p2.qvalue, q[0], "qvalue");
    assert_near(p7.qvalue, q[1], "qvalue");
    let threaded = motif_peaks(
        &hits,
        100,
        &PeaksOptions {
            nthreads: 2,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_eq!(res, threaded);
}

#[test]
fn peaks_local_windows_slide_along_the_sequence() {
    let centres = [18.0, 19.0, 20.0, 21.0, 22.0, 60.0, 70.0, 80.0, 90.0, 99.0];
    let hits: Vec<PeakHit> = centres
        .iter()
        .enumerate()
        .map(|(s, &c)| peak(0, s, c, 1.0))
        .collect();
    let opts = PeaksOptions {
        mode: PeaksMode::Local,
        min_window: 10,
        max_window: Some(10),
        position_step: 5,
        ..Default::default()
    };
    let res = motif_peaks(&hits, 100, &opts).unwrap();
    let p = &res[0];
    // Windows are centred at 5, 10, ..., 95 counting from 1, which is 19
    // tests; the first to hold all five clustered hits is centred at 19.
    let n_tests = (5..=95).step_by(5).count();
    assert_eq!(n_tests, 19);
    assert_eq!(p.best_window, Some(10));
    assert_eq!(p.best_centre, Some(19.0));
    assert_eq!(p.hits_in, 5);
    assert_near(
        p.pvalue,
        (n_tests as f64 * binom_upper(5, 10, 0.1)).min(1.0),
        "pvalue",
    );
    // A window as wide as an odd-length sequence has no legal centre, so
    // nothing is tested.
    let none = motif_peaks(
        &hits[..3],
        11,
        &PeaksOptions {
            min_window: 11,
            max_window: Some(11),
            ..opts
        },
    )
    .unwrap();
    assert_eq!(none[0].best_window, None);
    assert_eq!(none[0].best_centre, None);
    assert_eq!(none[0].expected_in, None);
    assert_eq!(none[0].enrichment, None);
    assert_eq!(none[0].pvalue, 1.0);
    assert_eq!((none[0].hits_in, none[0].hits_out), (0, 3));
}

#[test]
fn peaks_edge_cases_and_errors() {
    let o = PeaksOptions::default();
    assert!(motif_peaks(&[], 100, &o).unwrap().is_empty());
    let hits = [peak(0, 0, 3.0, 1.0)];
    assert_eq!(
        arg_msg(motif_peaks(&hits, 1, &o)),
        "the sequence length must be at least 2"
    );
    for bad in [
        PeaksOptions {
            min_window: 0,
            ..o.clone()
        },
        PeaksOptions {
            window_step: 0,
            ..o.clone()
        },
        PeaksOptions {
            position_step: 0,
            ..o.clone()
        },
    ] {
        assert_eq!(
            arg_msg(motif_peaks(&hits, 100, &bad)),
            "min_window, window_step and position_step must be at least 1"
        );
    }
    // The default widest window is half the length; a larger one is capped
    // at the length.
    assert_eq!(
        arg_msg(motif_peaks(&hits, 18, &o)),
        "the widest window (9) is narrower than min_window (10)"
    );
    assert_eq!(
        arg_msg(motif_peaks(
            &hits,
            8,
            &PeaksOptions {
                max_window: Some(50),
                ..o.clone()
            }
        )),
        "the widest window (8) is narrower than min_window (10)"
    );
    // A window covering the whole sequence holds every hit with
    // probability 1, so it is never significant.
    let all = motif_peaks(
        &hits,
        10,
        &PeaksOptions {
            max_window: Some(10),
            ..o
        },
    )
    .unwrap();
    assert_eq!(all[0].pvalue, 1.0);
    assert_eq!(all[0].best_window, None);
}

// motif_coocc

fn coocc_hit(motif: usize, sequence: usize, start: i64) -> CooccHit {
    CooccHit {
        motif,
        sequence,
        start,
    }
}

/// Three motifs over 70 sequences, so that presence spans two 64-bit words.
fn coocc_hits() -> Vec<CooccHit> {
    let mut h = Vec::new();
    for s in (0..11).chain([65, 66]) {
        h.push(coocc_hit(0, s, 10));
    }
    // Motif 1 shares sequences 5 to 10 and 66, at various distances.
    for (s, start) in [
        (5, 12),
        (6, 100),
        (7, 13),
        (8, 9),
        (9, 40),
        (10, 10),
        (66, 11),
        (69, 0),
    ] {
        h.push(coocc_hit(1, s, start));
    }
    h.push(coocc_hit(2, 30, 0));
    h
}

#[test]
fn coocc_tables_and_fisher_tests() {
    let hits = coocc_hits();
    let opts = CooccOptions {
        max_distance: Some(3),
        ..Default::default()
    };
    let res = motif_coocc(&hits, 3, 70, &opts).unwrap();
    assert_eq!(res.len(), 3);
    // The tested pair first, then the two untested pairs in input order.
    let ab = &res[0];
    assert_eq!((ab.a, ab.b), (0, 1));
    assert_eq!((ab.both, ab.a_only, ab.b_only, ab.neither), (7, 6, 1, 56));
    assert_near(ab.odds_ratio.unwrap(), 7.0 * 56.0 / 6.0, "odds ratio");
    assert_near(ab.pvalue.unwrap(), fisher_upper(7, 6, 1, 56), "pvalue");
    assert_near(ab.qvalue.unwrap(), ab.pvalue.unwrap(), "single q-value");
    // Nearest distances 2, 90, 3, 1, 30, 0 and 1; five are within 3.
    assert_eq!(ab.both_clustered, Some(5));
    assert_eq!(ab.median_distance, Some(1.0));
    for p in &res[1..] {
        assert_eq!(p.both, 0);
        assert_eq!((p.odds_ratio, p.pvalue, p.qvalue), (None, None, None));
        assert_eq!(p.both_clustered, Some(0));
        assert_eq!(p.median_distance, None);
    }
    assert_eq!((res[1].a, res[1].b), (0, 2));
    assert_eq!((res[2].a, res[2].b), (1, 2));
    assert_eq!((res[1].a_only, res[1].b_only, res[1].neither), (13, 1, 56));

    let pc = motif_coocc(
        &hits,
        3,
        70,
        &CooccOptions {
            pseudocount: 1,
            min_coocc: 0,
            ..Default::default()
        },
    )
    .unwrap();
    // Every pair is tested now, and the pseudocount goes into every cell
    // of the test while the counts stay as they are.
    assert!(pc.iter().all(|p| p.pvalue.is_some()));
    let ab = pc.iter().find(|p| (p.a, p.b) == (0, 1)).unwrap();
    assert_eq!(ab.both, 7);
    assert_near(ab.odds_ratio.unwrap(), 8.0 * 57.0 / (7.0 * 2.0), "odds");
    assert_near(ab.pvalue.unwrap(), fisher_upper(8, 7, 2, 57), "pvalue");
    assert_eq!(ab.both_clustered, None);
    let ps: Vec<f64> = pc.iter().map(|p| p.pvalue.unwrap()).collect();
    let q = bh(&ps);
    for (p, w) in pc.iter().zip(&q) {
        assert_near(p.qvalue.unwrap(), *w, "qvalue");
    }
    assert!(pc.windows(2).all(|w| w[0].qvalue <= w[1].qvalue));
    let threaded = motif_coocc(
        &hits,
        3,
        70,
        &CooccOptions {
            pseudocount: 1,
            min_coocc: 0,
            nthreads: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(pc, threaded);
}

#[test]
fn coocc_self_pairs() {
    // Duplicate hits, as palindromes give on both strands, count once.
    let hits = [
        coocc_hit(0, 0, 5),
        coocc_hit(0, 0, 5),
        coocc_hit(0, 0, 20),
        coocc_hit(0, 0, 9),
        coocc_hit(0, 1, 3),
        coocc_hit(0, 1, 3),
        coocc_hit(0, 2, 50),
        coocc_hit(0, 2, 40),
    ];
    let opts = CooccOptions {
        self_pairs: true,
        max_distance: Some(5),
        ..Default::default()
    };
    let res = motif_coocc(&hits, 1, 4, &opts).unwrap();
    assert_eq!(res.len(), 1);
    let p = &res[0];
    assert_eq!((p.a, p.b), (0, 0));
    assert_eq!((p.both, p.a_only, p.b_only, p.neither), (2, 1, 0, 1));
    assert_eq!((p.odds_ratio, p.pvalue, p.qvalue), (None, None, None));
    // Nearest distances are 4 and 10.
    assert_eq!(p.both_clustered, Some(1));
    assert_eq!(p.median_distance, Some(4.0));
    let wide = motif_coocc(
        &hits,
        1,
        4,
        &CooccOptions {
            max_distance: Some(10),
            ..opts.clone()
        },
    )
    .unwrap();
    assert_eq!(wide[0].both_clustered, Some(2));
    assert_eq!(wide[0].median_distance, Some(7.0));
    // With a second motif, self pairs come after the tested pair, and
    // without a distance there is nothing to summarize.
    let mut two = hits.to_vec();
    two.extend([coocc_hit(1, 0, 7), coocc_hit(1, 2, 0)]);
    let res = motif_coocc(
        &two,
        2,
        4,
        &CooccOptions {
            max_distance: None,
            ..opts
        },
    )
    .unwrap();
    let pairs: Vec<_> = res.iter().map(|p| (p.a, p.b)).collect();
    assert_eq!(pairs, [(0, 1), (0, 0), (1, 1)]);
    assert!(res.iter().all(|p| p.both_clustered.is_none()));
    assert_eq!((res[2].both, res[2].a_only, res[2].neither), (0, 2, 2));
    assert_near(res[0].pvalue.unwrap(), fisher_upper(2, 1, 0, 1), "pvalue");
}

#[test]
fn coocc_rejects_bad_input() {
    let o = CooccOptions::default();
    assert_eq!(
        arg_msg(motif_coocc(&[], 1, 5, &o)),
        "need at least 2 motifs (or self pairs)"
    );
    assert_eq!(
        arg_msg(motif_coocc(&[coocc_hit(2, 0, 0)], 2, 5, &o)),
        "a hit refers to motif 2 and sequence 0, but there are 2 motifs and 5 sequences"
    );
    assert_eq!(
        arg_msg(motif_coocc(&[coocc_hit(0, 5, 0)], 2, 5, &o)),
        "a hit refers to motif 0 and sequence 5, but there are 2 motifs and 5 sequences"
    );
    // No hits at all is not an error: nothing is shared, so nothing is tested.
    let none = motif_coocc(&[], 2, 5, &o).unwrap();
    assert_eq!(none[0].neither, 5);
    assert_eq!(none[0].pvalue, None);
}

// motif_proximity

fn ph(motif: usize, chrom: usize, start: i64, end: i64, score: f64) -> ProximityHit {
    ProximityHit {
        motif,
        range: GRange::new(chrom, start, end),
        score,
    }
}

fn stranded(chrom: usize, start: i64, end: i64, strand: GStrand) -> GRange {
    GRange {
        chrom,
        start,
        end,
        strand,
    }
}

/// One anchor at position 500 of a 1000-letter chromosome, and a second
/// chromosome without anchors.
fn proximity_setup() -> (
    Vec<ProximityHit>,
    Vec<GRange>,
    Vec<GRange>,
    Vec<Option<i64>>,
) {
    let hits = vec![
        // Middles 496, 501, 506, 102 and 802.
        ph(0, 0, 494, 499, 5.0),
        ph(0, 0, 499, 504, 5.0),
        ph(0, 0, 504, 509, 5.0),
        ph(0, 0, 100, 105, 5.0),
        ph(0, 0, 800, 805, 5.0),
        // On a chromosome without anchors, so outside the universe used.
        ph(0, 1, 10, 15, 9.0),
    ];
    let anchors = vec![stranded(0, 500, 501, GStrand::Plus)];
    let universe = vec![GRange::new(0, 0, 1000), GRange::new(1, 0, 500)];
    (hits, anchors, universe, vec![Some(1000), Some(500)])
}

#[test]
fn proximity_binomial_counts_hits_near_anchors() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let opts = ProximityOptions {
        dist_thresholds: Some(vec![50, 5, 50]),
        ..Default::default()
    };
    let res = motif_proximity(&hits, &anchors, &universe, &lens, &opts).unwrap();
    assert_eq!(res.len(), 1);
    let r = &res[0];
    assert_eq!(r.nhits, 5);
    // Within 5 of the anchor is 495..506, 11 of the 1000 positions, holding
    // the middles 496 and 501; within 50 is 101 positions holding three.
    let p5 = binom_upper(2, 5, 0.011);
    let p50 = binom_upper(3, 5, 0.101);
    assert!(p5 < p50);
    assert_eq!(r.best_dist, Some(5));
    assert_eq!((r.hits_near, r.hits_far), (Some(2), Some(3)));
    assert_near(r.frac_near.unwrap(), 0.4, "fraction");
    assert_near(r.expected_near.unwrap(), 5.0 * 0.011, "expected");
    assert_near(r.enrichment.unwrap(), 2.0 / 0.055, "enrichment");
    assert_near(r.log2_enrichment.unwrap(), (2.0f64 / 0.055).log2(), "log2");
    // Duplicate thresholds are tested once.
    assert_near(r.pvalue, 2.0 * p5, "pvalue");
    assert_near(r.qvalue, r.pvalue, "qvalue");
    assert_eq!(r.ks_d, None);
    // Signed distances from each middle to the anchor, which is 0 inside it
    // and counts the anchor position itself otherwise.
    assert_eq!(
        r.dist_to_nearest,
        [Some(-4.0), Some(1.0), Some(6.0), Some(-398.0), Some(302.0)]
    );
    // On the minus strand, upstream and downstream swap.
    let minus = [stranded(0, 500, 501, GStrand::Minus)];
    let m = motif_proximity(&hits, &minus, &universe, &lens, &opts).unwrap();
    assert_eq!(
        m[0].dist_to_nearest,
        [Some(4.0), Some(-1.0), Some(-6.0), Some(398.0), Some(-302.0)]
    );
    assert_eq!(m[0].pvalue, r.pvalue);
}

#[test]
fn proximity_directional_zones() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let base = ProximityOptions {
        dist_thresholds: Some(vec![5]),
        ..Default::default()
    };
    let run = |anchors: &[GRange], orientation| {
        motif_proximity(
            &hits,
            anchors,
            &universe,
            &lens,
            &ProximityOptions {
                orientation,
                ..base.clone()
            },
        )
        .unwrap()
        .remove(0)
    };
    let minus = [stranded(0, 500, 501, GStrand::Minus)];
    // Upstream of a + anchor is 495..500 (middle 496); downstream is 501..506
    // (middle 501). On the - strand they swap.
    for (a, o, middle_near) in [
        (&anchors[..], Orientation::Upstream, 496),
        (&anchors[..], Orientation::Downstream, 501),
        (&minus[..], Orientation::Upstream, 501),
        (&minus[..], Orientation::Downstream, 496),
    ] {
        let r = run(a, o);
        assert_eq!(r.hits_near, Some(1), "{o} {middle_near}");
        assert_near(r.expected_near.unwrap(), 5.0 * 5.0 / 1000.0, "expected");
        assert_near(r.pvalue, binom_upper(1, 5, 0.005), "pvalue");
    }
    let u = run(&anchors, Orientation::Undirected);
    assert_eq!(u.hits_near, Some(2));
    // An unstranded anchor among stranded ones counts as +.
    let mixed = [anchors[0], GRange::new(0, 100, 101)];
    let r = run(&mixed, Orientation::Downstream);
    // Downstream of 100..101 is 101..106, holding the middle 102.
    assert_eq!(r.hits_near, Some(2));
}

#[test]
fn proximity_collapses_overlapping_hits() {
    let (_, anchors, universe, lens) = proximity_setup();
    let hits = vec![
        stranded_hit(0, 494, 499, GStrand::Plus, 2.0),
        stranded_hit(0, 496, 501, GStrand::Minus, 3.0),
        stranded_hit(0, 497, 502, GStrand::Plus, 1.0),
        stranded_hit(0, 800, 805, GStrand::Plus, 1.0),
    ];
    let opts = ProximityOptions {
        dist_thresholds: Some(vec![5]),
        ..Default::default()
    };
    // Ignoring strand, the first three overlap and only the best is kept.
    let r = motif_proximity(&hits, &anchors, &universe, &lens, &opts).unwrap();
    assert_eq!(r[0].nhits, 2);
    assert_eq!(r[0].dist_to_nearest, [Some(-2.0), Some(302.0)]);
    // By strand, the two + hits still overlap each other.
    let r = motif_proximity(
        &hits,
        &anchors,
        &universe,
        &lens,
        &ProximityOptions {
            ignore_strand: false,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_eq!(r[0].nhits, 3);
    assert_eq!(r[0].dist_to_nearest, [Some(-4.0), Some(-2.0), Some(302.0)]);
    // Counting anchors, nothing is collapsed.
    let r = motif_proximity(
        &hits,
        &anchors,
        &universe,
        &lens,
        &ProximityOptions {
            count: ProximityCount::Anchors,
            ..opts
        },
    )
    .unwrap();
    assert_eq!(r[0].nhits, 1);
    assert_eq!(r[0].dist_to_nearest.len(), 4);
}

fn stranded_hit(chrom: usize, start: i64, end: i64, strand: GStrand, score: f64) -> ProximityHit {
    ProximityHit {
        motif: 0,
        range: stranded(chrom, start, end, strand),
        score,
    }
}

#[test]
fn proximity_counts_anchors_with_clipped_zones() {
    let (hits, _, universe, lens) = proximity_setup();
    let anchors = vec![
        GRange::new(0, 0, 1),
        GRange::new(0, 500, 501),
        GRange::new(0, 900, 901),
    ];
    let opts = ProximityOptions {
        count: ProximityCount::Anchors,
        dist_thresholds: Some(vec![5]),
        ..Default::default()
    };
    let r = motif_proximity(&hits, &anchors, &universe, &lens, &opts)
        .unwrap()
        .remove(0);
    assert_eq!(r.nhits, 3);
    // Five hits in 1000 positions. The zone of the anchor at 0 is clipped to
    // 0..6; the others are 11 wide. Only the anchor at 500 has a hit nearby.
    let lambda = 5.0 / 1000.0;
    let f = |w: f64| 1.0 - (-lambda * w).exp();
    let clipped = (f(6.0) + f(11.0) + f(11.0)) / 3.0;
    assert_eq!(r.hits_near, Some(1));
    assert_eq!(r.hits_far, Some(2));
    assert_near(r.expected_near.unwrap(), 3.0 * clipped, "expected");
    assert_near(r.pvalue, binom_upper(1, 3, clipped), "pvalue");
    // Without chromosome lengths, the zone runs past the start.
    let unclipped = motif_proximity(&hits, &anchors, &universe, &[], &opts)
        .unwrap()
        .remove(0);
    assert_near(
        unclipped.expected_near.unwrap(),
        3.0 * f(11.0),
        "unclipped expected",
    );
}

#[test]
fn proximity_skips_distances_covering_everything() {
    // Hits at every position of a ten-letter universe: any zone holds all of
    // them, so no distance can be tested.
    let hits: Vec<ProximityHit> = (0..10).map(|i| ph(0, 0, i, i + 1, 1.0)).collect();
    let anchors = [GRange::new(0, 5, 6)];
    let universe = [GRange::new(0, 0, 10)];
    for count in [ProximityCount::Hits, ProximityCount::Anchors] {
        let r = motif_proximity(
            &hits,
            &anchors,
            &universe,
            &[],
            &ProximityOptions {
                count,
                dist_thresholds: Some(vec![1000]),
                ..Default::default()
            },
        )
        .unwrap()
        .remove(0);
        assert_eq!(r.best_dist, None, "{count}");
        assert_eq!(r.expected_near, None, "{count}");
        assert_eq!(r.pvalue, 1.0, "{count}");
    }
    // With the default thresholds, the grid stops at half the universe, so
    // only 10 is tried, and it covers everything as well.
    let r = motif_proximity(&hits, &anchors, &universe, &[], &Default::default()).unwrap();
    assert_eq!(r[0].best_dist, None);
}

#[test]
fn proximity_default_thresholds_follow_the_universe() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let r = motif_proximity(&hits, &anchors, &universe, &lens, &Default::default())
        .unwrap()
        .remove(0);
    // The default grid for a 1000-letter universe runs up to 500, but at 500
    // the zone covers the whole universe and is not tested. The best of the
    // other five is the smallest distance.
    assert_eq!(r.best_dist, Some(10));
    assert_eq!(r.hits_near, Some(3));
    let grid = [10.0, 25.0, 50.0, 100.0, 250.0];
    let best = grid
        .iter()
        .map(|&d| {
            let k = [496, 501, 506, 102, 802]
                .iter()
                .filter(|&&m| (m as f64 - 500.0).abs() <= d)
                .count() as u64;
            binom_upper(k, 5, (2.0 * d + 1.0) / 1000.0)
        })
        .fold(1.0, f64::min);
    assert_near(r.pvalue, (5.0 * best).min(1.0), "pvalue");
}

#[test]
fn proximity_ks_test() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let opts = ProximityOptions {
        method: ProximityMethod::Ks,
        ..Default::default()
    };
    let r = motif_proximity(&hits, &anchors, &universe, &lens, &opts)
        .unwrap()
        .remove(0);
    assert!(r.ks_d.unwrap() > 0.0);
    assert!(r.pvalue < 0.05, "{}", r.pvalue);
    assert_eq!(r.best_dist, None);
    // An anchor beyond the universe, further from the hit than the
    // universe is wide, puts the hit where the null distribution is 1, so
    // the statistic is 0 and the P-value 1.
    let far = [ph(0, 0, 10, 11, 1.0)];
    let r = motif_proximity(&far, &[GRange::new(0, 5000, 5001)], &universe, &[], &opts)
        .unwrap()
        .remove(0);
    assert_eq!(r.dist_to_nearest, [Some(-4990.0)]);
    assert_eq!(r.ks_d, Some(0.0));
    assert_eq!(r.pvalue, 1.0);
}

#[test]
fn proximity_permutation_null() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let opts = ProximityOptions {
        method: ProximityMethod::Permutation,
        n_perm: 49,
        dist_thresholds: Some(vec![5]),
        seed: 4,
        ..Default::default()
    };
    let a = motif_proximity(&hits, &anchors, &universe, &lens, &opts).unwrap();
    // The P-value counts the observed statistic among the permutations.
    let k = (a[0].pvalue * 50.0).round();
    assert!((1.0..=50.0).contains(&k));
    assert_near(a[0].pvalue, k / 50.0, "permutation pvalue");
    assert_eq!(a[0].best_dist, Some(5));
    assert_eq!(a[0].hits_near, Some(2));
    assert_eq!(
        a,
        motif_proximity(&hits, &anchors, &universe, &lens, &opts).unwrap()
    );
    // An anchor wider than every segment of the universe cannot be placed.
    let wide = [GRange::new(0, 0, 50)];
    let small = [GRange::new(0, 0, 30)];
    assert_eq!(
        arg_msg(motif_proximity(
            &[ph(0, 0, 10, 12, 1.0)],
            &wide,
            &small,
            &[],
            &opts
        )),
        "could not place an anchor of width 50 in the universe; the universe is too small"
    );
}

#[test]
fn proximity_nearest_anchor_among_nested_ones() {
    // A long anchor reaching past a short one nested in it is the nearest
    // from either side, and from inside the short one's span it still
    // counts, as the hit lies within it.
    let anchors = [
        GRange::new(0, 0, 100),
        GRange::new(0, 10, 11),
        GRange::new(0, 50, 60),
    ];
    let universe = [GRange::new(0, 0, 1000)];
    let hits = [
        ph(0, 0, 150, 151, 1.0),
        ph(0, 0, 80, 81, 1.0),
        ph(0, 0, 55, 56, 1.0),
        ph(0, 0, 100, 101, 1.0),
    ];
    let r = motif_proximity(
        &hits,
        &anchors,
        &universe,
        &[],
        &ProximityOptions {
            dist_thresholds: Some(vec![5]),
            ..Default::default()
        },
    )
    .unwrap();
    // Distances count the anchor's own position, so the position just past
    // the long anchor is 1 away and 150 is 51 away; inside an anchor it is 0.
    assert_eq!(
        r[0].dist_to_nearest,
        [Some(51.0), Some(0.0), Some(0.0), Some(1.0)]
    );
}

#[test]
fn proximity_rejects_bad_input() {
    let (hits, anchors, universe, lens) = proximity_setup();
    let o = ProximityOptions::default();
    let run = |a: &[GRange], u: &[GRange], o: &ProximityOptions| {
        arg_msg(motif_proximity(&hits, a, u, &lens, o))
    };
    assert_eq!(
        run(
            &anchors,
            &universe,
            &ProximityOptions {
                method: ProximityMethod::Ks,
                count: ProximityCount::Anchors,
                ..o.clone()
            }
        ),
        "the KS test is not available when counting anchors; use the permutation test"
    );
    assert_eq!(
        run(
            &anchors,
            &universe,
            &ProximityOptions {
                method: ProximityMethod::Ks,
                orientation: Orientation::Upstream,
                ..o.clone()
            }
        ),
        "the KS test is not available for directional zones"
    );
    assert_eq!(
        run(
            &anchors,
            &universe,
            &ProximityOptions {
                n_perm: 0,
                ..o.clone()
            }
        ),
        "n_perm must be at least 1"
    );
    assert_eq!(
        run(
            &[GRange::new(0, 500, 501)],
            &universe,
            &ProximityOptions {
                orientation: Orientation::Downstream,
                ..o.clone()
            }
        ),
        "directional zones need stranded anchors"
    );
    assert_eq!(
        run(&[GRange::new(2, 5, 6)], &universe, &o),
        "the anchors and the universe share no chromosomes"
    );
    assert_eq!(
        run(
            &anchors,
            &universe,
            &ProximityOptions {
                dist_thresholds: Some(vec![]),
                ..o.clone()
            }
        ),
        "give at least one distance threshold"
    );
    assert_eq!(
        run(
            &anchors,
            &universe,
            &ProximityOptions {
                dist_thresholds: Some(vec![10, 0]),
                ..o.clone()
            }
        ),
        "distance thresholds must be positive"
    );
    // Hits that all fall outside the universe give no results.
    assert!(
        motif_proximity(&hits[5..], &anchors, &universe, &lens, &o)
            .unwrap()
            .is_empty()
    );
    // A hit whose middle lies outside, though it overlaps, is dropped too.
    let edge = [ph(0, 0, 998, 1004, 1.0)];
    assert!(
        motif_proximity(&edge, &anchors, &[GRange::new(0, 0, 1000)], &lens, &o)
            .unwrap()
            .is_empty()
    );
}

// ranges

#[test]
fn range_operations() {
    let r = GRange::new;
    let empty = r(0, 5, 3);
    assert!(empty.is_empty());
    assert_eq!(empty.width(), 0);
    assert!(r(0, 5, 5).is_empty());
    assert!(!r(0, 5, 6).is_empty());
    // Gaps below min_gap merge: a gap of 4 needs min_gap 5. Strands go.
    let x = [
        stranded(0, 10, 20, GStrand::Plus),
        stranded(0, 24, 30, GStrand::Minus),
        r(0, 15, 16),
        r(1, 0, 3),
        empty,
    ];
    assert_eq!(reduce(&x, 4), [r(0, 10, 20), r(0, 24, 30), r(1, 0, 3)]);
    assert_eq!(reduce(&x, 5), [r(0, 10, 30), r(1, 0, 3)]);
    assert!(reduce(&[empty], 1).is_empty());
    assert_eq!(total_width(&reduce(&x, 5)), 23.0);
    // Intersections skip chromosomes missing from either side.
    let a = reduce(&[r(0, 0, 10), r(2, 5, 15), r(3, 0, 4)], 1);
    let b = reduce(&[r(1, 0, 100), r(2, 0, 7), r(2, 9, 12), r(3, 4, 8)], 1);
    assert_eq!(intersect(&a, &b), [r(2, 5, 7), r(2, 9, 12)]);
    assert!(intersect(&a, &[]).is_empty());
    // Overlap is half-open: touching is not overlapping.
    assert!(overlaps_any(&r(2, 6, 7), &a));
    assert!(!overlaps_any(&r(2, 15, 20), &a));
    assert!(!overlaps_any(&r(3, 4, 9), &a));
    assert!(!overlaps_any(&r(2, 8, 8), &a));
    assert_eq!(distance(&r(0, 10, 20), &r(0, 0, 4)), 6);
    assert_eq!(distance(&r(0, 0, 4), &r(0, 10, 20)), 6);
}

#[test]
fn zones_on_both_strands_and_clipped() {
    let plus = stranded(0, 10, 12, GStrand::Plus);
    let minus = stranded(0, 10, 12, GStrand::Minus);
    let any = GRange::new(0, 10, 12);
    let z = |a: &GRange, o, len| {
        let g = zone(a, 3, o, len);
        (g.start, g.end)
    };
    use Orientation::*;
    assert_eq!(z(&plus, Undirected, None), (7, 15));
    assert_eq!(z(&minus, Undirected, None), (7, 15));
    assert_eq!(z(&plus, Upstream, None), (7, 10));
    assert_eq!(z(&any, Upstream, None), (7, 10));
    assert_eq!(z(&minus, Upstream, None), (12, 15));
    assert_eq!(z(&plus, Downstream, None), (12, 15));
    assert_eq!(z(&minus, Downstream, None), (7, 10));
    // Clipped to the chromosome: partly outside is trimmed, wholly outside
    // becomes empty at the nearest end.
    assert_eq!(z(&plus, Undirected, Some(13)), (7, 13));
    assert_eq!(z(&plus, Downstream, Some(12)), (12, 12));
    assert_eq!(z(&plus, Downstream, Some(11)), (11, 11));
    let start = stranded(0, 1, 2, GStrand::Minus);
    assert_eq!(z(&start, Undirected, Some(100)), (0, 5));
    assert_eq!(z(&start, Downstream, Some(100)), (0, 1));
    let zero = stranded(0, 0, 1, GStrand::Minus);
    assert_eq!(z(&zero, Downstream, Some(100)), (0, 0));
    assert_eq!(zone(&plus, 3, Upstream, None).strand, GStrand::Any);
}

// Background matching

#[test]
fn gc_fraction_counts_either_case() {
    assert_eq!(gc_fraction(""), 0.0);
    assert_eq!(gc_fraction("AT"), 0.0);
    assert_eq!(gc_fraction("gGcC"), 1.0);
    assert_eq!(gc_fraction("ACgtNn"), 2.0 / 6.0);
}

#[test]
fn match_axis_edge_cases() {
    // Equal-width bins over 0..1: the top break belongs to the last bin and
    // values outside go to the end bins.
    let ax = match_axis(
        &[0.0, 0.5, 1.0, -0.3, 1.7],
        &[0.25, 0.75],
        2,
        BinType::EqualWidth,
        Some((0.0, 1.0)),
    );
    assert_eq!(ax.n_bins, 2);
    assert_eq!(ax.bin_targets, [0, 1, 1, 0, 1]);
    assert_eq!(ax.bin_universe, [0, 1]);
    // Constant values collapse to one bin, either way.
    for bt in [BinType::EqualWidth, BinType::Quantile] {
        let ax = match_axis(&[3.0, 3.0], &[3.0; 4], 5, bt, None);
        assert_eq!(ax.n_bins, 1, "{bt}");
        assert_eq!(ax.bin_targets, [0, 0]);
        assert_eq!(ax.bin_universe, [0; 4]);
    }
    // A constant universe makes one quantile bin even if targets vary.
    let ax = match_axis(&[1.0, 9.0], &[3.0; 4], 5, BinType::Quantile, None);
    assert_eq!(ax.n_bins, 1);
    assert_eq!(ax.bin_targets, [0, 0]);
    // No bins asked for is one bin.
    let ax = match_axis(&[1.0, 9.0], &[3.0, 4.0], 0, BinType::EqualWidth, None);
    assert_eq!(ax.n_bins, 1);
    // Quantile bins of the universe, with targets outside it in the end bins.
    let ax = match_axis(
        &[-10.0, 2.5, 100.0],
        &[1.0, 2.0, 3.0, 4.0],
        2,
        BinType::Quantile,
        None,
    );
    assert_eq!(ax.n_bins, 2);
    assert_eq!(ax.bin_universe, [0, 0, 1, 1]);
    assert_eq!(ax.bin_targets, [0, 1, 1]);
}

#[test]
fn match_bkg_covariates_pick_their_group() {
    // Every universe sequence has the same GC and length; a covariate
    // splits them in two.
    let universe = vec!["ACGTACGTAC"; 20];
    let targets = vec!["ACGTACGTAC"; 3];
    let cov = Covariate {
        name: "signal".into(),
        targets: vec![1.0, 0.0, 1.0],
        universe: (0..20).map(|i| if i < 10 { 0.0 } else { 1.0 }).collect(),
        n_bins: 2,
        bin_type: BinType::EqualWidth,
    };
    let opts = MatchOptions {
        n_per_target: 3,
        seed: 7,
        ..Default::default()
    };
    let m = match_bkg(&targets, &universe, std::slice::from_ref(&cov), &opts).unwrap();
    assert_eq!(m.universe_index.len(), 9);
    assert_eq!(m.fallback_targets, 0);
    for (k, &i) in m.universe_index.iter().enumerate() {
        let group = if k / 3 == 1 { 0..10 } else { 10..20 };
        assert!(group.contains(&i), "draw {k} is {i}");
    }
    let uniq: std::collections::HashSet<_> = m.universe_index.iter().collect();
    assert_eq!(uniq.len(), 9);
    assert_eq!(m.gc_targets, [0.5; 3]);
    assert_eq!(
        m,
        match_bkg(&targets, &universe, std::slice::from_ref(&cov), &opts).unwrap()
    );
    // When a group runs out, the neighbouring group is searched before any
    // draw repeats.
    let greedy = MatchOptions {
        n_per_target: 6,
        ..opts.clone()
    };
    let g = match_bkg(&targets, &universe, std::slice::from_ref(&cov), &greedy).unwrap();
    let uniq: std::collections::HashSet<_> = g.universe_index.iter().collect();
    assert_eq!(uniq.len(), 18);
    assert_eq!(g.fallback_targets, 0);
    assert!(g.universe_index[12..].iter().any(|&i| i < 10));
}

#[test]
fn match_bkg_fallbacks() {
    // No universe sequence is anywhere near a GC-rich target, so its draws
    // come from the whole universe, with replacement.
    let universe = ["AAAAAAAAAA", "TTTTTTTTTT", "ATATATATAT"];
    let m = match_bkg(
        &["GGGGGGGGGG"],
        &universe,
        &[],
        &MatchOptions {
            n_per_target: 5,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(m.fallback_targets, 1);
    assert_eq!(m.universe_index.len(), 5);
    assert!(m.universe_index.iter().all(|&i| i < 3));
    // Without unique draws, a bin too small still needs a fallback, but its
    // first draws are distinct.
    let m = match_bkg(
        &["ACGTACGTAC"],
        &["ACGTACGTAC", "CAGTCAGTCA"],
        &[],
        &MatchOptions {
            n_per_target: 5,
            unique: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(m.fallback_targets, 1);
    assert_ne!(m.universe_index[0], m.universe_index[1]);
    assert!(m.universe_index.iter().all(|&i| i < 2));
}

#[test]
fn match_bkg_rejects_bad_input() {
    let t = ["ACGT"];
    let u = ["ACGT", "AAAA"];
    let o = MatchOptions::default();
    let none: [&str; 0] = [];
    assert_eq!(
        arg_msg(match_bkg(&none, &u, &[], &o)),
        "no target sequences given"
    );
    assert_eq!(
        arg_msg(match_bkg(&t, &none, &[], &o)),
        "the universe is empty"
    );
    assert_eq!(
        arg_msg(match_bkg(&t, &["ACGT", ""], &[], &o)),
        "the universe must not contain empty sequences"
    );
    assert_eq!(
        arg_msg(match_bkg(
            &t,
            &u,
            &[],
            &MatchOptions {
                n_per_target: 0,
                ..o.clone()
            }
        )),
        "n_per_target must be at least 1"
    );
    for bad in [
        MatchOptions {
            n_bins_gc: 1,
            ..o.clone()
        },
        MatchOptions {
            n_bins_length: 1,
            ..o.clone()
        },
    ] {
        assert_eq!(
            arg_msg(match_bkg(&t, &u, &[], &bad)),
            "every axis needs at least 2 bins"
        );
    }
    let cov = |targets: Vec<f64>, universe: Vec<f64>, n_bins| Covariate {
        name: "cpg".into(),
        targets,
        universe,
        n_bins,
        bin_type: BinType::Quantile,
    };
    assert_eq!(
        arg_msg(match_bkg(&t, &u, &[cov(vec![1.0], vec![1.0, 2.0], 1)], &o)),
        "every axis needs at least 2 bins"
    );
    assert_eq!(
        arg_msg(match_bkg(&t, &u, &[cov(vec![1.0], vec![1.0], 2)], &o)),
        "covariate 'cpg' needs one value per target and per universe sequence"
    );
    assert_eq!(
        arg_msg(match_bkg(&t, &u, &[cov(vec![], vec![1.0, 2.0], 2)], &o)),
        "covariate 'cpg' needs one value per target and per universe sequence"
    );
    for v in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            arg_msg(match_bkg(&t, &u, &[cov(vec![1.0], vec![v, 2.0], 2)], &o)),
            "covariate 'cpg' has missing or infinite values"
        );
    }
}

#[test]
fn genome_windows_respect_lengths_and_exclusions() {
    // The first chromosome is too short for the window; the third is
    // negative, so never drawn.
    let lens = [500, 1000, -50];
    let w = genome_windows(&lens, &[600], 200, &[], 3).unwrap();
    assert!(!w.is_empty() && w.len() < 200);
    for g in &w {
        assert_eq!(g.chrom, 1);
        assert_eq!(g.width(), 600);
        assert!(g.start >= 0 && g.end <= 1000);
    }
    assert_eq!(w, genome_windows(&lens, &[600], 200, &[], 3).unwrap());
    assert_ne!(w, genome_windows(&lens, &[600], 200, &[], 4).unwrap());
    // Widths are drawn from the targets; excluded windows are dropped, and
    // twice as many are tried so that enough remain.
    let excl = [GRange::new(0, 100, 300), GRange::new(0, 290, 350)];
    let w = genome_windows(&[1000], &[3, 7], 300, &excl, 1).unwrap();
    assert_eq!(w.len(), 300);
    assert!(w.iter().all(|g| g.width() == 3 || g.width() == 7));
    assert!(w.iter().any(|g| g.width() == 3) && w.iter().any(|g| g.width() == 7));
    assert!(w.iter().all(|g| g.end <= 100 || g.start >= 350));
    let w = genome_windows(&[20], &[10], 50, &[GRange::new(0, 10, 11)], 2).unwrap();
    // Only the window ending where the excluded position starts fits.
    assert!(!w.is_empty());
    assert!(w.iter().all(|g| *g == GRange::new(0, 0, 10)));
}

#[test]
fn genome_windows_reject_bad_input() {
    assert_eq!(
        arg_msg(genome_windows(&[100], &[], 10, &[], 1)),
        "target sequences must not be empty"
    );
    assert_eq!(
        arg_msg(genome_windows(&[100], &[5, 0], 10, &[], 1)),
        "target sequences must not be empty"
    );
    assert_eq!(
        arg_msg(genome_windows(&[100], &[5, 5, 5], 2, &[], 1)),
        "n_candidates must be at least the number of targets"
    );
    assert_eq!(
        arg_msg(genome_windows(&[0, -3], &[5], 10, &[], 1)),
        "no usable chromosomes"
    );
    assert_eq!(
        arg_msg(genome_windows(
            &[100],
            &[5],
            10,
            &[GRange::new(0, 0, 100)],
            1
        )),
        "could not sample enough windows outside the excluded regions (got 0, need 1)"
    );
}

// motif_finder

/// Sequences with a planted word in each, between random flanks.
fn planted(n: usize, word: &str, seed: u64) -> Vec<String> {
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut s: String = (0..60)
                .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
                .collect();
            if !word.is_empty() {
                let at = rng.random_range(0..60 - word.len());
                s.replace_range(at..at + word.len(), word);
            }
            s
        })
        .collect()
}

#[test]
fn finder_recovers_a_planted_word() {
    let pos = planted(40, "GTCACG", 1);
    let neg = planted(40, "", 2);
    let opts = FinderOptions {
        min_width: 8,
        max_width: 8,
        nmotifs: 2,
        hit_pvalue: 1e-3,
        ..Default::default()
    };
    let found = motif_finder(&pos, Some(&neg), &opts).unwrap();
    assert!(!found.is_empty());
    let f = &found[0];
    assert_eq!(f.rank, 1);
    assert_eq!(f.motif.name, "motif_1");
    assert_eq!(f.motif.alphabet, Alphabet::Dna);
    assert_eq!((f.n_pos, f.n_neg), (40, 40));
    assert_eq!(f.seqs_pos, 40);
    // The random flanks of the eight-wide search are trimmed away.
    let cons = f.motif.consensus.clone().unwrap();
    assert!(cons.contains("GTCACG") || cons.contains("CGTGAC"), "{cons}");
    assert!(f.motif.width() < 8, "{cons}");
    assert!(f.qvalue <= opts.qvalue);
    assert_eq!(f.motif.eval, Some(f.qvalue));
    assert_near(f.log10_pvalue, f.pvalue.log10(), "log10 pvalue");
    assert_near(f.log10_qvalue, f.qvalue.log10(), "log10 qvalue");
    for (i, g) in found.iter().enumerate() {
        assert_eq!(g.rank, i + 1);
        assert_eq!(g.motif.name, format!("motif_{}", i + 1));
    }
    assert!(found.windows(2).all(|w| w[0].pvalue <= w[1].pvalue));
    // Threads do not change the result.
    let threaded = motif_finder(
        &pos,
        Some(&neg),
        &FinderOptions {
            nthreads: 2,
            min_width: 7,
            ..opts.clone()
        },
    )
    .unwrap();
    let single = motif_finder(
        &pos,
        Some(&neg),
        &FinderOptions {
            min_width: 7,
            ..opts.clone()
        },
    )
    .unwrap();
    assert_eq!(threaded, single);
}

#[test]
fn finder_default_negatives_rna_and_missing_letters() {
    // Without negatives the positives are shuffled; the result is seeded.
    let pos: Vec<String> = planted(30, "GTCACG", 3)
        .into_iter()
        .map(|s| s.replace('T', "U"))
        .collect();
    let opts = FinderOptions {
        min_width: 6,
        max_width: 6,
        nmotifs: 1,
        hit_pvalue: 1e-3,
        ..Default::default()
    };
    let a = motif_finder(&pos, None, &opts).unwrap();
    assert!(!a.is_empty());
    assert_eq!(a, motif_finder(&pos, None, &opts).unwrap());
    for f in &a {
        assert_eq!(f.n_neg, 30);
        assert_eq!(f.motif.alphabet, Alphabet::Rna);
    }
    // Sequences missing a letter have its background floored, and with no
    // pseudocount the scores of unseen letters are floored too.
    let no_t: Vec<String> = planted(30, "GCCAGG", 4)
        .into_iter()
        .map(|s| s.replace('T', "A"))
        .collect();
    let neg: Vec<String> = planted(30, "", 5)
        .into_iter()
        .map(|s| s.replace('T', "C"))
        .collect();
    let b = motif_finder(
        &no_t,
        Some(&neg),
        &FinderOptions {
            pseudocount: 0,
            qvalue: 1.0,
            ..opts.clone()
        },
    )
    .unwrap();
    assert!(!b.is_empty());
    for f in &b {
        assert_eq!(f.motif.pseudocount, 0.0);
        let bkg = &f.motif.bkg;
        assert_near(bkg[3], 1e-3 / (1.0 + 1e-3), "floored T");
        assert_near(bkg.iter().sum::<f64>(), 1.0, "background sum");
    }
    // Without a strand to search beyond the given one, the result may
    // differ but stays well formed.
    let c = motif_finder(
        &no_t,
        Some(&neg),
        &FinderOptions {
            rc: false,
            ..opts.clone()
        },
    )
    .unwrap();
    assert!(c.iter().all(|f| f.n_pos == 30 && f.qvalue <= opts.qvalue));
    // Sequences with no shared word find nothing.
    let none = motif_finder(&planted(5, "", 6), Some(&planted(5, "", 7)), &opts).unwrap();
    assert!(none.is_empty());
}

#[test]
fn finder_deduplicates_across_widths_and_skips_short_sequences() {
    // Two words, each found at several widths; the copies covering the
    // same sites are dropped, leaving one motif per word.
    let mut pos = planted(40, "GTCACG", 8);
    for (k, s) in planted(40, "TTAGGC", 9).into_iter().enumerate() {
        pos[k].push_str(&s);
    }
    let neg = planted(40, "", 10);
    let opts = FinderOptions {
        min_width: 6,
        max_width: 8,
        nmotifs: 2,
        hit_pvalue: 1e-3,
        ..Default::default()
    };
    let found = motif_finder(&pos, Some(&neg), &opts).unwrap();
    let cons: Vec<String> = found
        .iter()
        .map(|f| f.motif.consensus.clone().unwrap())
        .collect();
    let has = |w: &str, rc: &str| {
        cons.iter()
            .filter(|c| c.contains(w) || c.contains(rc))
            .count()
    };
    assert_eq!(has("GTCACG", "CGTGAC"), 1, "{cons:?}");
    assert_eq!(has("TTAGGC", "GCCTAA"), 1, "{cons:?}");
    // Sequences shorter than every width count towards the totals but
    // change nothing else.
    let mut with_short = pos.clone();
    with_short.extend(["ACG".to_string(), "".to_string()]);
    let mut neg_short = neg.clone();
    neg_short.push("GT".to_string());
    let short = motif_finder(&with_short, Some(&neg_short), &opts).unwrap();
    assert_eq!(short.len(), found.len());
    for (a, b) in found.iter().zip(&short) {
        assert_eq!(a.motif.consensus, b.motif.consensus);
        assert_eq!((b.n_pos, b.n_neg), (42, 41));
        assert_eq!((a.seqs_pos, a.sites_pos), (b.seqs_pos, b.sites_pos));
    }
    // With nothing as long as the narrowest width, there is no word to
    // start from.
    let tiny = motif_finder(&["ACG", "GTA"], Some(&["ACGT"][..]), &opts).unwrap();
    assert!(tiny.is_empty());
}

#[test]
fn finder_rejects_bad_input() {
    let pos = ["ACGTACGTAC"];
    let o = FinderOptions::default();
    let none: [&str; 0] = [];
    assert_eq!(arg_msg(motif_finder(&none, None, &o)), "no sequences given");
    let widths = "widths must satisfy 3 <= min_width <= max_width <= 30";
    for (min_width, max_width) in [(2, 10), (6, 31), (10, 9)] {
        assert_eq!(
            arg_msg(motif_finder(
                &pos,
                None,
                &FinderOptions {
                    min_width,
                    max_width,
                    ..o.clone()
                }
            )),
            widths
        );
    }
    assert_eq!(
        arg_msg(motif_finder(
            &pos,
            None,
            &FinderOptions {
                nmotifs: 0,
                ..o.clone()
            }
        )),
        "nmotifs must be at least 1"
    );
    let probs = "hit_pvalue must be in (0, 1); stop_pvalue and qvalue in (0, 1]; \
                 dedup_overlap in [0, 1]";
    for bad in [
        FinderOptions {
            hit_pvalue: 1.0,
            ..o.clone()
        },
        FinderOptions {
            hit_pvalue: 0.0,
            ..o.clone()
        },
        FinderOptions {
            stop_pvalue: 0.0,
            ..o.clone()
        },
        FinderOptions {
            stop_pvalue: 1.5,
            ..o.clone()
        },
        FinderOptions {
            qvalue: 0.0,
            ..o.clone()
        },
        FinderOptions {
            dedup_overlap: 1.1,
            ..o.clone()
        },
        FinderOptions {
            dedup_overlap: -0.1,
            ..o.clone()
        },
    ] {
        assert_eq!(arg_msg(motif_finder(&pos, None, &bad)), probs);
    }
    assert_eq!(
        arg_msg(motif_finder(&pos, Some(&none), &o)),
        "no negative sequences"
    );
    assert_eq!(
        arg_msg(motif_finder(
            &pos,
            None,
            &FinderOptions {
                shuffle_k: 0,
                ..o.clone()
            }
        )),
        "k must be at least 1"
    );
}

// Genomes on disk

/// Compress text in BGZF blocks of at most `block` bytes, with the empty
/// end-of-file block.
fn bgzip(text: &[u8], block: usize) -> Vec<u8> {
    let mut out = Vec::new();
    let put = |data: &[u8], out: &mut Vec<u8>| {
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        let deflated = enc.finish().unwrap();
        let bsize = (deflated.len() + 25) as u16;
        out.extend_from_slice(&[
            0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
        ]);
        out.extend_from_slice(&bsize.to_le_bytes());
        out.extend_from_slice(&deflated);
        let mut crc = flate2::Crc::new();
        crc.update(data);
        out.extend_from_slice(&crc.sum().to_le_bytes());
        out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    };
    for chunk in text.chunks(block) {
        put(chunk, &mut out);
    }
    put(&[], &mut out);
    out
}

/// A 2bit file of one or more sequences, in either byte order and either
/// version, with N runs and lower case runs as its blocks.
fn twobit(chroms: &[(&str, &str)], big: bool, version: u32) -> Vec<u8> {
    let u32b = |x: u32| {
        if big {
            x.to_be_bytes()
        } else {
            x.to_le_bytes()
        }
    };
    let runs = |seq: &str, f: &dyn Fn(u8) -> bool| {
        let mut out: Vec<(u32, u32)> = Vec::new();
        for (i, b) in seq.bytes().enumerate() {
            if f(b) {
                match out.last_mut() {
                    Some(r) if (r.0 + r.1) as usize == i => r.1 += 1,
                    _ => out.push((i as u32, 1)),
                }
            }
        }
        out
    };
    let mut records = Vec::new();
    for (_, seq) in chroms {
        let mut r = Vec::new();
        r.extend_from_slice(&u32b(seq.len() as u32));
        for blocks in [
            runs(seq, &|b| b == b'N' || b == b'n'),
            runs(seq, &|b| b.is_ascii_lowercase()),
        ] {
            r.extend_from_slice(&u32b(blocks.len() as u32));
            for b in &blocks {
                r.extend_from_slice(&u32b(b.0));
            }
            for b in &blocks {
                r.extend_from_slice(&u32b(b.1));
            }
        }
        r.extend_from_slice(&u32b(0));
        let mut packed = vec![0u8; seq.len().div_ceil(4)];
        for (i, b) in seq.bytes().enumerate() {
            let code = match b.to_ascii_uppercase() {
                b'C' => 1,
                b'A' => 2,
                b'G' => 3,
                _ => 0,
            };
            packed[i / 4] |= code << (6 - 2 * (i % 4));
        }
        r.extend_from_slice(&packed);
        records.push(r);
    }
    let mut out = Vec::new();
    out.extend_from_slice(&u32b(0x1A412743));
    out.extend_from_slice(&u32b(version));
    out.extend_from_slice(&u32b(chroms.len() as u32));
    out.extend_from_slice(&u32b(0));
    let off_bytes = if version == 1 { 8 } else { 4 };
    let index_len: usize = chroms.iter().map(|c| 1 + c.0.len() + off_bytes).sum();
    let mut offset = 16 + index_len;
    for ((name, _), r) in chroms.iter().zip(&records) {
        out.push(name.len() as u8);
        out.extend_from_slice(name.as_bytes());
        if version == 1 {
            let o = offset as u64;
            out.extend_from_slice(&if big {
                o.to_be_bytes()
            } else {
                o.to_le_bytes()
            });
        } else {
            out.extend_from_slice(&u32b(offset as u32));
        }
        offset += r.len();
    }
    for r in records {
        out.extend_from_slice(&r);
    }
    out
}

#[test]
fn twobit_byte_orders_and_versions() {
    let dir = tempdir();
    let chroms = [("chrA", "ACGTnnNNacgtTTGCAg"), ("b", "GATTACA")];
    for big in [false, true] {
        for version in [0, 1] {
            let path = dir.join(format!("g-{big}-{version}.2bit"));
            std::fs::write(&path, twobit(&chroms, big, version)).unwrap();
            let g = Genome::open(&path).unwrap();
            assert_eq!(g.names(), ["chrA", "b"]);
            assert_eq!(g.lengths(), [18, 7]);
            assert_eq!(g.len(), 2);
            assert!(!g.is_empty());
            // N runs read as N, and soft-masked runs, N or not, in lower
            // case.
            assert_eq!(g.fetch("chrA", 0, 18).unwrap(), "ACGTnnNNacgtTTGCAg");
            assert_eq!(g.fetch("chrA", 5, 9).unwrap(), "nNNa");
            assert_eq!(g.fetch("chrA", 17, 18).unwrap(), "g");
            assert_eq!(g.fetch_at(1, 1, 4).unwrap(), "ATT");
            assert_eq!(
                g.fetch_many(
                    &[
                        GRange::new(1, 3, 7),
                        GRange::new(0, 9, 13),
                        GRange::new(1, 0, 0)
                    ],
                    2
                )
                .unwrap(),
                ["TACA", "cgtT", ""]
            );
        }
    }
    let path = dir.join("v2.2bit");
    std::fs::write(&path, twobit(&chroms, false, 2)).unwrap();
    let e = Genome::open(&path).unwrap_err().to_string();
    assert!(e.contains("unknown 2bit version 2"), "{e}");
}

#[test]
fn genome_fetch_errors() {
    let dir = tempdir();
    let path = dir.join("g.fa");
    std::fs::write(&path, ">chr1 first\nACGTA\nCG\n>chr2\nTTTT\n").unwrap();
    let g = Genome::open(&path).unwrap();
    assert_eq!(g.names(), ["chr1", "chr2"]);
    assert_eq!(g.lengths(), [7, 4]);
    assert_eq!(g.position("chr2"), Some(1));
    assert_eq!(g.position("chr3"), None);
    assert_eq!(g.fetch("chr1", 3, 7).unwrap(), "TACG");
    assert_eq!(arg_msg(g.fetch_at(2, 0, 1)), "sequence 2 is out of range");
    assert_eq!(
        arg_msg(g.fetch("chr2", 3, 5)),
        "positions 3 to 5 are outside chr2, which has 4 letters"
    );
    assert_eq!(
        arg_msg(g.fetch("chr2", 3, 2)),
        "positions 3 to 2 are outside chr2, which has 4 letters"
    );
    assert_eq!(
        arg_msg(g.fetch("chrX", 0, 1)),
        "no sequence called 'chrX' in the genome"
    );
    for (s, e) in [(-1, 3), (3, 2)] {
        assert_eq!(
            arg_msg(g.fetch_many(&[GRange::new(0, 0, 1), GRange::new(0, s, e)], 1)),
            format!("positions {s} to {e} are not a valid range")
        );
    }
    assert_eq!(
        arg_msg(g.fetch_many(&[GRange::new(5, 0, 1)], 1)),
        "sequence 5 is out of range"
    );
    // Strands are ignored.
    assert_eq!(
        g.fetch_many(&[stranded(1, 0, 2, GStrand::Minus)], 1)
            .unwrap(),
        ["TT"]
    );
    let missing = dir.join("missing.fa");
    let e = arg_msg(Genome::open(&missing));
    assert!(
        e.starts_with(&format!(
            "could not read the genome {}: ",
            missing.display()
        )),
        "{e}"
    );
}

#[test]
fn fasta_layout_errors() {
    let dir = tempdir();
    let open = |name: &str, text: &str| {
        let path = dir.join(name);
        std::fs::write(&path, text).unwrap();
        Genome::open(&path)
    };
    let e = arg_msg(open("before.fa", "ACGT\n>a\nAC\n"));
    assert!(e.contains("sequence before the first FASTA header"), "{e}");
    // A line longer than the first, or ending differently, breaks the
    // layout.
    for (i, text) in [">a\nACGT\nACGTA\n", ">a\nACGT\r\nACGT\n"]
        .iter()
        .enumerate()
    {
        let e = arg_msg(open(&format!("ragged{i}.fa"), text));
        assert!(
            e.contains("the lines of a differ in length"),
            "{text:?}: {e}"
        );
    }
    // An empty file holds no sequences.
    let g = open("empty.fa", "").unwrap();
    assert!(g.is_empty());
    assert_eq!(g.len(), 0);
    // A header with no sequence is a sequence of length 0; blank lines
    // between records are skipped.
    let g = open("blank.fa", ">a\n\n>b\nAC\n\n").unwrap();
    assert_eq!(g.lengths(), [0, 2]);
    assert_eq!(g.fetch("a", 0, 0).unwrap(), "");
    assert_eq!(g.fetch("b", 0, 2).unwrap(), "AC");
}

#[test]
fn fai_files_are_trusted_and_checked() {
    let dir = tempdir();
    let path = dir.join("g.fa");
    std::fs::write(&path, ">a\nACGT\nAC\n>b\nGGG\n").unwrap();
    // Blank lines in the index are skipped.
    std::fs::write(dir.join("g.fa.fai"), "a\t6\t3\t4\t5\n\nb\t3\t14\t3\t4\n").unwrap();
    let g = Genome::open(&path).unwrap();
    assert_eq!(g.lengths(), [6, 3]);
    assert_eq!(g.fetch("a", 2, 6).unwrap(), "GTAC");
    assert_eq!(g.fetch("b", 0, 3).unwrap(), "GGG");
    std::fs::write(dir.join("g.fa.fai"), "a\tsix\t3\t4\t5\n").unwrap();
    let e = arg_msg(Genome::open(&path));
    assert!(e.contains("bad line in"), "{e}");
    std::fs::write(dir.join("g.fa.fai"), "a\t6\t3\n").unwrap();
    assert!(arg_msg(Genome::open(&path)).contains("bad line in"));
}

#[test]
fn bgzip_edge_cases() {
    let dir = tempdir();
    let text = b">a\nACGTACGTAC\nGT\n>b\nTTTT\n";
    let path = dir.join("g.fa.gz");
    std::fs::write(&path, bgzip(text, 7)).unwrap();
    let g = Genome::open(&path).unwrap();
    assert_eq!(g.fetch("a", 0, 12).unwrap(), "ACGTACGTACGT");
    assert_eq!(g.fetch("b", 1, 4).unwrap(), "TTT");
    // Pieces spread over several blocks, read in any order.
    assert_eq!(
        g.fetch_many(
            &[
                GRange::new(1, 0, 4),
                GRange::new(0, 9, 12),
                GRange::new(0, 0, 2)
            ],
            1
        )
        .unwrap(),
        ["TTTT", "CGT", "AC"]
    );
    // A .fai next to it is used as it is; one claiming more sequence than
    // the file holds fails when that part is read.
    std::fs::write(
        dir.join("g.fa.gz.fai"),
        "a\t12\t3\t10\t11\nb\t40\t20\t4\t5\n",
    )
    .unwrap();
    let g = Genome::open(&path).unwrap();
    assert_eq!(g.lengths(), [12, 40]);
    assert_eq!(g.fetch("a", 8, 12).unwrap(), "ACGT");
    let e = arg_msg(g.fetch("b", 0, 40));
    assert!(e.contains("the file ends before the sequence does"), "{e}");
    // A truncated .gzi is an error.
    for gzi in [&[1u8, 0, 0][..], &5u64.to_le_bytes()[..]] {
        std::fs::write(dir.join("g.fa.gz.gzi"), gzi).unwrap();
        let e = arg_msg(Genome::open(&path));
        assert!(e.contains("the .gzi index is truncated"), "{e}");
    }
    // Gzip with an extra field that is not BGZF's is not bgzip.
    let mut fake = vec![
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'X', b'Y', 2, 0, 0, 0,
    ];
    fake.extend_from_slice(&[0; 10]);
    let fake_path = dir.join("fake.fa.gz");
    std::fs::write(&fake_path, fake).unwrap();
    let e = arg_msg(Genome::open(&fake_path));
    assert!(e.contains("not with bgzip"), "{e}");
}

#[test]
fn fai_with_zero_line_bases_is_an_error() {
    let dir = tempdir();
    let path = dir.join("g.fa");
    std::fs::write(&path, ">a\nACGT\n").unwrap();
    std::fs::write(dir.join("g.fa.fai"), "a\t4\t3\t0\t0\n").unwrap();
    // The index claims four letters on lines of no letters; reading them
    // should fail with an error rather than panic.
    let r = std::panic::catch_unwind(|| Genome::open(&path).and_then(|g| g.fetch("a", 0, 2)));
    assert!(matches!(r, Ok(Err(Error::Argument(_)))), "{r:?}");
}

#[test]
fn corrupt_bgzf_block_size_is_an_error() {
    let dir = tempdir();
    let path = dir.join("corrupt.fa.gz");
    // A BGZF header whose block size field says the block is one byte long.
    let mut data = vec![
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0, 0,
    ];
    data.extend_from_slice(&[0; 10]);
    std::fs::write(&path, data).unwrap();
    let r = std::panic::catch_unwind(|| Genome::open(&path));
    assert!(matches!(r, Ok(Err(Error::Argument(_)))), "{r:?}");
}
