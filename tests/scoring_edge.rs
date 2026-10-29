//! Edge cases and exact checks for scoring, P-values, column metrics, motif
//! comparison, clustering, the statistical helpers and scanning. Expected
//! values come from exhaustive enumeration of small motifs or are worked out
//! by hand.

use std::collections::{BTreeMap, BTreeSet};
use std::slice::from_ref;

use universalmotif::alphabet::AMINO_ACIDS;
use universalmotif::cdf::ScoreDistribution;
use universalmotif::cluster::{Linkage, Tree, TreeOptions, hclust, motif_tree};
use universalmotif::compare::{
    CompareOptions, NullModel, Overhang, anchor_layout, compare_matrix, compare_motifs,
};
use universalmotif::convert::convert_type;
use universalmotif::create::{CreateOptions, from_strings, random};
use universalmotif::metrics::{
    Average, Column, ColumnMetric, average, average_ic, compare_columns,
};
use universalmotif::ops::{add_gap, reverse_complement, ungap};
use universalmotif::packed::EncodedSequences;
use universalmotif::pvalue::{
    PvalueOptions, log_pvalue_text, min_pvalues, motif_pvalue, pvalue_to_score,
    score_distributions, scoring_matrix,
};
use universalmotif::scan::{
    Hit, MatchValue, NoOverlaps, OverlapPriority, ScanOptions, motif_matches, scan_encoded,
    scan_fasta, scan_sequences, scan_table,
};
use universalmotif::scoring::{
    ThresholdType, get_matches, get_scores, motif_range, motif_score, prob_match, prob_match_bkg,
    score_match,
};
use universalmotif::shuffle::create_sequences;
use universalmotif::stats::{
    dbinom, dhyper, fisher_greater, ln_fisher_greater, ln_p_adjust_bh, median, p_adjust_bh,
    p_adjust_bonferroni, pbinom_upper, phyper,
};
use universalmotif::{Alphabet, Error, Matrix, Motif, MotifType, Result, Strand};

const DNA: [char; 4] = ['A', 'C', 'G', 'T'];

fn site_opts() -> CreateOptions {
    CreateOptions {
        nsites: Some(20.0),
        pseudocount: 1.0,
        ..Default::default()
    }
}

fn sites(s: &[&str]) -> Motif {
    from_strings(s, None, &site_opts()).unwrap()
}

fn sites_in(s: &[&str], alphabet: Alphabet) -> Motif {
    from_strings(s, Some(alphabet), &site_opts()).unwrap()
}

fn pwm_motif(cols: &[&[f64]], alphabet: Alphabet, bkg: Option<Vec<f64>>) -> Motif {
    Motif::builder(Matrix::from_columns(cols).unwrap(), alphabet)
        .motif_type(Some(MotifType::Pwm))
        .bkg(bkg)
        .build()
        .unwrap()
}

/// A probability matrix kept exactly as given, zeros included.
fn raw_ppm(cols: &[&[f64]], alphabet: Alphabet) -> Motif {
    Motif::builder(Matrix::from_columns(cols).unwrap(), alphabet)
        .motif_type(Some(MotifType::Ppm))
        .pseudocount(0.0)
        .nsites(Some(10.0))
        .build()
        .unwrap()
}

/// A PWM entry on the integer score grid, truncated toward zero.
fn grid(x: f64) -> i64 {
    (x * 1000.0) as i32 as i64
}

/// Every word of `width` letters as row indices, the last position varying
/// fastest.
fn all_words(nrow: usize, width: usize) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new()];
    for _ in 0..width {
        out = out
            .into_iter()
            .flat_map(|w: Vec<usize>| {
                (0..nrow).map(move |r| {
                    let mut v = w.clone();
                    v.push(r);
                    v
                })
            })
            .collect();
    }
    out
}

fn word_text(rows: &[usize], alphabet: &Alphabet) -> String {
    rows.iter().map(|&r| alphabet.letters()[r]).collect()
}

/// The grid score and background probability of every word of a PWM.
fn enumerate(pwm: &Matrix, bkg: &[f64]) -> Vec<(i64, f64)> {
    all_words(pwm.nrow(), pwm.ncol())
        .iter()
        .map(|w| {
            let s = w
                .iter()
                .enumerate()
                .map(|(j, &r)| grid(pwm.get(r, j)))
                .sum();
            let p = w.iter().map(|&r| bkg[r]).product();
            (s, p)
        })
        .collect()
}

/// The probability of a grid score of at least `g`.
fn tail(all: &[(i64, f64)], g: i64) -> f64 {
    all.iter().filter(|x| x.0 >= g).map(|x| x.1).sum()
}

fn near(got: f64, want: f64, tol: f64) -> bool {
    if got.is_infinite() || want.is_infinite() || got.is_nan() || want.is_nan() {
        return got == want || (got.is_nan() && want.is_nan());
    }
    (got - want).abs() <= tol * f64::max(1.0, want.abs())
}

#[track_caller]
fn assert_near(got: f64, want: f64, tol: f64, ctx: &str) {
    assert!(near(got, want, tol), "{ctx}: got {got}, want {want}");
}

/// The message of an argument error.
#[track_caller]
fn arg_msg<T: std::fmt::Debug>(r: Result<T>) -> String {
    match r {
        Err(Error::Argument(m)) => m,
        other => panic!("expected an argument error, got {other:?}"),
    }
}

fn rc_word(w: &str, rna: bool) -> String {
    w.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' if rna => 'U',
            'A' => 'T',
            'T' | 'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            o => o,
        })
        .collect()
}

fn pv(m: &Motif, score: f64) -> f64 {
    motif_pvalue(from_ref(m), &[(0, score)], &PvalueOptions::default()).unwrap()[0]
}

fn threshold(m: &Motif, p: f64) -> f64 {
    pvalue_to_score(from_ref(m), &[(0, p)], &PvalueOptions::default()).unwrap()[0]
}

// P-values and the score distribution.

#[test]
fn pvalues_equal_exhaustive_enumeration() {
    let bkg = vec![0.1, 0.2, 0.3, 0.4];
    let m = pwm_motif(
        &[
            &[1.2345, -0.5, 0.25, -2.0009],
            &[-1.0, 0.75, 0.3333, 0.0],
            &[0.5, 0.5, -0.4999, -3.1],
            &[2.0, -1.5, 0.001, -0.0015],
        ],
        Alphabet::Dna,
        Some(bkg.clone()),
    );
    let all = enumerate(&m.matrix, &bkg);
    let lo = all.iter().map(|x| x.0).min().unwrap();
    let hi = all.iter().map(|x| x.0).max().unwrap();

    // Every grid score, and every score halfway between two grid points,
    // which must count as the next grid point up.
    let mut queries = Vec::new();
    let mut want = Vec::new();
    for g in lo - 3..=hi + 3 {
        queries.push((0, g as f64 / 1000.0));
        want.push(tail(&all, g));
        queries.push((0, (g as f64 + 0.5) / 1000.0));
        want.push(tail(&all, g + 1));
    }
    let motifs = [m.clone()];
    let got = motif_pvalue(&motifs, &queries, &PvalueOptions::default()).unwrap();
    for ((q, g), w) in queries.iter().zip(&got).zip(&want) {
        assert_near(*g, *w, 1e-12, &format!("score {}", q.1));
    }
    // More threads give the same answers.
    let opts = PvalueOptions {
        nthreads: 0,
        ..Default::default()
    };
    assert_eq!(motif_pvalue(&motifs, &queries, &opts).unwrap(), got);

    // The distribution itself, and its upper tail.
    let dist = ScoreDistribution::new(&m.matrix, &bkg);
    assert_eq!(dist.score_range(), (lo as f64 / 1000.0, hi as f64 / 1000.0));
    assert_eq!(
        motif_range(&m, false).unwrap(),
        (lo as f64 / 1000.0, hi as f64 / 1000.0)
    );
    for (q, g) in queries.iter().zip(&got) {
        assert_eq!(dist.pvalue(q.1), *g);
    }
    let from = (lo + hi) / 2;
    let t = dist.tail(from as f64 / 1000.0);
    assert_eq!(t.len(), (hi - from + 1) as usize);
    assert!(!t.is_empty());
    for g in from..=hi + 3 {
        let s = g as f64 / 1000.0;
        assert_eq!(t.pvalue(s), dist.pvalue(s), "tail at {s}");
    }
    assert!(dist.tail((hi + 1) as f64 / 1000.0).is_empty());
    assert!(dist.tail(-1e6).len() >= (hi - lo + 1) as usize);
    // A background of zeros gives no probability to any score in range.
    let none = ScoreDistribution::new(&m.matrix, &[0.0; 4]);
    assert_eq!(none.pvalue(lo as f64 / 1000.0), 0.0);
    assert_eq!(none.pvalue(-f64::MAX), 1.0);
}

#[test]
fn score_thresholds_invert_pvalues() {
    let bkg = vec![0.1, 0.2, 0.3, 0.4];
    let m = pwm_motif(
        &[
            &[1.2345, -0.5, 0.25, -2.0009],
            &[-1.0, 0.75, 0.3333, 0.0],
            &[0.5, 0.5, -0.4999, -3.1],
        ],
        Alphabet::Dna,
        Some(bkg.clone()),
    );
    let all = enumerate(&m.matrix, &bkg);
    let lo = all.iter().map(|x| x.0).min().unwrap();
    let hi = all.iter().map(|x| x.0).max().unwrap();

    // The best word: the top letter of each column, ties included.
    let best: f64 = m
        .matrix
        .columns()
        .map(|c| {
            let top = c.iter().map(|&v| grid(v)).max().unwrap();
            c.iter()
                .zip(&bkg)
                .filter(|(v, _)| grid(**v) == top)
                .map(|(_, b)| b)
                .sum::<f64>()
        })
        .product();
    let min_p = min_pvalues(from_ref(&m), &PvalueOptions::default()).unwrap()[0];
    assert_near(min_p, best, 1e-12, "smallest P-value");
    assert_near(
        pv(&m, hi as f64 / 1000.0),
        best,
        1e-12,
        "P-value of the best",
    );

    for p in [0.37, 0.0123, 0.00271, 4.1e-4, best * 1.5] {
        let s = threshold(&m, p);
        let Some(g) = (lo..=hi).find(|&g| tail(&all, g) <= p) else {
            assert!(p < best);
            assert_eq!(s, f64::INFINITY, "threshold for {p}");
            continue;
        };
        assert_eq!(s, g as f64 / 1000.0, "threshold for {p}");
        assert!(pv(&m, s) <= p);
        assert!(pv(&m, s - 0.001) > p, "one grid point lower exceeds {p}");
    }
    // At P = 1 or above every score qualifies; below the best word's
    // P-value none does.
    for p in [1.0, 2.0] {
        assert_eq!(threshold(&m, p), lo as f64 / 1000.0);
    }
    assert_eq!(threshold(&m, best * 0.5), f64::INFINITY);
    assert_eq!(threshold(&m, 0.0), f64::INFINITY);
    assert_eq!(pv(&m, f64::MAX), 0.0);
    assert_eq!(pv(&m, -f64::MAX), 1.0);
}

#[test]
fn pvalues_of_other_alphabets_equal_enumeration() {
    let aa = sites_in(&["MKV", "MRV", "LKV", "MKI"], Alphabet::Aa);
    let custom = sites_in(
        &["xyzzy", "xyxzy", "zyzzy"],
        Alphabet::custom("xyz".chars()).unwrap(),
    );
    for m in [aa, custom] {
        let sm = scoring_matrix(&m, None).unwrap();
        assert!(!sm.pseudocount_added);
        assert_eq!(sm.bkg, m.bkg);
        let all = enumerate(&sm.pwm, &sm.bkg);
        let scores: BTreeSet<i64> = all.iter().map(|x| x.0).collect();
        let queries: Vec<(usize, f64)> = scores.iter().map(|&g| (0, g as f64 / 1000.0)).collect();
        let got = motif_pvalue(from_ref(&m), &queries, &PvalueOptions::default()).unwrap();
        for (&g, p) in scores.iter().zip(got) {
            assert_near(p, tail(&all, g), 1e-12, &format!("{} at {g}", m.alphabet));
        }
    }
}

#[test]
fn backgrounds_given_in_the_options() {
    let m = sites(&["TGACTCA", "TGAGTCA", "TGACTAA"]);
    let bkg = vec![0.3, 0.2, 0.2, 0.3];
    let sm = scoring_matrix(&m, Some(&bkg)).unwrap();
    assert_eq!(sm.bkg, bkg);
    assert!(!sm.pseudocount_added);
    let all = enumerate(&sm.pwm, &bkg);
    let (lo, hi) = (
        all.iter().map(|x| x.0).min().unwrap(),
        all.iter().map(|x| x.0).max().unwrap(),
    );
    let gs: Vec<i64> = (0..20).map(|i| lo + (hi - lo) * i / 19).collect();
    let queries: Vec<(usize, f64)> = gs.iter().map(|&g| (0, g as f64 / 1000.0)).collect();
    let one = PvalueOptions {
        bkg: vec![bkg.clone()],
        ..Default::default()
    };
    let got = motif_pvalue(from_ref(&m), &queries, &one).unwrap();
    for (g, p) in gs.iter().zip(&got) {
        assert_near(*p, tail(&all, *g), 1e-12, &format!("score {g}"));
    }

    // One background per motif answers each motif as on its own.
    let other = sites(&["CACGTG", "CACATG"]);
    let per = PvalueOptions {
        bkg: vec![bkg.clone(), vec![0.25; 4]],
        nthreads: 2,
    };
    let pair = [m.clone(), other.clone()];
    let q = [(0, 3.0), (1, 2.0), (0, 5.0)];
    let both = motif_pvalue(&pair, &q, &per).unwrap();
    assert_eq!(
        both[0],
        motif_pvalue(from_ref(&m), &[(0, 3.0)], &one).unwrap()[0]
    );
    assert_eq!(
        both[2],
        motif_pvalue(from_ref(&m), &[(0, 5.0)], &one).unwrap()[0]
    );
    assert_eq!(both[1], pv(&other, 2.0));
    let dists = score_distributions(&pair, &per).unwrap();
    assert_eq!(dists[0].pvalue(3.0), both[0]);
    assert_eq!(dists[1].pvalue(2.0), both[1]);
    let shared = score_distributions(&pair, &one).unwrap();
    assert_eq!(shared[0], dists[0]);
    assert_eq!(
        shared[1].pvalue(2.0),
        motif_pvalue(&pair, &[(1, 2.0)], &one).unwrap()[0]
    );
    assert_eq!(
        score_distributions(&pair, &PvalueOptions::default()).unwrap()[1].pvalue(2.0),
        both[1]
    );
    let thr = pvalue_to_score(&pair, &[(1, 1e-3), (0, 1e-3)], &per).unwrap();
    assert_eq!(thr[0], threshold(&other, 1e-3));
    assert_eq!(
        thr[1],
        pvalue_to_score(from_ref(&m), &[(0, 1e-3)], &one).unwrap()[0]
    );

    // A finite PWM against its own background is used as it is.
    let w = convert_type(&m, MotifType::Pwm, None, false).unwrap();
    let own = scoring_matrix(&w, None).unwrap();
    assert_eq!(own.pwm, w.matrix);
    // Zeros without a pseudocount are given one.
    let z = raw_ppm(
        &[&[1.0, 0.0, 0.0, 0.0], &[0.5, 0.5, 0.0, 0.0]],
        Alphabet::Dna,
    );
    let sz = scoring_matrix(&z, None).unwrap();
    assert!(sz.pseudocount_added);
    assert!(sz.pwm.as_slice().iter().all(|v| v.is_finite()));
}

#[test]
fn pvalue_arguments_are_checked() {
    let m = sites(&["TGACTCA"]);
    let motifs = [m.clone()];
    let opts = PvalueOptions::default();
    assert_eq!(
        arg_msg(scoring_matrix(&m, Some(&[0.5, 0.5]))),
        "the background must have 4 values, one per letter"
    );
    assert_eq!(
        arg_msg(motif_pvalue(&motifs, &[(1, 0.0)], &opts)),
        "query refers to motif 1, but there are 1"
    );
    for v in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            arg_msg(motif_pvalue(&motifs, &[(0, v)], &opts)),
            "scores and P-values must be finite numbers"
        );
        assert_eq!(
            arg_msg(pvalue_to_score(&motifs, &[(0, v)], &opts)),
            "scores and P-values must be finite numbers"
        );
    }
    let three = [m.clone(), m.clone(), m.clone()];
    let two_bkg = PvalueOptions {
        bkg: vec![vec![0.25; 4]; 2],
        ..Default::default()
    };
    let msg = "give one background for all motifs or one per motif";
    assert_eq!(arg_msg(motif_pvalue(&three, &[(0, 1.0)], &two_bkg)), msg);
    assert_eq!(arg_msg(score_distributions(&three, &two_bkg)), msg);
    assert_eq!(arg_msg(min_pvalues(&three, &two_bkg)), msg);
    let short = PvalueOptions {
        bkg: vec![vec![0.5; 2]],
        ..Default::default()
    };
    assert_eq!(
        arg_msg(motif_pvalue(&motifs, &[(0, 1.0)], &short)),
        "the background must have 4 values, one per letter"
    );
    assert_eq!(
        arg_msg(score_distributions(&motifs, &short)),
        "the background must have 4 values, one per letter"
    );
    assert!(motif_pvalue(&motifs, &[], &opts).unwrap().is_empty());
    assert!(min_pvalues(&[], &opts).unwrap().is_empty());
}

#[test]
fn log_pvalues_from_text() {
    let ln10 = std::f64::consts::LN_10;
    assert_near(log_pvalue_text("0.5").unwrap(), 0.5f64.ln(), 1e-15, "0.5");
    assert_near(
        log_pvalue_text(" 1e-3 ").unwrap(),
        -3.0 * ln10,
        1e-14,
        "1e-3",
    );
    assert_near(
        log_pvalue_text("1e-400").unwrap(),
        -400.0 * ln10,
        1e-14,
        "1e-400",
    );
    assert_near(
        log_pvalue_text("2.5E-500").unwrap(),
        2.5f64.ln() - 500.0 * ln10,
        1e-14,
        "2.5E-500",
    );
    assert_eq!(log_pvalue_text("0").unwrap(), f64::NEG_INFINITY);
    assert_eq!(log_pvalue_text("1").unwrap(), 0.0);
    for bad in ["abc", "-1", "NaN", "", "1e"] {
        assert_eq!(
            arg_msg(log_pvalue_text(bad)),
            format!("'{bad}' is not a P-value")
        );
    }
}

// Scores of single matches.

#[test]
fn scores_of_every_word_agree() {
    let m = sites(&["TGAC", "TGAG", "AGAC"]);
    let pwm = convert_type(&m, MotifType::Pwm, None, false).unwrap();
    let words: Vec<String> = all_words(4, 4)
        .iter()
        .map(|w| word_text(w, &Alphabet::Dna))
        .collect();
    let all = enumerate(&pwm.matrix, &m.bkg);
    let scored = score_match(&m, &words, false).unwrap();
    for ((w, s), (g, _)) in words.iter().zip(&scored).zip(&all) {
        assert_eq!(*s, *g as f64 / 1000.0, "{w}");
    }
    // Lower case DNA scores the same.
    let lower: Vec<String> = words.iter().map(|w| w.to_lowercase()).collect();
    assert_eq!(score_match(&m, &lower, false).unwrap(), scored);

    // get_scores lists the same multiset, highest first.
    let mut sorted = scored.clone();
    sorted.sort_by(|a, b| b.total_cmp(a));
    assert_eq!(get_scores(&m, false).unwrap(), sorted);
    assert_eq!(
        motif_range(&m, false).unwrap(),
        (*sorted.last().unwrap(), sorted[0])
    );

    // get_matches gives exactly the words at or above each reachable score.
    let reachable: BTreeSet<i64> = all.iter().map(|x| x.0).collect();
    let mut wrong = Vec::new();
    for &g in &reachable {
        let s = g as f64 / 1000.0;
        let got = get_matches(&m, s, false).unwrap();
        let set: BTreeSet<String> = got.iter().cloned().collect();
        assert_eq!(set.len(), got.len(), "duplicates at {s}");
        let want: BTreeSet<String> = words
            .iter()
            .zip(&all)
            .filter(|(_, x)| x.0 >= g)
            .map(|(w, _)| w.clone())
            .collect();
        if set != want {
            wrong.push(s);
        }
    }
    assert!(
        wrong.is_empty(),
        "get_matches differs from enumeration at {wrong:?}"
    );
    assert!(get_matches(&m, f64::INFINITY, false).unwrap().is_empty());
}

#[test]
fn match_scores_outside_the_range() {
    let m = sites(&["TGAC", "TGAG"]);
    let (lo, hi) = motif_range(&m, false).unwrap();
    assert_eq!(
        arg_msg(get_matches(&m, hi + 0.01, false)),
        format!("the score is above the highest possible score, {hi}")
    );
    assert_eq!(
        arg_msg(get_matches(&m, lo - 0.01, false)),
        format!("the score is below the lowest possible score, {lo}")
    );
    assert_eq!(get_matches(&m, lo, false).unwrap().len(), 256);
}

#[test]
fn matches_between_grid_points() {
    // Words score 1 (AA), 0.75 (CA) or less.
    let m = pwm_motif(
        &[&[1.0, 0.75, -1.0, -2.0], &[0.0, -1.0, -1.0, -1.0]],
        Alphabet::Dna,
        None,
    );
    assert_eq!(score_match(&m, &["AA", "CA"], false).unwrap(), [1.0, 0.75]);
    // Only AA scores at least 0.7505, as scanning and P-values agree.
    assert_eq!(pv(&m, 0.7505), pv(&m, 1.0));
    assert_eq!(get_matches(&m, 0.7505, false).unwrap(), ["AA"]);
}

#[test]
fn motif_score_thresholds() {
    let m = sites(&["TGACTCA", "TGAGTCA"]);
    let (lo, hi) = motif_range(&m, false).unwrap();
    assert!(lo < 0.0 && hi > 0.0);
    let t = [0.0, 0.25, 0.5, 1.0];
    let total = motif_score(&m, &t, ThresholdType::Total, false).unwrap();
    let zero = motif_score(&m, &t, ThresholdType::FromZero, false).unwrap();
    for (i, &f) in t.iter().enumerate() {
        assert_near(total[i], lo + f * (hi - lo), 1e-12, "total");
        assert_near(zero[i], f * hi, 1e-12, "from zero");
    }
    for bad in [-0.1, 1.1, f64::NAN] {
        assert_eq!(
            arg_msg(motif_score(&m, &[0.5, bad], ThresholdType::Total, false)),
            "thresholds must be between 0 and 1"
        );
    }
    for tt in [ThresholdType::Total, ThresholdType::FromZero] {
        assert_eq!(tt.as_str().parse::<ThresholdType>().unwrap(), tt);
        assert_eq!(tt.to_string(), tt.as_str());
    }
    assert_eq!(ThresholdType::default(), ThresholdType::Total);
    assert_eq!(
        arg_msg("zero".parse::<ThresholdType>()),
        "threshold type must be 'total' or 'fromzero'"
    );
}

#[test]
fn nonfinite_weights_when_allowed() {
    // Columns: A or C, anything, T only.
    let m = raw_ppm(
        &[&[0.5, 0.5, 0.0, 0.0], &[0.25; 4], &[0.0, 0.0, 0.0, 1.0]],
        Alphabet::Dna,
    );
    assert_eq!(motif_range(&m, true).unwrap(), (f64::NEG_INFINITY, 3.0));
    // Otherwise a pseudocount is added first.
    let (lo, hi) = motif_range(&m, false).unwrap();
    assert!(lo.is_finite() && hi.is_finite() && lo < hi);

    let s = score_match(&m, &["AAT", "CGT", "GAT", "AAA"], true).unwrap();
    assert_eq!(s, [3.0, 3.0, f64::NEG_INFINITY, f64::NEG_INFINITY]);
    let finite = score_match(&m, &["GAT"], false).unwrap()[0];
    assert!(finite.is_finite() && finite < 3.0);

    let scores = get_scores(&m, true).unwrap();
    assert_eq!(scores.len(), 64);
    assert_eq!(scores.iter().filter(|s| s.is_infinite()).count(), 56);
    assert!(scores[..8].iter().all(|&s| s == 3.0));
    assert!(get_scores(&m, false).unwrap().iter().all(|s| s.is_finite()));

    // Only the eight words avoiding zero probabilities ever match.
    let want: BTreeSet<String> = ['A', 'C']
        .iter()
        .flat_map(|a| DNA.iter().map(move |b| format!("{a}{b}T")))
        .collect();
    for score in [3.0, 0.0, -1000.0] {
        let got: BTreeSet<String> = get_matches(&m, score, true).unwrap().into_iter().collect();
        assert_eq!(got, want, "at {score}");
    }
    assert_eq!(
        motif_score(&m, &[0.0, 1.0], ThresholdType::FromZero, true).unwrap(),
        [0.0, 3.0]
    );
}

#[test]
fn words_are_checked_against_the_motif() {
    let m = sites(&["TGAC"]);
    assert_eq!(
        arg_msg(score_match(&m, &["TGA"], false)),
        "the match has 3 letters but the motif has 4 positions"
    );
    assert_eq!(
        arg_msg(score_match(&m, &["TGAZ"], false)),
        "found letters in 'TGAZ' that are not in the motif's alphabet"
    );
    assert_eq!(
        arg_msg(prob_match(&m, &["TGACT"], false)),
        "the match has 5 letters but the motif has 4 positions"
    );
    // Custom alphabets match case exactly.
    let c = sites_in(&["aBa", "aBB"], Alphabet::custom("aB".chars()).unwrap());
    assert!(score_match(&c, &["aBa"], false).is_ok());
    assert_eq!(
        arg_msg(score_match(&c, &["ABA"], false)),
        "found letters in 'ABA' that are not in the motif's alphabet"
    );
    // Amino acid words match in either case.
    let aa = sites_in(&["MKV"], Alphabet::Aa);
    assert_eq!(
        score_match(&aa, &["mkv"], false).unwrap(),
        score_match(&aa, &["MKV"], false).unwrap()
    );
    // Too many words to list.
    let wide = sites_in(&["MKVLAAGG"], Alphabet::Aa);
    assert_eq!(
        arg_msg(get_scores(&wide, false)),
        "the motif has 25600000000 possible scores, too many to list; use a narrower motif"
    );
}

#[test]
fn word_probabilities() {
    let bkg = vec![0.1, 0.2, 0.3, 0.4];
    let m = Motif::builder(
        Matrix::from_columns(&[[0.25; 4], [0.25; 4]]).unwrap(),
        Alphabet::Dna,
    )
    .bkg(Some(bkg.clone()))
    .build()
    .unwrap();
    let words: Vec<String> = all_words(4, 2)
        .iter()
        .map(|w| word_text(w, &Alphabet::Dna))
        .collect();
    let p = prob_match(&m, &words, false).unwrap();
    assert_near(p.iter().sum(), 1.0, 1e-12, "probabilities sum to 1");
    assert_near(p[0], 0.01, 1e-15, "AA");
    assert_near(p[6], 0.2 * 0.3, 1e-15, "CG");

    // A zero in the background gets the pseudocount spread over the letters,
    // scaled by nsites.
    let zbkg = vec![0.5, 0.5, 0.0, 0.0];
    let z = Motif::builder(
        Matrix::from_columns(&[[0.25; 4], [0.25; 4]]).unwrap(),
        Alphabet::Dna,
    )
    .bkg(Some(zbkg))
    .nsites(Some(10.0))
    .pseudocount(1.0)
    .build()
    .unwrap();
    let add = 0.1 / 4.0;
    let got = prob_match(&z, &["AG", "GT"], false).unwrap();
    assert_near(got[0], (0.5 + add) * add, 1e-15, "AG");
    assert_near(got[1], add * add, 1e-15, "GT");
    assert_eq!(prob_match(&z, &["AG", "AC"], true).unwrap(), [0.0, 0.25]);
    // Without a pseudocount or nsites, 1 and 100 stand in.
    let mut parts = z.into_parts();
    parts.pseudocount = 0.0;
    parts.nsites = None;
    let z0 = Motif::from_parts(parts).unwrap();
    let add = 0.01 / 4.0;
    assert_near(
        prob_match(&z0, &["GG"], false).unwrap()[0],
        add * add,
        1e-15,
        "GG",
    );

    let b = [('A', 0.1), ('C', 0.2), ('G', 0.3), ('T', 0.4)];
    let got = prob_match_bkg(&b, &["ACGT", "TT", ""]).unwrap();
    assert_near(got[0], 0.0024, 1e-15, "ACGT");
    assert_near(got[1], 0.16, 1e-15, "TT");
    assert_eq!(got[2], 1.0);
    assert_eq!(
        arg_msg(prob_match_bkg(&b, &["ACNT"])),
        "the letter 'N' is not in the background"
    );
}

// Column metrics.

const ALL_METRICS: [ColumnMetric; 12] = [
    ColumnMetric::Pcc,
    ColumnMetric::Eucl,
    ColumnMetric::Sw,
    ColumnMetric::Kl,
    ColumnMetric::Weucl,
    ColumnMetric::Allr,
    ColumnMetric::Bhat,
    ColumnMetric::Hell,
    ColumnMetric::Wpcc,
    ColumnMetric::Seucl,
    ColumnMetric::Man,
    ColumnMetric::AllrLl,
];

fn cmp(x: &[f64], y: &[f64], metric: ColumnMetric) -> f64 {
    let u = [0.25; 4];
    compare_columns(
        &Column::new(x).with_bkg(&u),
        &Column::new(y).with_bkg(&u),
        metric,
    )
    .unwrap()
}

#[test]
fn column_metrics_by_hand() {
    let x = [0.5, 0.5, 0.0, 0.0];
    let u = [0.25; 4];
    let half = 0.5f64.sqrt();
    assert_near(cmp(&x, &u, ColumnMetric::Eucl), 0.5, 1e-15, "EUCL");
    assert_near(cmp(&x, &u, ColumnMetric::Seucl), 0.25, 1e-15, "SEUCL");
    assert_near(cmp(&x, &u, ColumnMetric::Sw), 1.75, 1e-15, "SW");
    assert_near(cmp(&x, &u, ColumnMetric::Man), 1.0, 1e-15, "MAN");
    assert_near(cmp(&x, &u, ColumnMetric::Bhat), half, 1e-15, "BHAT");
    assert_near(
        cmp(&x, &u, ColumnMetric::Hell),
        (1.0 - half).sqrt(),
        1e-15,
        "HELL",
    );
    // With uniform weights, WEUCL is EUCL scaled by the root of 2 * 0.25.
    assert_near(cmp(&x, &u, ColumnMetric::Weucl), 0.5 * half, 1e-15, "WEUCL");
    // A zero facing a non-zero probability makes KL infinite, but zeros
    // facing zeros add nothing.
    assert_eq!(cmp(&x, &u, ColumnMetric::Kl), f64::INFINITY);
    let y = [0.5, 0.25, 0.25, 0.0];
    assert_eq!(cmp(&y, &y, ColumnMetric::Kl), 0.0);
    assert_eq!(cmp(&x, &y, ColumnMetric::Kl), f64::INFINITY);
    let a = [0.4, 0.3, 0.2, 0.1];
    let b = [0.1, 0.2, 0.3, 0.4];
    let kl_ab = 0.5
        * (0.4 * 4f64.ln()
            + 0.1 * 0.25f64.ln()
            + 0.3 * 1.5f64.ln()
            + 0.2 * (2.0f64 / 3.0).ln()
            + 0.2 * (2.0f64 / 3.0).ln()
            + 0.3 * 1.5f64.ln()
            + 0.1 * 0.25f64.ln()
            + 0.4 * 4f64.ln());
    assert_near(cmp(&a, &b, ColumnMetric::Kl), kl_ab, 1e-14, "KL");

    // Correlations.
    let p = [0.7, 0.1, 0.1, 0.1];
    assert_near(
        cmp(&x, &p, ColumnMetric::Pcc),
        1.0 / 3f64.sqrt(),
        1e-14,
        "PCC",
    );
    assert_near(
        cmp(&x, &p, ColumnMetric::Wpcc),
        1.0 / 3f64.sqrt(),
        1e-14,
        "WPCC",
    );
    assert_near(cmp(&a, &b, ColumnMetric::Pcc), -1.0, 1e-14, "PCC reversed");
    assert_near(cmp(&a, &a, ColumnMetric::Pcc), 1.0, 1e-14, "PCC self");
    let w = [0.4, 0.4, 0.1, 0.1];
    let wpcc = compare_columns(
        &Column::new(&x).with_bkg(&w),
        &Column::new(&p).with_bkg(&w),
        ColumnMetric::Wpcc,
    )
    .unwrap();
    assert_near(wpcc, 1.0 / 6f64.sqrt(), 1e-14, "weighted WPCC");
    let weucl = compare_columns(
        &Column::new(&x).with_bkg(&w),
        &Column::new(&p).with_bkg(&w),
        ColumnMetric::Weucl,
    )
    .unwrap();
    assert_near(weucl, 0.164f64.sqrt(), 1e-14, "weighted WEUCL");
    // A flat column has no variance; it counts as perfectly correlated.
    assert_eq!(cmp(&x, &u, ColumnMetric::Pcc), 1.0);
    assert_eq!(cmp(&x, &u, ColumnMetric::Wpcc), 1.0);

    // ALLR, with site counts of 10 and 30.
    let allr = |x: &[f64], y: &[f64], metric| {
        compare_columns(
            &Column::new(x).with_bkg(&u).with_nsites(10.0),
            &Column::new(y).with_bkg(&u).with_nsites(30.0),
            metric,
        )
        .unwrap()
    };
    let want = {
        let left: f64 = a
            .iter()
            .zip(&b)
            .map(|(a, b)| b * 30.0 * (a / 0.25).ln())
            .sum();
        let right: f64 = a
            .iter()
            .zip(&b)
            .map(|(a, b)| a * 10.0 * (b / 0.25).ln())
            .sum();
        (left + right) / 40.0
    };
    assert!(want > -2.0);
    assert_near(allr(&a, &b, ColumnMetric::Allr), want, 1e-14, "ALLR");
    assert_near(allr(&a, &b, ColumnMetric::AllrLl), want, 1e-14, "ALLR_LL");
    // A zero probability in a column facing a non-zero one sends ALLR to
    // minus infinity; ALLR_LL stops at -2.
    assert_eq!(allr(&x, &p, ColumnMetric::Allr), f64::NEG_INFINITY);
    assert_eq!(allr(&x, &p, ColumnMetric::AllrLl), -2.0);
    let (s, t) = ([0.97, 0.01, 0.01, 0.01], [0.01, 0.01, 0.01, 0.97]);
    let far = allr(&s, &t, ColumnMetric::Allr);
    assert!(far.is_finite() && far < -2.0, "{far}");
    assert_eq!(allr(&s, &t, ColumnMetric::AllrLl), -2.0);

    // A column whose first value is negative scores 0 with every metric.
    let neg = [-1.0, 0.5, 0.5, 0.0];
    for metric in ALL_METRICS {
        assert_eq!(cmp(&neg, &a, metric), 0.0, "{metric}");
        assert_eq!(cmp(&a, &neg, metric), 0.0, "{metric}");
    }
}

#[test]
fn column_metric_identities() {
    let cols: [[f64; 4]; 4] = [
        [0.4, 0.3, 0.2, 0.1],
        [0.1, 0.2, 0.3, 0.4],
        [0.7, 0.1, 0.1, 0.1],
        [0.05, 0.05, 0.6, 0.3],
    ];
    for metric in ALL_METRICS {
        for x in &cols {
            let own = cmp(x, x, metric);
            let want = match metric {
                ColumnMetric::Pcc | ColumnMetric::Wpcc | ColumnMetric::Bhat => Some(1.0),
                ColumnMetric::Sw => Some(2.0),
                ColumnMetric::Allr | ColumnMetric::AllrLl => None,
                _ => Some(0.0),
            };
            if let Some(w) = want {
                assert_near(own, w, 1e-12, &format!("{metric} of a column with itself"));
            }
            for y in &cols {
                assert_near(
                    cmp(x, y, metric),
                    cmp(y, x, metric),
                    1e-12,
                    &format!("{metric} is symmetric"),
                );
            }
        }
        assert_eq!(metric.as_str().parse::<ColumnMetric>().unwrap(), metric);
        assert_eq!(metric.to_string(), metric.as_str());
    }
    // The Hellinger distance follows from the Bhattacharyya coefficient.
    for x in &cols {
        for y in &cols {
            let h = cmp(x, y, ColumnMetric::Hell);
            let b = cmp(x, y, ColumnMetric::Bhat);
            assert_near(h * h, 1.0 - b, 1e-12, "HELL from BHAT");
        }
    }
    assert_eq!(
        arg_msg("pcc".parse::<ColumnMetric>()),
        "unknown metric 'pcc'; use one of PCC, EUCL, SW, KL, WEUCL, ALLR, BHAT, HELL, WPCC, \
         SEUCL, MAN, ALLR_LL"
    );
}

#[test]
fn column_metric_arguments_are_checked() {
    let a = [0.4, 0.3, 0.2, 0.1];
    let three = [0.5, 0.25, 0.25];
    let one = [1.0];
    let u = [0.25; 4];
    assert_eq!(
        arg_msg(compare_columns(
            &Column::new(&a),
            &Column::new(&three),
            ColumnMetric::Pcc
        )),
        "x and y must have the same length, not 4 and 3"
    );
    assert_eq!(
        arg_msg(compare_columns(
            &Column::new(&one),
            &Column::new(&one),
            ColumnMetric::Eucl
        )),
        "columns should have at least 2 entries"
    );
    for metric in [
        ColumnMetric::Wpcc,
        ColumnMetric::Weucl,
        ColumnMetric::Allr,
        ColumnMetric::AllrLl,
    ] {
        assert_eq!(
            arg_msg(compare_columns(
                &Column::new(&a).with_bkg(&u),
                &Column::new(&a),
                metric
            )),
            "incorrect background vector length"
        );
    }
    // Metrics without a background do not need one.
    assert!(compare_columns(&Column::new(&a), &Column::new(&a), ColumnMetric::Kl).is_ok());
    for metric in [ColumnMetric::Allr, ColumnMetric::AllrLl] {
        assert_eq!(
            arg_msg(compare_columns(
                &Column::new(&a).with_bkg(&u).with_nsites(1.0),
                &Column::new(&a).with_bkg(&u),
                metric
            )),
            "nsites1/nsites2 should be greater than 1"
        );
    }
}

#[test]
fn averages() {
    let v = [4.0, 1.0, -1.0, 0.0];
    assert_eq!(average(&v, Average::ArithmeticMean), 1.0);
    // Only positive values enter the logarithms, but all count.
    assert_near(
        average(&v, Average::GeometricMean),
        2f64.sqrt(),
        1e-15,
        "g.mean",
    );
    assert_eq!(average(&[1.0, -2.0], Average::GeometricMean), 0.0);
    assert_eq!(average(&v, Average::Median), 0.5);
    assert_eq!(average(&[3.0, 1.0, 2.0], Average::Median), 2.0);
    assert_eq!(average(&[7.0], Average::Median), 7.0);
    assert!(average(&[], Average::ArithmeticMean).is_nan());
    assert!(average(&[], Average::FisherZ).is_nan());
    assert_near(average(&[0.5], Average::FisherZ), 0.5, 1e-15, "fzt of one");
    assert_near(
        average(&[0.2, 0.4], Average::FisherZ),
        ((0.2f64.tanh() + 0.4f64.tanh()) / 2.0).atanh(),
        1e-15,
        "fzt",
    );
    for how in [
        Average::ArithmeticMean,
        Average::GeometricMean,
        Average::Median,
        Average::FisherZ,
    ] {
        assert_eq!(how.as_str().parse::<Average>().unwrap(), how);
        assert_eq!(how.to_string(), how.as_str());
    }
    assert_eq!(
        arg_msg("mean".parse::<Average>()),
        "average must be one of a.mean, g.mean, median, fzt"
    );

    // Columns of 2 and 0 bits.
    let m = raw_ppm(&[&[1.0, 0.0, 0.0, 0.0], &[0.25; 4]], Alphabet::Dna);
    let ic = |how| average_ic(from_ref(&m), how).unwrap()[0];
    assert_near(ic(Average::ArithmeticMean), 1.0, 1e-12, "mean IC");
    assert_near(ic(Average::Median), 1.0, 1e-12, "median IC");
    assert_near(ic(Average::GeometricMean), 2f64.sqrt(), 1e-12, "g.mean IC");
}

#[test]
fn median_of_no_values() {
    assert!(average(&[], Average::Median).is_nan());
}

// Statistics.

fn choose(n: u64, k: u64) -> f64 {
    if k > n {
        return 0.0;
    }
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

/// The natural log of n!, summed directly.
fn ln_fact(n: u64) -> f64 {
    (2..=n).map(|k| (k as f64).ln()).sum()
}

fn ln_choose(n: u64, k: u64) -> f64 {
    ln_fact(n) - ln_fact(k) - ln_fact(n - k)
}

#[test]
fn binomial_probabilities() {
    for p in [0.3f64, 0.8, 0.5] {
        let n = 20;
        let pmf: Vec<f64> = (0..=n)
            .map(|x| choose(n, x) * p.powi(x as i32) * (1.0 - p).powi((n - x) as i32))
            .collect();
        for x in 0..=n {
            assert_near(dbinom(x, n, p), pmf[x as usize], 1e-12, "dbinom");
        }
        for k in 0..=n + 1 {
            let want: f64 = pmf.iter().skip(k as usize).sum();
            let got = pbinom_upper(k, n, p);
            assert!(
                (got - want).abs() < 1e-14 + 1e-12 * want,
                "P(X >= {k}), p = {p}: {got} vs {want}"
            );
        }
    }
    assert_eq!(dbinom(0, 0, 0.4), 1.0);
    assert_eq!(dbinom(0, 5, 0.0), 1.0);
    assert_eq!(dbinom(1, 5, 0.0), 0.0);
    assert_eq!(dbinom(5, 5, 1.0), 1.0);
    assert_eq!(dbinom(4, 5, 1.0), 0.0);
    assert_eq!(dbinom(6, 5, 0.5), 0.0);
    assert_eq!(pbinom_upper(0, 5, 0.0), 1.0);
    assert_eq!(pbinom_upper(1, 5, 0.0), 0.0);
    assert_eq!(pbinom_upper(3, 5, 1.0), 1.0);
    assert_eq!(pbinom_upper(6, 5, 0.5), 0.0);
}

#[test]
fn hypergeometric_probabilities() {
    let (r, b, n) = (7u64, 5u64, 6u64);
    let pmf = |x: u64| choose(r, x) * choose(b, n - x.min(n)) / choose(r + b, n);
    for x in 0..=8u64 {
        let want = if x > n { 0.0 } else { pmf(x) };
        assert_near(dhyper(x, r, b, n), want, 1e-12, &format!("dhyper {x}"));
    }
    assert_eq!(dhyper(0, 3, 4, 0), 1.0);
    assert_eq!(dhyper(4, 3, 4, 5), 0.0);
    assert_eq!(dhyper(0, 3, 1, 3), 0.0);
    for x in -2..=8i64 {
        let lower: f64 = (0..=x.max(-1))
            .filter(|&i| i >= 0 && i as u64 <= n)
            .map(|i| pmf(i as u64))
            .sum();
        assert_near(
            phyper(x, r, b, n, true),
            lower,
            1e-12,
            &format!("lower {x}"),
        );
        assert_near(
            phyper(x, r, b, n, false),
            1.0 - lower,
            1e-12,
            &format!("upper {x}"),
        );
    }
    assert!(phyper(1, 2, 2, 5, true).is_nan());
    // Every ball drawn, or no black ones to draw.
    assert_eq!(phyper(3, 3, 4, 7, true), 1.0);
    assert_eq!(phyper(3, 3, 4, 7, false), 0.0);
    assert_eq!(phyper(5, 7, 0, 5, true), 1.0);

    // [[3, 1], [1, 3]]: three or more of four in the top-left cell.
    assert_near(fisher_greater(3, 1, 1, 3), 17.0 / 70.0, 1e-14, "fisher");
    assert_near(
        ln_fisher_greater(3, 1, 1, 3),
        (17.0f64 / 70.0).ln(),
        1e-14,
        "ln fisher",
    );
    assert_eq!(fisher_greater(0, 4, 4, 0), 1.0);

    // Too small for a double: a perfect split of 1400 items, and a near one.
    let got = ln_fisher_greater(700, 0, 0, 700);
    assert!(fisher_greater(700, 0, 0, 700) < 1e-300);
    assert_near(got, -ln_choose(1400, 700), 1e-10, "perfect split");
    let (a, b2, c, d) = (900u64, 5u64, 5u64, 900u64);
    assert!(fisher_greater(a, b2, c, d) < 1e-300);
    let (rr, nn, kk) = (a + c, a + b2 + c + d, a + b2);
    let terms: Vec<f64> = (a..=rr.min(kk))
        .map(|x| ln_choose(rr, x) + ln_choose(nn - rr, kk - x) - ln_choose(nn, kk))
        .collect();
    let top = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let want = top + terms.iter().map(|t| (t - top).exp()).sum::<f64>().ln();
    assert_near(ln_fisher_greater(a, b2, c, d), want, 1e-10, "near split");
}

#[test]
fn pvalue_adjustment() {
    let p = [0.01, 0.04, 0.03, 0.005, 0.5, 0.04];
    let bh = p_adjust_bh(&p);
    let n = p.len() as f64;
    // Each q-value is the smallest p_(j) n / j over ranks j at or above its own.
    let mut sorted: Vec<f64> = p.to_vec();
    sorted.sort_by(f64::total_cmp);
    for (i, &pi) in p.iter().enumerate() {
        let rank = sorted.iter().rposition(|&s| s == pi).unwrap();
        let want = (rank..sorted.len())
            .map(|j| sorted[j] * n / (j + 1) as f64)
            .fold(1.0, f64::min);
        assert_near(bh[i], want, 1e-15, &format!("BH {i}"));
    }
    let ln_p: Vec<f64> = p.iter().map(|x| x.ln()).collect();
    for (q, b) in ln_p_adjust_bh(&ln_p, n).iter().zip(&bh) {
        assert_near(*q, b.ln(), 1e-14, "log BH");
    }
    // Counting tests that were not given scales q-values up, capped at 1.
    let m = 100.0;
    let lq = ln_p_adjust_bh(&ln_p, m);
    for (i, &pi) in p.iter().enumerate() {
        let rank = sorted.iter().rposition(|&s| s == pi).unwrap();
        let want = (rank..sorted.len())
            .map(|j| sorted[j] * m / (j + 1) as f64)
            .fold(1.0, f64::min);
        assert_near(lq[i], want.ln(), 1e-14, &format!("BH of {m} tests"));
    }
    let bonf = p_adjust_bonferroni(&[0.01, 0.3, 0.5]);
    assert_near(bonf[0], 0.03, 1e-15, "Bonferroni");
    assert_near(bonf[1], 0.9, 1e-15, "Bonferroni");
    assert_eq!(bonf[2], 1.0);
    assert!(p_adjust_bh(&[]).is_empty());
    assert_eq!(median(&[]), None);
    assert_eq!(median(&[3.0, 1.0, 2.0]), Some(2.0));
    assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), Some(2.5));
}

// Comparisons.

fn compare_set() -> Vec<Motif> {
    vec![
        sites(&["TGACTCA", "TGAGTCA"]),
        sites(&["TGACGTCA"]),
        sites(&["CACGTGAC"]),
        sites(&["GATTACAGG", "GATAACAGG"]),
        sites(&["AAAACCCC"]),
    ]
}

#[test]
fn comparison_arguments_are_checked() {
    let set = compare_set();
    let opts = CompareOptions::default();
    assert_eq!(
        arg_msg(compare_motifs(&[], &[], &opts)),
        "no motifs to compare"
    );
    let mixed = [set[0].clone(), sites_in(&["MKVL"], Alphabet::Aa)];
    assert_eq!(
        arg_msg(compare_motifs(&mixed, &[0], &opts)),
        "all motifs must share the same alphabet"
    );
    let zero = CompareOptions {
        min_overlap: 0,
        ..Default::default()
    };
    assert_eq!(
        arg_msg(compare_motifs(&set, &[0], &zero)),
        "min_overlap must be at least 1"
    );
    assert_eq!(
        arg_msg(compare_motifs(&set, &[5], &opts)),
        "motif index out of range"
    );
    assert_eq!(
        arg_msg(compare_matrix(&set[..1], &opts, true)),
        "need at least 2 motifs for an all-against-all comparison"
    );
    assert_eq!(
        arg_msg(anchor_layout(&set, 5, &opts)),
        "anchor index out of range"
    );
    let aa = [
        sites_in(&["MKVLAG"], Alphabet::Aa),
        sites_in(&["MKILAG"], Alphabet::Aa),
    ];
    let parametric = CompareOptions {
        null: NullModel::Parametric,
        ..Default::default()
    };
    match compare_motifs(&aa, &[0], &parametric) {
        Err(Error::Unsupported(m)) => assert_eq!(
            m,
            "the parametric null is only available for DNA and RNA motifs; use the empirical \
             null for AA motifs"
        ),
        other => panic!("expected an unsupported error, got {other:?}"),
    }
    assert!(compare_motifs(&set, &[], &opts).unwrap().is_empty());
}

#[test]
fn self_comparisons_and_reverse_complements() {
    let set = compare_set();
    for (null, complete, rc) in [
        (NullModel::Empirical, false, true),
        (NullModel::Empirical, true, true),
        (NullModel::Empirical, true, false),
        (NullModel::Parametric, false, true),
        (NullModel::Parametric, true, true),
    ] {
        {
            let opts = CompareOptions {
                null,
                complete,
                rc,
                ..Default::default()
            };
            let ctx = format!("{null:?}, complete {complete}, rc {rc}");
            let rows = compare_motifs(&set, &[0, 3], &opts).unwrap();
            assert_eq!(rows.len(), 2 * set.len());
            for r in &rows {
                assert!((0.0..=1.0).contains(&r.pvalue), "{ctx}: {r:?}");
                assert!((-1.0..=1.0 + 1e-12).contains(&r.score), "{ctx}: {r:?}");
                assert!(r.qvalue >= r.pvalue, "{ctx}: {r:?}");
            }
            for q in [0usize, 3] {
                let own: Vec<_> = rows.iter().filter(|r| r.query == q).collect();
                let s = &own[q];
                assert_eq!(s.target, q);
                assert_near(s.score, 1.0, 1e-12, &ctx);
                assert_eq!(
                    (s.offset, s.strand, s.overlap),
                    (0, Strand::Plus, set[q].width())
                );
                assert_eq!(s.query_consensus, s.target_consensus);
                assert_eq!(Some(&s.query_consensus), set[q].consensus.as_ref());
                // The motif itself is the best match of its query.
                for o in &own {
                    assert!(o.pvalue >= s.pvalue, "{ctx}: {o:?}");
                }
                // q-values are the BH adjustment of the query's P-values.
                let p: Vec<f64> = own.iter().map(|r| r.pvalue).collect();
                let q_want = p_adjust_bh(&p);
                for (r, w) in own.iter().zip(q_want) {
                    assert_eq!(r.qvalue, w);
                }
            }
        }
    }

    // A motif's reverse complement aligns on the minus strand, fully.
    let m = set[3].clone();
    let rc = reverse_complement(&m, false).unwrap();
    let pair = [m.clone(), rc.clone()];
    let row = &compare_motifs(&pair, &[0], &CompareOptions::default()).unwrap()[1];
    assert_eq!((row.strand, row.overlap, row.offset), (Strand::Minus, 9, 0));
    assert_near(row.score, 1.0, 1e-12, "reverse complement");
    assert_eq!(Some(&row.target_consensus), rc.consensus.as_ref());
    assert_eq!(Some(&row.query_consensus), m.consensus.as_ref());
    // A piece of the reverse complement: the query's consensus is still read
    // in its own orientation, over the columns it shares with the target.
    let piece = sites(&[&rc_word("TTACAG", false)]);
    let pair = [m.clone(), piece.clone()];
    let row = &compare_motifs(&pair, &[0], &CompareOptions::default()).unwrap()[1];
    assert_eq!((row.strand, row.overlap), (Strand::Minus, 6));
    assert_eq!(row.query_consensus, m.consensus.as_ref().unwrap()[2..8]);
    assert_eq!(Some(&row.target_consensus), piece.consensus.as_ref());
    let pair = [m.clone(), rc.clone()];
    let fwd_only = CompareOptions {
        rc: false,
        ..Default::default()
    };
    let row = &compare_motifs(&pair, &[0], &fwd_only).unwrap()[1];
    assert_eq!(row.strand, Strand::Plus);
    assert!(row.score < 0.9);
}

#[test]
fn comparison_matrices_are_symmetric() {
    let dna = compare_set();
    let aa = vec![
        sites_in(&["MKVLAG", "MKILAG"], Alphabet::Aa),
        sites_in(&["WWKVLA"], Alphabet::Aa),
        sites_in(&["PPGGSTY"], Alphabet::Aa),
    ];
    let custom_alph = Alphabet::custom("xyz".chars()).unwrap();
    let custom = vec![
        sites_in(&["xxyzzyx", "xxyzzyy"], custom_alph.clone()),
        sites_in(&["zzyxxyz"], custom_alph.clone()),
        sites_in(&["xyzxyzx"], custom_alph),
    ];
    for (name, set) in [("DNA", dna), ("AA", aa), ("custom", custom)] {
        for rc in [true, false] {
            let opts = CompareOptions {
                rc,
                min_overlap: 4,
                ..Default::default()
            };
            let full = compare_matrix(&set, &opts, true).unwrap();
            let fast = compare_matrix(&set, &opts, false).unwrap();
            let n = set.len();
            for i in 0..n {
                assert_near(full.scores[i][i], 1.0, 1e-12, name);
                for j in 0..n {
                    let ctx = format!("{name}, rc {rc}, [{i}][{j}]");
                    assert_near(full.scores[i][j], full.scores[j][i], 1e-12, &ctx);
                    assert_eq!(full.scores[i][j], fast.scores[i][j], "{ctx}");
                    assert!(fast.pvalues[i][j].is_nan() && fast.qvalues[i][j].is_nan());
                    assert!((0.0..=1.0).contains(&full.pvalues[i][j]), "{ctx}");
                }
            }
            if !set[0].alphabet.is_nucleotide() {
                let rows = compare_motifs(&set, &[0], &opts).unwrap();
                assert!(rows.iter().all(|r| r.strand == Strand::Plus));
            }
        }
    }
}

#[test]
fn anchored_layouts() {
    let a = sites(&["GATTACAGG"]);
    let b = sites(&["CCGATTACAGG"]);
    // The reverse complement of the anchor with one more column.
    let c = sites(&["CCTGTAATCT"]);
    let d = sites(&["CGC"]);
    let motifs = [a, b, c, d];
    let layout = anchor_layout(&motifs, 0, &CompareOptions::default()).unwrap();
    let got: Vec<(usize, bool)> = layout.placements.iter().map(|p| (p.left, p.rc)).collect();
    // The anchor starts two columns in, after B's extra columns; C is
    // flipped and starts one column in; D, narrower than the minimum
    // overlap, is aligned over its whole width, as comparisons allow.
    assert_eq!(got, [(2, false), (0, false), (1, true), (5, false)]);
    assert_eq!(layout.width, 11);

    let alone = anchor_layout(&motifs[..1], 0, &CompareOptions::default()).unwrap();
    assert_eq!((alone.placements[0].left, alone.width), (0, 9));
}

// Comparison by Euclidean distance and with background overhangs.

/// PPM columns of a motif.
fn ppm_columns(m: &Motif) -> Vec<Vec<f64>> {
    convert_type(m, MotifType::Ppm, None, false)
        .unwrap()
        .matrix
        .columns()
        .map(<[f64]>::to_vec)
        .collect()
}

fn rc_columns(cols: &[Vec<f64>]) -> Vec<Vec<f64>> {
    cols.iter()
        .rev()
        .map(|c| c.iter().rev().copied().collect())
        .collect()
}

/// A column score as comparisons use it: the Pearson correlation of columns
/// centred on 1/n, 0 if either is uniform, or the negated Euclidean
/// distance.
fn column_score(a: &[f64], b: &[f64], eucl: bool) -> f64 {
    if eucl {
        return -compare_columns(&Column::new(a), &Column::new(b), ColumnMetric::Eucl).unwrap();
    }
    let centre = 1.0 / a.len() as f64;
    let ca: Vec<f64> = a.iter().map(|x| x - centre).collect();
    let cb: Vec<f64> = b.iter().map(|x| x - centre).collect();
    let ss = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>();
    if ss(&ca) < 1e-12 || ss(&cb) < 1e-12 {
        return 0.0;
    }
    ca.iter().zip(&cb).map(|(x, y)| x * y).sum::<f64>() / (ss(&ca) * ss(&cb)).sqrt()
}

/// Every alignment of a target of width `wt` to a query of width `wq`
/// overlapping at least `mo` columns, as (offset, query start, target
/// start, overlap).
fn layouts(wq: usize, wt: usize, mo: usize) -> Vec<(i64, usize, usize, usize)> {
    let (wq, wt) = (wq as i64, wt as i64);
    let mo = (mo as i64).min(wq).min(wt);
    (mo - wt..=wq - mo)
        .filter_map(|d| {
            let t_start = (-d).max(0);
            let len = (wq - d).min(wt) - t_start;
            (len >= mo).then_some((d, (d + t_start) as usize, t_start as usize, len as usize))
        })
        .collect()
}

/// The terms an alignment's score is made of, by target column: each
/// target column facing the query is scored against the query column it
/// faces, counted twice with background overhangs, and each one hanging
/// over is scored against the background, with background overhangs only.
/// Returns, for each target column, the query column it faces if any, and
/// the sum of the query's overhanging columns' background scores.
#[allow(clippy::type_complexity)]
fn alignment_terms(
    q: &[Vec<f64>],
    wt: usize,
    (q_start, t_start, len): (usize, usize, usize),
    bkg: Option<&[f64]>,
    eucl: bool,
) -> (Vec<Option<usize>>, f64) {
    let facing = (0..wt)
        .map(|j| {
            (t_start..t_start + len)
                .contains(&j)
                .then(|| q_start + j - t_start)
        })
        .collect();
    let q_out = bkg.map_or(0.0, |b| {
        (0..q.len())
            .filter(|i| !(q_start..q_start + len).contains(i))
            .map(|i| column_score(&q[i], b, eucl))
            .sum()
    });
    (facing, q_out)
}

/// The score of target column `t` in an alignment, with its weight.
fn target_term(
    q: &[Vec<f64>],
    t: &[f64],
    facing: Option<usize>,
    bkg: Option<&[f64]>,
    eucl: bool,
) -> Option<(f64, usize)> {
    match (facing, bkg) {
        (Some(i), None) => Some((column_score(&q[i], t, eucl), 1)),
        (Some(i), Some(_)) => Some((2.0 * column_score(&q[i], t, eucl), 2)),
        (None, Some(b)) => Some((column_score(t, b, eucl), 1)),
        (None, None) => None,
    }
}

const UNIFORM: [f64; 4] = [0.25; 4];

#[test]
fn alignments_by_method_and_overhang() {
    let mut set = compare_set();
    // Informative columns around a shared core, and nearly uniform ones.
    set.push(sites(&["CCGATTACAGGTT"]));
    set.push(raw_ppm(
        &[
            &[0.3, 0.2, 0.25, 0.25],
            &[0.9, 0.0, 0.1, 0.0],
            &[0.0, 0.0, 0.0, 1.0],
            &[0.0, 0.0, 0.0, 1.0],
            &[0.25, 0.25, 0.3, 0.2],
        ],
        Alphabet::Dna,
    ));
    let cols: Vec<Vec<Vec<f64>>> = set.iter().map(ppm_columns).collect();
    let skewed = [0.3, 0.2, 0.2, 0.3];
    for (eucl, overhang, bkg) in [
        (false, Overhang::Ignore, UNIFORM),
        (true, Overhang::Ignore, UNIFORM),
        (false, Overhang::Background, UNIFORM),
        (true, Overhang::Background, UNIFORM),
        (false, Overhang::Background, skewed),
        (true, Overhang::Background, skewed),
    ] {
        let opts = CompareOptions {
            method: if eucl {
                ColumnMetric::Eucl
            } else {
                ColumnMetric::Pcc
            },
            overhang,
            bkg,
            min_overlap: 3,
            ..Default::default()
        };
        let ctx = format!("EUCL {eucl}, {overhang:?}, bkg {bkg:?}");
        let b = (overhang == Overhang::Background).then_some(&bkg[..]);
        let queries: Vec<usize> = (0..set.len()).collect();
        let rows = compare_motifs(&set, &queries, &opts).unwrap();
        let matrix = compare_matrix(&set, &opts, false).unwrap();
        for r in &rows {
            let ctx = format!("{ctx}, {} vs {}", r.query, r.target);
            assert_eq!(matrix.scores[r.query][r.target], r.score, "{ctx}");
            let (wq, wt) = (cols[r.query].len(), cols[r.target].len());
            // What each alignment is chosen by and the score it reports.
            let judge = |q: &[Vec<f64>], (q_start, t_start, len): (usize, usize, usize)| {
                let (facing, q_out) = alignment_terms(q, wt, (q_start, t_start, len), b, eucl);
                let total: f64 = (0..wt)
                    .filter_map(|j| target_term(q, &cols[r.target][j], facing[j], b, eucl))
                    .map(|x| x.0)
                    .sum::<f64>()
                    + q_out;
                let (chosen, mean) = match overhang {
                    Overhang::Background => (total, total / (wq + wt) as f64),
                    _ if eucl => (total / len as f64, total / len as f64),
                    _ => (total, total / len as f64),
                };
                (chosen, if eucl { -mean } else { mean })
            };
            let plus = cols[r.query].clone();
            let minus = rc_columns(&plus);
            let mut best = f64::NEG_INFINITY;
            let mut own = None;
            for (q, strand) in [(&plus, Strand::Plus), (&minus, Strand::Minus)] {
                for (d, q_start, t_start, len) in layouts(wq, wt, 3) {
                    let (chosen, mean) = judge(q, (q_start, t_start, len));
                    best = best.max(chosen);
                    if (d, strand, len) == (r.offset, r.strand, r.overlap) {
                        own = Some((chosen, mean));
                    }
                }
            }
            let (chosen, mean) = own.unwrap_or_else(|| panic!("{ctx}: no such alignment"));
            assert_near(chosen, best, 1e-12, &ctx);
            assert_near(r.score, mean, 1e-12, &ctx);
            if r.query == r.target {
                assert_eq!(
                    (r.offset, r.strand, r.overlap),
                    (0, Strand::Plus, wq),
                    "{ctx}"
                );
                if eucl {
                    assert_eq!(r.score.to_bits(), 0.0f64.to_bits(), "{ctx}");
                }
            }
            assert!((0.0..=1.0).contains(&r.pvalue), "{ctx}: {r:?}");
        }
    }
}

#[test]
fn background_overhangs_with_correlations() {
    let set = compare_set();
    let plain = compare_motifs(&set, &[0, 3], &CompareOptions::default()).unwrap();
    let padded = compare_motifs(
        &set,
        &[0, 3],
        &CompareOptions {
            overhang: Overhang::Background,
            ..Default::default()
        },
    )
    .unwrap();
    // A uniform background correlates 0 with anything, which leaves the
    // alignments as they are and averages the score over both motifs.
    for (a, b) in plain.iter().zip(&padded) {
        assert_eq!(
            (a.offset, a.strand, a.overlap),
            (b.offset, b.strand, b.overlap)
        );
        let n = (set[a.query].width() + set[a.target].width()) as f64;
        assert_near(b.score, a.score * 2.0 * a.overlap as f64 / n, 1e-12, "PCC");
        assert!(b.pvalue <= (a.pvalue * 1.5).min(1.0) && b.pvalue >= a.pvalue / 1.5);
    }
}

#[test]
fn euclidean_distances_penalize_informative_overhangs() {
    let core = sites(&["GATTACA"]);
    let flanked = sites(&["CCGATTACAGG"]);
    let flat = raw_ppm(
        &[
            &[0.25, 0.25, 0.25, 0.25],
            &[0.25, 0.25, 0.25, 0.25],
            &[0.9, 0.0, 0.1, 0.0],
            &[0.0, 0.0, 0.0, 1.0],
            &[0.0, 0.0, 0.0, 1.0],
            &[0.25, 0.25, 0.25, 0.25],
        ],
        Alphabet::Dna,
    );
    let set = [core, flanked, flat];
    let score = |overhang, t: usize| {
        let opts = CompareOptions {
            method: ColumnMetric::Eucl,
            overhang,
            ..Default::default()
        };
        compare_motifs(&set, &[0], &opts).unwrap()[t].score
    };
    // Overlaps alone match the core perfectly either way.
    assert!(score(Overhang::Ignore, 1) < 0.05);
    // Uniform overhangs cost nothing, informative ones cost their distance
    // from the background.
    let flanked_bkg = score(Overhang::Background, 1);
    assert!(
        flanked_bkg > 3.0 * score(Overhang::Ignore, 1),
        "{flanked_bkg}"
    );
    let cols = ppm_columns(&set[1]);
    let outside: f64 = [0, 1, 9, 10]
        .iter()
        .map(|&j| -column_score(&cols[j], &UNIFORM, true))
        .sum();
    assert!(flanked_bkg >= outside / 18.0 && flanked_bkg < outside / 18.0 + 0.05);
}

#[test]
fn comparison_options_are_checked() {
    let set = compare_set();
    let both = CompareOptions {
        complete: true,
        overhang: Overhang::Background,
        ..Default::default()
    };
    assert_eq!(
        arg_msg(compare_motifs(&set, &[0], &both)),
        "complete scores and background overhangs cannot be combined"
    );
    let bad_bkg = CompareOptions {
        overhang: Overhang::Background,
        bkg: [0.5, 0.5, 0.5, 0.5],
        ..Default::default()
    };
    assert_eq!(
        arg_msg(compare_matrix(&set, &bad_bkg, false)),
        "bkg must be four non-negative values summing to 1"
    );
    let sw = CompareOptions {
        method: ColumnMetric::Sw,
        ..Default::default()
    };
    match anchor_layout(&set, 0, &sw) {
        Err(Error::Unsupported(m)) => {
            assert_eq!(m, "motifs can only be compared by PCC or EUCL, not SW")
        }
        other => panic!("expected an unsupported error, got {other:?}"),
    }
    // Amino acids are compared with a uniform background.
    let aa = [
        sites_in(&["MKVLAG"], Alphabet::Aa),
        sites_in(&["KVLAGW"], Alphabet::Aa),
    ];
    let rows = compare_motifs(&aa, &[0], &bad_bkg).unwrap();
    assert_eq!((rows[1].offset, rows[1].overlap), (1, 5));
}

/// The P-value of every alignment against the exact distribution of its
/// score under the empirical null, where each target column is drawn from
/// the columns of all the motifs. Column scores are rounded to 51 bins
/// over their range for correlations and to 201 for distances, as the
/// comparisons do, which the P-values must match exactly; without rounding,
/// the P-values must lie between the exact probabilities of the score plus
/// and minus the most the rounding can move it.
#[test]
fn comparison_pvalues_match_enumeration() {
    let set = vec![
        sites(&["TGACTCA", "TGAGTCA"]),
        sites(&["TGACGTC"]),
        sites(&["CACGTGA", "CACATGA", "CACGTGT"]),
        sites(&["GATTAC"]),
        sites(&["AACCC", "AACCG"]),
    ];
    let pool: Vec<Vec<f64>> = set.iter().flat_map(ppm_columns).collect();
    let cols: Vec<Vec<Vec<f64>>> = set.iter().map(ppm_columns).collect();
    let skewed = [0.35, 0.15, 0.15, 0.35];
    let mut checked = 0;
    for (eucl, overhang, bkg) in [
        (false, Overhang::Ignore, UNIFORM),
        (true, Overhang::Ignore, UNIFORM),
        (false, Overhang::Background, skewed),
        (true, Overhang::Background, UNIFORM),
        (true, Overhang::Background, skewed),
    ] {
        let opts = CompareOptions {
            method: if eucl {
                ColumnMetric::Eucl
            } else {
                ColumnMetric::Pcc
            },
            overhang,
            bkg,
            min_overlap: 4,
            ..Default::default()
        };
        let (lo, hi, steps) = if eucl {
            (-std::f64::consts::SQRT_2, 0.0, 200)
        } else {
            (-1.0, 1.0, 50)
        };
        let width = (hi - lo) / steps as f64;
        let bin = |s: f64| ((s - lo) / width + 0.5).floor() as usize;
        let b = (overhang == Overhang::Background).then_some(&bkg[..]);
        for q in 0..set.len() {
            let pairs: Vec<(usize, usize)> = (0..set.len()).map(|t| (q, t)).collect();
            let al = universalmotif::compare::align(&set, &pairs, &opts).unwrap();
            let rows = compare_motifs(&set, &[q], &opts).unwrap();
            for (a, r) in al.iter().zip(&rows) {
                let ctx = format!("EUCL {eucl}, {overhang:?}, {bkg:?}, {q} vs {}", a.target);
                let qq = if a.strand == Strand::Minus {
                    rc_columns(&cols[q])
                } else {
                    cols[q].clone()
                };
                let wt = cols[a.target].len();
                let t_start = (-a.offset).max(0) as usize;
                let (facing, _) =
                    alignment_terms(&qq, wt, (a.query_start, t_start, a.overlap), b, eucl);
                let mut observed = 0.0;
                let mut units = 0;
                // Binned and exact distributions of the target's part.
                let mut binned = vec![1.0];
                let mut exact: Vec<(f64, f64)> = vec![(0.0, 1.0)];
                for (j, &f) in facing.iter().enumerate() {
                    let Some((s, w)) = target_term(&qq, &cols[a.target][j], f, b, eucl) else {
                        continue;
                    };
                    observed += s;
                    units += w;
                    // The distinct scores of a pool column here, with their
                    // bins and probabilities.
                    let mut terms: Vec<(f64, usize, f64)> = Vec::new();
                    for t in &pool {
                        let (s, w) = target_term(&qq, t, f, b, eucl).unwrap();
                        let p = 1.0 / pool.len() as f64;
                        match terms.iter_mut().find(|x| (x.0 - s).abs() < 1e-12) {
                            Some(x) => x.2 += p,
                            None => terms.push((s, bin(s / w as f64) * w, p)),
                        }
                    }
                    let mut next = vec![0.0; binned.len() + steps * w];
                    for (k, &p) in binned.iter().enumerate() {
                        for &(_, idx, pt) in &terms {
                            next[k + idx] += p * pt;
                        }
                    }
                    binned = next;
                    let mut merged: BTreeMap<i64, (f64, f64)> = BTreeMap::new();
                    for &(v, p) in &exact {
                        for &(s, _, pt) in &terms {
                            let e = merged.entry(((v + s) * 1e9).round() as i64).or_default();
                            *e = (v + s, e.1 + p * pt);
                        }
                    }
                    exact = merged.into_values().collect();
                }
                let top = binned.iter().rposition(|&p| p > 0.0).unwrap();
                let per_col = ((observed / units as f64 - lo) / (hi - lo)).clamp(0.0, 1.0);
                let k = ((per_col * steps as f64 * units as f64 + 0.5).floor() as usize).min(top);
                let want: f64 = binned[k..].iter().sum();
                let got = r.pvalue / a.n_tested as f64;
                if r.pvalue >= 1.0 {
                    assert!(want * a.n_tested as f64 >= 1.0 - 1e-9, "{ctx}: {want}");
                    continue;
                }
                assert_near(got, want, 1e-9, &ctx);
                let reach = exact.last().unwrap().0;
                let slack = (units + 1) as f64 * width / 2.0 + 1e-9;
                let at_least = |x: f64| exact.iter().filter(|e| e.0 >= x).map(|e| e.1).sum::<f64>();
                let (low, high) = (
                    at_least(observed + slack),
                    at_least(observed.min(reach) - slack),
                );
                assert!(
                    low - 1e-12 <= got && got <= high + 1e-12,
                    "{ctx}: {low} {got} {high}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 60, "{checked}");
}

// Clustering.

/// Four points on a line at 0, 1, 3 and 7.
fn line_distances() -> Vec<Vec<f64>> {
    let x = [0.0f64, 1.0, 3.0, 7.0];
    x.iter()
        .map(|a| x.iter().map(|b| (a - b).abs()).collect())
        .collect()
}

const ALL_LINKAGES: [Linkage; 8] = [
    Linkage::WardD,
    Linkage::Single,
    Linkage::Complete,
    Linkage::Average,
    Linkage::Mcquitty,
    Linkage::Median,
    Linkage::Centroid,
    Linkage::WardD2,
];

#[test]
fn linkage_heights_by_hand() {
    let d = line_distances();
    let want: [(Linkage, [f64; 3]); 8] = [
        (Linkage::Single, [1.0, 2.0, 4.0]),
        (Linkage::Complete, [1.0, 3.0, 7.0]),
        (Linkage::Average, [1.0, 2.5, 17.0 / 3.0]),
        (Linkage::Mcquitty, [1.0, 2.5, 5.25]),
        (Linkage::Centroid, [1.0, 2.25, 5.0]),
        (Linkage::Median, [1.0, 2.25, 4.5625]),
        (Linkage::WardD, [1.0, 3.0, 7.5]),
        (
            Linkage::WardD2,
            [1.0, (25.0f64 / 3.0).sqrt(), (578.0f64 / 12.0).sqrt()],
        ),
    ];
    for (method, heights) in want {
        let t = hclust(&d, method).unwrap();
        for (g, w) in t.height.iter().zip(heights) {
            assert_near(*g, w, 1e-12, method.as_str());
        }
        // Every method joins 1 and 2, then 3, then 4.
        assert_eq!(t.merge, [[-1, -2], [-3, 1], [-4, 2]], "{method}");
        assert_eq!(t.order, [3, 2, 0, 1], "{method}");
        assert_eq!(
            t.linkage(),
            [
                [0.0, 1.0, t.height[0], 2.0],
                [2.0, 4.0, t.height[1], 3.0],
                [3.0, 5.0, t.height[2], 4.0],
            ]
        );
    }
    for m in ALL_LINKAGES {
        assert_eq!(m.as_str().parse::<Linkage>().unwrap(), m);
        assert_eq!(m.to_string(), m.as_str());
    }
    assert_eq!(Linkage::default(), Linkage::Complete);
    assert_eq!(
        arg_msg("ward".parse::<Linkage>()),
        "unknown clustering method 'ward' (use ward.D, single, complete, average, mcquitty, \
         median, centroid or ward.D2)"
    );
}

#[test]
fn cutting_trees() {
    let t = hclust(&line_distances(), Linkage::Complete).unwrap();
    assert_eq!(t.cutree(1).unwrap(), [0, 0, 0, 0]);
    assert_eq!(t.cutree(2).unwrap(), [0, 0, 0, 1]);
    assert_eq!(t.cutree(3).unwrap(), [0, 0, 1, 2]);
    assert_eq!(t.cutree(4).unwrap(), [0, 1, 2, 3]);
    for k in [0, 5] {
        assert_eq!(arg_msg(t.cutree(k)), "k must be between 1 and 4");
    }
    // Cutting at a height keeps the merges at or below it.
    for (h, k) in [
        (0.5, 4),
        (1.0, 3),
        (2.9, 3),
        (3.0, 2),
        (6.99, 2),
        (7.0, 1),
        (100.0, 1),
    ] {
        assert_eq!(t.cutree_height(h).unwrap(), t.cutree(k).unwrap(), "h = {h}");
    }
    // Groups are numbered by their first leaf.
    let shuffled = Tree {
        merge: vec![[-2, -4], [-1, -3], [1, 2]],
        height: vec![1.0, 2.0, 3.0],
        order: vec![1, 3, 0, 2],
    };
    assert_eq!(shuffled.cutree(2).unwrap(), [0, 1, 0, 1]);
    assert_eq!(shuffled.cutree(3).unwrap(), [0, 1, 2, 1]);

    // Centroid linkage can merge below an earlier merge.
    let tri = vec![
        vec![0.0, 1.0, 1.1],
        vec![1.0, 0.0, 1.1],
        vec![1.1, 1.1, 0.0],
    ];
    let inv = hclust(&tri, Linkage::Centroid).unwrap();
    assert_near(inv.height[1], 0.85, 1e-12, "inversion");
    assert_eq!(
        arg_msg(inv.cutree_height(1.0)),
        "the tree heights are not sorted increasingly"
    );
    assert_eq!(inv.cutree(2).unwrap(), [0, 0, 1]);
}

#[test]
fn hclust_inputs_and_newick() {
    assert_eq!(
        arg_msg(hclust(&[vec![0.0]], Linkage::Single)),
        "need at least 2 objects to cluster"
    );
    assert_eq!(
        arg_msg(hclust(&[vec![0.0, 1.0], vec![1.0]], Linkage::Single)),
        "the distance matrix must be square"
    );
    assert_eq!(
        arg_msg(hclust(
            &[vec![0.0, 1.0], vec![f64::INFINITY, 0.0]],
            Linkage::Single
        )),
        "distances must be finite"
    );
    // Only the lower triangle is read.
    let t = hclust(&[vec![0.0, f64::NAN], vec![2.0, 0.0]], Linkage::Average).unwrap();
    assert_eq!(t.height, [2.0]);

    let t = hclust(&line_distances(), Linkage::Complete).unwrap();
    let labels: Vec<String> = ["a b", "c,d", "(e)", "f:g"].map(String::from).to_vec();
    assert_eq!(
        t.newick(&labels),
        "(f-g:3.5,(e:1.5,(a_b:0.5,c-d:0.5):1):2);"
    );
    // Leaves without labels are numbered from 1.
    assert_eq!(
        t.newick(&labels[..2]),
        "(4:3.5,(3:1.5,(a_b:0.5,c-d:0.5):1):2);"
    );
    let big = hclust(&[vec![0.0, 0.0], vec![3e9, 0.0]], Linkage::Single).unwrap();
    assert_eq!(big.newick(&[]), "(1:1500000000,2:1500000000);");
    let tiny = hclust(&[vec![0.0, 0.0], vec![3e-7, 0.0]], Linkage::Single).unwrap();
    assert_eq!(tiny.newick(&[]), "(1:1.5e-07,2:1.5e-07);");
}

#[test]
fn motif_trees_follow_the_comparisons() {
    let motifs = vec![
        sites(&["TGACTCAT", "TGAGTCAT"]),
        sites(&["CACGTGCC"]),
        sites(&["TGACTCAG"]),
        sites(&["CACGTGCA", "CACATGCA"]),
    ];
    let copts = CompareOptions {
        min_overlap: 6,
        ..Default::default()
    };
    let scores = compare_matrix(&motifs, &copts, false).unwrap().scores;
    let dist: Vec<Vec<f64>> = scores
        .iter()
        .map(|r| r.iter().map(|s| (1.0 - s) / 2.0).collect())
        .collect();
    for method in ALL_LINKAGES {
        let opts = TreeOptions {
            method,
            nthreads: 2,
            ..Default::default()
        };
        let tree = motif_tree(&motifs, &opts).unwrap();
        assert_eq!(tree, hclust(&dist, method).unwrap(), "{method}");
        // The two families split apart.
        assert_eq!(tree.cutree(2).unwrap(), [0, 1, 0, 1], "{method}");
        assert_eq!(tree.order.len(), 4);
        let mut leaves = tree.order.clone();
        leaves.sort();
        assert_eq!(leaves, [0, 1, 2, 3]);
    }
    let pair = motif_tree(&motifs[..2], &TreeOptions::default()).unwrap();
    assert_near(pair.height[0], dist[0][1], 1e-12, "two motifs");
    assert!(motif_tree(&motifs[..1], &TreeOptions::default()).is_err());
}

// Scanning.

/// Check every hit against `score_match` and `motif_pvalue`, and that every
/// window reaching the threshold is a hit. Returns the number of hits.
fn check_scan(motifs: &[Motif], seqs: &[String], opts: &ScanOptions) -> usize {
    let table = scan_table(motifs, seqs, opts).unwrap();
    let alphabet = &motifs[0].alphabet;
    let rna = *alphabet == Alphabet::Rna;
    let exact = matches!(alphabet, Alphabet::Custom(_));
    let both = alphabet.is_nucleotide() && opts.rc;
    let mut got = Vec::new();
    for i in 0..table.len() {
        let m = &motifs[table.motif[i] as usize];
        let (si, start, end) = (
            table.sequence[i] as usize,
            table.start[i] as usize,
            table.end[i] as usize,
        );
        let word = table.matched(i, seqs);
        let fwd = &seqs[si][start..end];
        let want_word = if table.minus[i] {
            rc_word(fwd, rna)
        } else if exact {
            fwd.to_string()
        } else {
            fwd.to_ascii_uppercase()
        };
        assert_eq!(word, want_word);
        assert_eq!(end - start, m.width());
        let score = score_match(m, &[&word], false).unwrap()[0];
        assert_eq!(table.score[i], score, "{word}");
        let p = pv(m, score);
        assert_eq!(table.pvalue[i], p, "{word}");
        assert!(p <= opts.pvalue);
        let pwm = convert_type(m, MotifType::Pwm, None, false).unwrap();
        let best: f64 = pwm
            .matrix
            .columns()
            .map(|c| c.iter().copied().fold(f64::NEG_INFINITY, f64::max))
            .sum();
        assert_near(table.score_pct[i], 100.0 * score / best, 1e-12, "score_pct");
        got.push((table.motif[i] as usize, si, start, table.minus[i]));
    }
    let mut sorted = got.clone();
    sorted.sort();
    assert_eq!(got, sorted, "hits are ordered");

    let mut want = Vec::new();
    for (mi, m) in motifs.iter().enumerate() {
        let thr = threshold(m, opts.pvalue);
        let w = m.width();
        let known = |c: char| {
            alphabet.letters().iter().any(|&l| {
                if exact {
                    l == c
                } else {
                    l.eq_ignore_ascii_case(&c)
                }
            })
        };
        for (si, s) in seqs.iter().enumerate() {
            for start in 0..(s.len() + 1).saturating_sub(w) {
                let word = &s[start..start + w];
                if !word.chars().all(known) {
                    continue;
                }
                if score_match(m, &[word], false).unwrap()[0] >= thr {
                    want.push((mi, si, start, false));
                }
                if both && score_match(m, &[rc_word(word, rna)], false).unwrap()[0] >= thr {
                    want.push((mi, si, start, true));
                }
            }
        }
    }
    want.sort();
    assert_eq!(got, want, "hits are every window at or above the threshold");
    table.len()
}

fn random_seqs(letters: &[char], n: usize, len: usize, seed: u64) -> Vec<String> {
    create_sequences(letters, n, len, None, seed, 1).unwrap()
}

#[test]
fn dna_scans_agree_with_match_scores() {
    let motifs = vec![
        sites(&["TGACTCA", "TGAGTCA", "TGACTAA"]),
        sites(&["CACGTG"]),
        sites(&["AAAA"]),
    ];
    let mut seqs = random_seqs(&DNA, 12, 150, 11);
    seqs[0] = seqs[0].to_lowercase();
    seqs[1].replace_range(20..24, "NNNN");
    seqs[2].replace_range(5..12, "TGACTCA");
    seqs.push("TGA".into());
    seqs.push(String::new());
    for rc in [true, false] {
        for p in [0.05, 1e-3] {
            let opts = ScanOptions {
                pvalue: p,
                rc,
                nthreads: 2,
                ..Default::default()
            };
            let n = check_scan(&motifs, &seqs, &opts);
            assert!(n > 0, "rc {rc}, P {p}");
        }
    }
    // RNA reports U for A on the minus strand.
    let rna: Vec<String> = seqs.iter().map(|s| s.replace(['T', 't'], "U")).collect();
    let rm = vec![sites_in(&["UGACUCA", "UGAGUCA"], Alphabet::Rna)];
    let opts = ScanOptions {
        pvalue: 0.05,
        ..Default::default()
    };
    assert!(check_scan(&rm, &rna, &opts) > 0);
}

#[test]
fn other_alphabet_scans_agree_with_match_scores() {
    let opts = ScanOptions {
        pvalue: 0.01,
        ..Default::default()
    };
    // Amino acids, in either case and with unknown letters.
    let aa = vec![
        sites_in(&["MKVLA", "MRVLA", "MKILA"], Alphabet::Aa),
        sites_in(&["WW"], Alphabet::Aa),
    ];
    let mut seqs = random_seqs(&AMINO_ACIDS, 6, 300, 5);
    seqs[0] = seqs[0].to_lowercase();
    seqs[1].replace_range(30..33, "XXX");
    assert!(check_scan(&aa, &seqs, &opts) > 0);

    // Custom letters match case exactly.
    let alph = Alphabet::custom("ABa".chars()).unwrap();
    let c = vec![sites_in(&["aAB", "aAB", "aBB"], alph)];
    let seqs = vec!["aABaABxaab".to_string(), "AABAAB".into(), "aBBaAB".into()];
    let opts3 = ScanOptions {
        pvalue: 0.05,
        ..Default::default()
    };
    // Only aAB is likely enough, at 0 and 3 of the first and 3 of the last.
    assert_eq!(check_scan(&c, &seqs, &opts3), 3);

    // More than 64 letters are scored one position per lookup.
    let letters: Vec<char> = (0..70u8).map(|i| (b'0' + i) as char).collect();
    let big = Alphabet::custom(letters.iter().copied()).unwrap();
    let ropts = CreateOptions {
        nsites: Some(50.0),
        pseudocount: 1.0,
        ..Default::default()
    };
    let motifs = vec![
        random(3, big.clone(), &ropts, 7).unwrap(),
        random(2, big, &ropts, 8).unwrap(),
    ];
    let seqs = random_seqs(&letters, 4, 500, 9);
    assert!(check_scan(&motifs, &seqs, &opts) > 0);
}

#[test]
fn scan_arguments_are_checked() {
    let m = sites(&["TGACTCA"]);
    let seqs = ["TGACTCAGG"];
    let opts = ScanOptions::default();
    assert_eq!(
        arg_msg(scan_table(&[], &seqs, &opts)),
        "no motifs to scan with"
    );
    for p in [0.0, 1.0, f64::NAN, -0.5] {
        let o = ScanOptions {
            pvalue: p,
            ..Default::default()
        };
        assert_eq!(
            arg_msg(scan_table(from_ref(&m), &seqs, &o)),
            "the P-value cutoff must be between 0 and 1"
        );
    }
    let mixed = [m.clone(), sites_in(&["MKV"], Alphabet::Aa)];
    assert_eq!(
        arg_msg(scan_sequences(&mixed, &seqs, &opts)),
        "all motifs must share the same alphabet"
    );
    let rna = EncodedSequences::new(&["UGACUCA"], &Alphabet::Rna, 1).unwrap();
    let msg = "the sequences were encoded as RNA but the motifs are DNA";
    assert_eq!(arg_msg(scan_encoded(from_ref(&m), &rna, &opts)), msg);
    assert_eq!(
        arg_msg(motif_matches(from_ref(&m), &rna, &opts, MatchValue::Counts)),
        msg
    );
    assert_eq!(
        arg_msg(motif_matches(&[], &rna, &opts, MatchValue::Counts)),
        "no motifs to scan with"
    );
    for v in [MatchValue::Matches, MatchValue::Counts, MatchValue::Scores] {
        assert_eq!(v.as_str().parse::<MatchValue>().unwrap(), v);
        assert_eq!(v.to_string(), v.as_str());
    }
    assert_eq!(
        arg_msg("count".parse::<MatchValue>()),
        "unknown value 'count'; use one of matches, counts and scores"
    );
    // A cutoff below the motif's smallest P-value gives no hits.
    let tiny = ScanOptions {
        pvalue: 1e-12,
        ..Default::default()
    };
    assert!(scan_table(from_ref(&m), &seqs, &tiny).unwrap().is_empty());
}

#[test]
fn overlapping_hits_are_removed() {
    let motifs = vec![sites(&["TGACTCA", "TGAGTCA"]), sites(&["GACTC"])];
    let seqs = random_seqs(&DNA, 20, 200, 21);
    let base = ScanOptions {
        pvalue: 0.02,
        ..Default::default()
    };
    let all = scan_sequences(&motifs, &seqs, &base).unwrap();
    for by in [OverlapPriority::Score, OverlapPriority::Pvalue] {
        for (ignore_strand, ignore_motif) in [(false, false), (true, false), (true, true)] {
            let opts = ScanOptions {
                no_overlaps: Some(NoOverlaps {
                    by,
                    ignore_strand,
                    ignore_motif,
                }),
                ..base.clone()
            };
            let kept = scan_sequences(&motifs, &seqs, &opts).unwrap();
            assert!(!kept.is_empty() && kept.len() < all.len());
            let same_group = |a: &Hit, b: &Hit| {
                a.sequence == b.sequence
                    && (ignore_motif || a.motif == b.motif)
                    && (ignore_strand || a.strand == b.strand)
            };
            let overlap = |a: &Hit, b: &Hit| same_group(a, b) && a.start < b.end && b.start < a.end;
            let better = |a: &Hit, b: &Hit| match by {
                OverlapPriority::Score => a.score >= b.score,
                OverlapPriority::Pvalue => a.pvalue <= b.pvalue,
            };
            for (i, a) in kept.iter().enumerate() {
                assert!(all.contains(a));
                for b in &kept[i + 1..] {
                    assert!(!overlap(a, b), "{a:?} overlaps {b:?}");
                }
            }
            // Each dropped hit overlaps a kept one at least as good.
            for d in all.iter().filter(|h| !kept.contains(h)) {
                assert!(
                    kept.iter().any(|k| overlap(k, d) && better(k, d)),
                    "{d:?} was dropped for nothing"
                );
            }
        }
    }
}

#[test]
fn fasta_scans_match_the_table() {
    let motifs = vec![sites(&["TGACTCA", "TGAGTCA"]), sites(&["CACGTG"])];
    let seqs = random_seqs(&DNA, 5, 300, 31);
    let dir = std::env::temp_dir();
    let path = dir.join(format!("scoring_edge_{}.fa", std::process::id()));
    let text: String = seqs
        .iter()
        .enumerate()
        .map(|(i, s)| format!(">seq{i}\n{s}\n"))
        .collect();
    std::fs::write(&path, text).unwrap();
    let opts = ScanOptions {
        pvalue: 1e-2,
        ..Default::default()
    };
    let fasta = scan_fasta(&motifs, &path, &opts);
    let none = scan_fasta(
        from_ref(&motifs[0]),
        &path,
        &ScanOptions {
            pvalue: 1e-12,
            ..Default::default()
        },
    );
    std::fs::remove_file(&path).unwrap();
    let fasta = fasta.unwrap();
    let table = scan_table(&motifs, &seqs, &opts).unwrap();
    assert!(!fasta.is_empty());
    assert_eq!(fasta.len(), table.len());
    assert_eq!(fasta.table, table);
    assert_eq!(fasta.names, ["seq0", "seq1", "seq2", "seq3", "seq4"]);
    assert_eq!(fasta.lengths, [300; 5]);
    for i in 0..table.len() {
        assert_eq!(fasta.matched(i), table.matched(i, &seqs));
    }
    let none = none.unwrap();
    assert!(none.is_empty());
    assert_eq!(none.len(), 0);
    assert_eq!(none.names.len(), 5);
}

#[test]
fn gapped_scans_agree_with_match_scores() {
    let core = sites(&["TTGACGTCAA", "TTGACGTCTA"]);
    let gapped = add_gap(&core, Some(&[5]), &[1], &[3]).unwrap();
    let flat = ungap(&gapped, true);
    let mut seqs = random_seqs(&DNA, 8, 160, 41);
    seqs[0].replace_range(10..22, "TTGACAAGTCAA");
    seqs[1].replace_range(40..50, "TTGACGTCAA");
    for correct_gaps in [false, true] {
        let opts = ScanOptions {
            pvalue: 1e-3,
            correct_gaps,
            ..Default::default()
        };
        let t = scan_table(from_ref(&gapped), &seqs, &opts).unwrap();
        assert!(!t.is_empty());
        let combos = if correct_gaps { 3.0 } else { 1.0 };
        let thr = threshold(&flat, opts.pvalue / combos);
        let mut got = Vec::new();
        for i in 0..t.len() {
            let gaps = t.gap_sizes(i).unwrap();
            assert_eq!(gaps.len(), 1);
            let g = gaps[0];
            assert!((1..=3).contains(&g));
            assert_eq!((t.end[i] - t.start[i]) as usize, 10 + g);
            let matched = t.matched(i, &seqs);
            assert_eq!(&matched[5..5 + g], &".".repeat(g));
            let letters = matched.replace('.', "");
            let score = score_match(&flat, &[&letters], false).unwrap()[0];
            assert_eq!(t.score[i], score, "{matched}");
            assert_eq!(t.pvalue[i], (pv(&flat, score) * combos).min(1.0));
            got.push((t.start[i] as usize, t.minus[i], g));
        }
        let mut want = Vec::new();
        for (si, s) in seqs.iter().enumerate() {
            for g in 1..=3 {
                for start in 0..(s.len() + 1).saturating_sub(10 + g) {
                    let win = &s[start..start + 10 + g];
                    for minus in [false, true] {
                        let o = if minus {
                            rc_word(win, false)
                        } else {
                            win.to_string()
                        };
                        let letters = format!("{}{}", &o[..5], &o[5 + g..]);
                        if score_match(&flat, &[&letters], false).unwrap()[0] >= thr {
                            want.push((si, start, minus, g));
                        }
                    }
                }
            }
        }
        let mut got_full: Vec<(usize, usize, bool, usize)> = (0..t.len())
            .map(|i| (t.sequence[i] as usize, got[i].0, got[i].1, got[i].2))
            .collect();
        got_full.sort();
        want.sort();
        assert_eq!(got_full, want, "correct_gaps {correct_gaps}");
    }
    // Without gaps, the motif scans as its ungapped self.
    let off = ScanOptions {
        use_gaps: false,
        ..Default::default()
    };
    let a = scan_table(from_ref(&gapped), &seqs, &off).unwrap();
    let b = scan_table(from_ref(&flat), &seqs, &off).unwrap();
    assert_eq!(a.start, b.start);
    assert_eq!(a.score, b.score);
    assert!(!a.is_empty());
    assert_eq!(a.gap_sizes(0), None);
}

#[test]
fn unusable_gaps_are_refused() {
    let m = sites(&["ACGTACGTAC"]);
    let seqs = ["ACGTACGTACGTACGT"];
    let opts = ScanOptions::default();
    let backwards = add_gap(&m, Some(&[6, 3]), &[1], &[2]).unwrap();
    assert_eq!(
        arg_msg(scan_table(from_ref(&backwards), &seqs, &opts)),
        "the gap locations of motif 'motif' must be in increasing order"
    );
    let many = add_gap(&m, Some(&[2, 4, 6]), &[0], &[40]).unwrap();
    assert_eq!(
        arg_msg(scan_table(from_ref(&many), &seqs, &opts)),
        "motif 'motif' has 68921 combinations of gap sizes, more than 65536 can be scanned"
    );
    // Just within the limit: 16 * 16 * 16 * 16 combinations.
    let most = add_gap(&m, Some(&[2, 4, 6, 8]), &[0], &[15]).unwrap();
    assert!(scan_table(from_ref(&most), &seqs, &opts).is_ok());
}
