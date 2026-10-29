//! Randomized property tests: invariants checked over many seeded random
//! motifs and sequences, to catch what hand-picked cases miss.
//!
//! Every property runs over a range of seeds. A failure names the property
//! and the seed, and setting `PROPERTY_SEED` to that seed reruns just that
//! case. Inputs are printed with the failure.

// Brute-force checks read most clearly with indices.
#![allow(clippy::needless_range_loop)]

use std::collections::{BTreeMap, HashMap};
use std::panic::AssertUnwindSafe;

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use universalmotif::cluster::{Linkage, TreeOptions, hclust, motif_tree};
use universalmotif::compare::{CompareOptions, compare_matrix};
use universalmotif::complexity::{ComplexityMethod, complexity};
use universalmotif::convert::{convert_type, position_icscore};
use universalmotif::io::{cisbp, homer, jaspar, matrix as matrix_io, meme, transfac, uniprobe};
use universalmotif::kmer::{BkgOptions, count_kmers, get_bkg, kmers, sequence_letters};
use universalmotif::mask::{mask_ranges, mask_seqs};
use universalmotif::merge::{MergeOptions, merge_motifs};
use universalmotif::ops::{
    TrimFrom, add_gap, normalize, reverse_complement, round_motif, switch_alph, trim_motif, ungap,
};
use universalmotif::packed::EncodedSequences;
use universalmotif::pvalue::{PvalueOptions, min_pvalues, motif_pvalue, pvalue_to_score};
use universalmotif::scan::{
    HitTable, NoOverlaps, OverlapPriority, ScanOptions, scan_encoded, scan_table,
};
use universalmotif::scoring::{motif_range, score_match};
use universalmotif::shuffle::{ShuffleMethod, ShuffleOptions, shuffle_motifs, shuffle_sequences};
use universalmotif::stats;
use universalmotif::{Alphabet, GapInfo, Matrix, Motif, MotifType};

// ---------------------------------------------------------------------------
// Harness

type Check = Result<(), String>;

macro_rules! ensure {
    ($cond:expr, $($arg:tt)*) => {
        {
            let ok: bool = $cond;
            if !ok {
                return Err(format!($($arg)*));
            }
        }
    };
}

trait Ctx<T> {
    fn ctx(self, what: &str) -> Result<T, String>;
}

impl<T, E: std::fmt::Display> Ctx<T> for Result<T, E> {
    fn ctx(self, what: &str) -> Result<T, String> {
        self.map_err(|e| format!("{what}: {e}"))
    }
}

/// Run a property over `cases` seeds, or only over `PROPERTY_SEED` if set.
/// Each property gets its own stream of seeds, derived from its name.
fn run(name: &str, cases: u64, f: impl Fn(&mut StdRng) -> Check) {
    let base = name.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ b as u64).wrapping_mul(0x0100_0000_01b3)
    });
    let seeds: Vec<u64> = match std::env::var("PROPERTY_SEED") {
        Ok(s) => vec![s.parse().expect("PROPERTY_SEED must be a number")],
        Err(_) => (0..cases).collect(),
    };
    for seed in seeds {
        let mut rng = StdRng::seed_from_u64(base ^ seed);
        let rerun = format!("rerun with PROPERTY_SEED={seed}");
        match std::panic::catch_unwind(AssertUnwindSafe(|| f(&mut rng))) {
            Ok(Ok(())) => {}
            Ok(Err(e)) => panic!("property '{name}' failed for seed {seed} ({rerun}):\n{e}"),
            Err(p) => {
                let msg = p
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default();
                panic!("property '{name}' panicked for seed {seed} ({rerun}): {msg}");
            }
        }
    }
}

fn close(a: f64, b: f64, abs: f64, rel: f64) -> bool {
    if a == b {
        return true;
    }
    (a - b).abs() <= abs + rel * a.abs().max(b.abs())
}

fn desc(m: &Motif) -> String {
    let cols: Vec<Vec<f64>> = m.matrix.columns().map(|c| c.to_vec()).collect();
    format!(
        "motif {:?}: {} {} width {} pseudocount {} nsites {:?} bkg {:?} gap {:?}\ncolumns {:?}",
        m.name,
        m.alphabet,
        m.motif_type,
        m.width(),
        m.pseudocount,
        m.nsites,
        m.bkg,
        m.gap,
        cols
    )
}

fn descs(ms: &[Motif]) -> String {
    ms.iter().map(desc).collect::<Vec<_>>().join("\n")
}

// ---------------------------------------------------------------------------
// Generators

const POOL: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

fn pick<T: Clone>(rng: &mut StdRng, v: &[T]) -> T {
    v[rng.random_range(0..v.len())].clone()
}

fn custom_alphabet(rng: &mut StdRng, lo: usize, hi: usize, pool: &[u8]) -> Alphabet {
    let k = rng.random_range(lo..=hi);
    let mut idx: Vec<usize> = (0..pool.len()).collect();
    for i in 0..k {
        let j = rng.random_range(i..idx.len());
        idx.swap(i, j);
    }
    Alphabet::custom(idx[..k].iter().map(|&i| pool[i] as char)).unwrap()
}

fn gen_alphabet(rng: &mut StdRng) -> Alphabet {
    match rng.random_range(0..10) {
        0..=3 => Alphabet::Dna,
        4..=5 => Alphabet::Rna,
        6..=7 => Alphabet::Aa,
        _ => custom_alphabet(rng, 2, 8, POOL),
    }
}

/// Background frequencies: uniform a third of the time, otherwise skewed but
/// with no letter above 0.9.
fn gen_bkg(rng: &mut StdRng, n: usize) -> Vec<f64> {
    if rng.random_bool(0.35) {
        return vec![1.0 / n as f64; n];
    }
    let raw: Vec<f64> = (0..n).map(|_| rng.random_range(0.2..1.0)).collect();
    let s: f64 = raw.iter().sum();
    raw.iter().map(|v| v / s).collect()
}

/// A probability column: uniform, one-hot, sparse with zeros, or skewed.
fn gen_ppm_col(rng: &mut StdRng, n: usize) -> Vec<f64> {
    let kind = rng.random_range(0..12);
    let mut v: Vec<f64> = match kind {
        0 => vec![1.0; n],
        1 => {
            let mut v = vec![0.0; n];
            v[rng.random_range(0..n)] = 1.0;
            v
        }
        2..=4 => {
            let mut v: Vec<f64> = (0..n)
                .map(|_| {
                    if rng.random_bool(0.4) {
                        0.0
                    } else {
                        rng.random::<f64>().powi(2) + 1e-3
                    }
                })
                .collect();
            if v.iter().all(|&x| x == 0.0) {
                v[rng.random_range(0..n)] = 1.0;
            }
            v
        }
        _ => {
            let p = rng.random_range(1..5);
            (0..n).map(|_| rng.random::<f64>().powi(p) + 1e-4).collect()
        }
    };
    let s: f64 = v.iter().sum();
    v.iter_mut().for_each(|x| *x /= s);
    v
}

fn gen_ppm(rng: &mut StdRng, n: usize, w: usize) -> Matrix {
    let cols: Vec<Vec<f64>> = (0..w).map(|_| gen_ppm_col(rng, n)).collect();
    Matrix::from_columns(&cols).unwrap()
}

fn gen_nsites(rng: &mut StdRng) -> Option<f64> {
    match rng.random_range(0..8) {
        0 | 1 => None,
        2 => Some(pick(rng, &[1.0, 2.0, 3.0, 7.0])),
        _ => Some(rng.random_range(1..=200) as f64),
    }
}

fn gen_pseudo(rng: &mut StdRng) -> f64 {
    pick(rng, &[0.0, 0.0, 0.01, 0.5, 1.0, 1.0, 3.7])
}

/// Gaps for a motif of width `w`: one or two, at increasing locations.
fn gen_gaps(rng: &mut StdRng, w: usize) -> GapInfo {
    if w < 2 {
        return GapInfo::default();
    }
    let ngaps = rng.random_range(1..=(w - 1).min(2));
    let mut locs: Vec<usize> = Vec::new();
    while locs.len() < ngaps {
        let l = rng.random_range(1..w);
        if !locs.contains(&l) {
            locs.push(l);
        }
    }
    locs.sort_unstable();
    let mingap: Vec<usize> = (0..ngaps).map(|_| rng.random_range(0..=2)).collect();
    let maxgap: Vec<usize> = mingap
        .iter()
        .map(|&lo| lo.max(1) + rng.random_range(0..=2))
        .collect();
    GapInfo {
        isgapped: true,
        gaploc: locs,
        mingap,
        maxgap,
    }
}

/// Probabilities to whole counts summing to `total`, by largest remainder.
fn counts(col: &[f64], total: f64) -> Vec<f64> {
    let scaled: Vec<f64> = col.iter().map(|p| p * total).collect();
    let mut out: Vec<f64> = scaled.iter().map(|v| v.floor()).collect();
    let mut short = (total - out.iter().sum::<f64>()).round() as usize;
    let mut order: Vec<usize> = (0..col.len()).collect();
    order.sort_by(|&a, &b| (scaled[b] - out[b]).total_cmp(&(scaled[a] - out[a])));
    for &i in order.iter().cycle() {
        if short == 0 {
            break;
        }
        out[i] += 1.0;
        short -= 1;
    }
    out
}

#[derive(Clone, Copy)]
struct Gen {
    /// Allow raw PWM, ICM and CWM matrices not derived from probabilities.
    raw: bool,
    /// Probability of gaps.
    gaps: f64,
}

const PLAIN: Gen = Gen {
    raw: true,
    gaps: 0.0,
};

/// A random motif, avoiding ICMs with a column of zeros, which no
/// conversion accepts (see `bug_icm_zero_column_round_trip`).
fn gen_motif(rng: &mut StdRng, alph: &Alphabet, w: usize, ty: MotifType, g: Gen) -> Motif {
    loop {
        if let Some(m) = gen_motif_any(rng, alph, w, ty, g)
            && !zero_icm_column(&m)
        {
            return m;
        }
    }
}

/// A random motif, or `None` if converting it to an ICM failed (see
/// `bug_uniform_column_negative_ic`).
fn gen_motif_any(
    rng: &mut StdRng,
    alph: &Alphabet,
    w: usize,
    ty: MotifType,
    g: Gen,
) -> Option<Motif> {
    let n = alph.len();
    let ppm = gen_ppm(rng, n, w);
    let bkg = gen_bkg(rng, n);
    let pseudo = gen_pseudo(rng);
    let nsites = gen_nsites(rng);
    let gap = if rng.random_bool(g.gaps) {
        gen_gaps(rng, w)
    } else {
        GapInfo::default()
    };
    let name = format!("m{}", rng.random_range(0..100_000));
    let build = |m: Matrix, t: MotifType, ns: Option<f64>| {
        Motif::builder(m, alph.clone())
            .name(name.clone())
            .motif_type(Some(t))
            .bkg(Some(bkg.clone()))
            .pseudocount(pseudo)
            .nsites(ns)
            .gap(gap.clone())
            .build()
            .unwrap()
    };
    let from_ppm = |t: MotifType, rel: bool| {
        let m = build(ppm.clone(), MotifType::Ppm, nsites);
        convert_type(&m, t, None, rel).ok()
    };
    let raw = g.raw && rng.random_bool(0.4);
    Some(match ty {
        MotifType::Ppm => build(ppm.clone(), MotifType::Ppm, nsites),
        MotifType::Pcm => {
            let total = nsites.unwrap_or_else(|| rng.random_range(1..=60) as f64);
            let cols: Vec<Vec<f64>> = ppm.columns().map(|c| counts(c, total)).collect();
            build(Matrix::from_columns(&cols).unwrap(), MotifType::Pcm, nsites)
        }
        MotifType::Pwm if raw => {
            let v: Vec<f64> = (0..n * w).map(|_| rng.random_range(-6.0..3.0)).collect();
            build(Matrix::from_vec(n, w, v).unwrap(), MotifType::Pwm, nsites)
        }
        MotifType::Pwm => from_ppm(MotifType::Pwm, false)?,
        MotifType::Icm if raw => {
            let cols: Vec<Vec<f64>> = (0..w)
                .map(|_| {
                    let mut c: Vec<f64> = (0..n)
                        .map(|_| {
                            if rng.random_bool(0.3) {
                                0.0
                            } else {
                                rng.random_range(0.0..2.0)
                            }
                        })
                        .collect();
                    if c.iter().all(|&x| x == 0.0) {
                        c[0] = 0.5;
                    }
                    c
                })
                .collect();
            build(Matrix::from_columns(&cols).unwrap(), MotifType::Icm, nsites)
        }
        MotifType::Icm => from_ppm(MotifType::Icm, rng.random_bool(0.3))?,
        _ => {
            let v: Vec<f64> = (0..n * w)
                .map(|_| {
                    if rng.random_bool(0.2) {
                        0.0
                    } else {
                        rng.random_range(-1.0..1.0)
                    }
                })
                .collect();
            build(Matrix::from_vec(n, w, v).unwrap(), MotifType::Cwm, nsites)
        }
    })
}

fn gen_type(rng: &mut StdRng, cwm: bool) -> MotifType {
    let n = if cwm { 5 } else { 4 };
    MotifType::ALL[rng.random_range(0..n)]
}

/// A letter outside the alphabet.
fn unknown_letter(rng: &mut StdRng, alph: &Alphabet) -> char {
    let opts: Vec<char> = match alph {
        Alphabet::Dna => vec!['N', 'n', 'U', '-'],
        Alphabet::Rna => vec!['N', 'n', 'T', '-'],
        Alphabet::Aa => vec!['X', 'x', '*', '-'],
        Alphabet::Custom(l) => POOL
            .iter()
            .map(|&b| b as char)
            .chain(['-', '.'])
            .filter(|c| !l.contains(c))
            .collect(),
    };
    pick(rng, &opts)
}

/// A random sequence of the alphabet's letters, in either case unless the
/// alphabet is custom, with some letters outside it.
fn gen_seq(rng: &mut StdRng, alph: &Alphabet, len: usize, unknown: f64) -> String {
    let letters = alph.letters();
    let custom = matches!(alph, Alphabet::Custom(_));
    let lower = rng.random_bool(0.5);
    (0..len)
        .map(|_| {
            if rng.random_bool(unknown) {
                unknown_letter(rng, alph)
            } else {
                let c = pick(rng, letters);
                if !custom && lower && rng.random_bool(0.3) {
                    c.to_ascii_lowercase()
                } else {
                    c
                }
            }
        })
        .collect()
}

fn gen_seqs(rng: &mut StdRng, alph: &Alphabet, count: usize, max_len: usize) -> Vec<String> {
    let unknown = pick(rng, &[0.0, 0.02, 0.1]);
    (0..count)
        .map(|_| {
            let len = match rng.random_range(0..6) {
                0 => rng.random_range(0..4),
                _ => rng.random_range(0..=max_len),
            };
            gen_seq(rng, alph, len, unknown)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Brute-force helpers

/// The row of a sequence letter: exact for custom alphabets, in either case
/// otherwise.
fn letter_index(alph: &Alphabet, c: char) -> Option<usize> {
    match alph {
        Alphabet::Custom(l) => l.iter().position(|&x| x == c),
        _ => alph
            .letters()
            .iter()
            .position(|x| x.eq_ignore_ascii_case(&c)),
    }
}

fn complement(c: char, rna: bool) -> char {
    match c.to_ascii_uppercase() {
        'A' => {
            if rna {
                'U'
            } else {
                'T'
            }
        }
        'T' | 'U' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

fn rc_chars(s: &[char], rna: bool) -> Vec<char> {
    s.iter().rev().map(|&c| complement(c, rna)).collect()
}

/// Every word of `k` letters, first letter varying slowest.
fn all_words(letters: &[char], k: usize) -> Vec<String> {
    let mut out = vec![String::new()];
    for _ in 0..k {
        out = out
            .iter()
            .flat_map(|p| letters.iter().map(move |l| format!("{p}{l}")))
            .collect();
    }
    out
}

/// One way of placing a gapped motif: the offset of each motif column in the
/// window, and the window's width.
struct Placement {
    offsets: Vec<usize>,
    span: usize,
}

fn placements(m: &Motif) -> Vec<Placement> {
    let w = m.width();
    let g = &m.gap;
    if !g.isgapped || g.gaploc.is_empty() {
        return vec![Placement {
            offsets: (0..w).collect(),
            span: w,
        }];
    }
    let mut bounds = vec![0];
    bounds.extend(&g.gaploc);
    bounds.push(w);
    let mut combos: Vec<Vec<usize>> = vec![vec![]];
    for (&lo, &hi) in g.mingap.iter().zip(&g.maxgap) {
        combos = combos
            .iter()
            .flat_map(|c| {
                (lo..=hi).map(move |s| {
                    let mut c = c.clone();
                    c.push(s);
                    c
                })
            })
            .collect();
    }
    combos
        .into_iter()
        .map(|sizes| {
            let mut offsets = Vec::with_capacity(w);
            let mut at = 0;
            for (k, b) in bounds.windows(2).enumerate() {
                for _ in b[0]..b[1] {
                    offsets.push(at);
                    at += 1;
                }
                if k < sizes.len() {
                    at += sizes[k];
                }
            }
            Placement { offsets, span: at }
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
struct HitKey {
    motif: usize,
    seq: usize,
    start: usize,
    end: usize,
    minus: bool,
    matched: String,
    score: f64,
}

fn sort_keys(v: &mut [HitKey]) {
    v.sort_by(|a, b| {
        (a.motif, a.seq, a.start, a.end, a.minus, &a.matched)
            .cmp(&(b.motif, b.seq, b.start, b.end, b.minus, &b.matched))
    });
}

fn table_keys(t: &HitTable, seqs: &[String]) -> Vec<HitKey> {
    let mut v: Vec<HitKey> = (0..t.len())
        .map(|i| HitKey {
            motif: t.motif[i] as usize,
            seq: t.sequence[i] as usize,
            start: t.start[i] as usize,
            end: t.end[i] as usize,
            minus: t.minus[i],
            matched: t.matched(i, seqs),
            score: t.score[i],
        })
        .collect();
    sort_keys(&mut v);
    v
}

/// Every window of every sequence scoring at least the motif's threshold
/// for `cutoff`, found by scoring each window with `score_match()`.
fn brute_hits(
    mi: usize,
    m: &Motif,
    seqs: &[String],
    cutoff: f64,
    rc: bool,
    use_gaps: bool,
) -> Result<Vec<HitKey>, String> {
    let thr = pvalue_to_score(
        std::slice::from_ref(m),
        &[(0, cutoff)],
        &PvalueOptions::default(),
    )
    .ctx("pvalue_to_score")?[0];
    if !thr.is_finite() {
        return Ok(Vec::new());
    }
    let alph = &m.alphabet;
    let custom = matches!(alph, Alphabet::Custom(_));
    let rna = *alph == Alphabet::Rna;
    let strands = if alph.is_nucleotide() && rc { 2 } else { 1 };
    let places = if use_gaps {
        placements(m)
    } else {
        vec![Placement {
            offsets: (0..m.width()).collect(),
            span: m.width(),
        }]
    };
    let mut cands: Vec<HitKey> = Vec::new();
    let mut words: Vec<String> = Vec::new();
    for (si, s) in seqs.iter().enumerate() {
        let chars: Vec<char> = s.chars().collect();
        for pl in &places {
            if pl.span > chars.len() {
                continue;
            }
            for p in 0..=chars.len() - pl.span {
                let win = &chars[p..p + pl.span];
                if !win.iter().all(|&c| letter_index(alph, c).is_some()) {
                    continue;
                }
                for strand in 0..strands {
                    let oriented: Vec<char> = if strand == 1 {
                        rc_chars(win, rna)
                    } else if custom {
                        win.to_vec()
                    } else {
                        win.iter().map(|c| c.to_ascii_uppercase()).collect()
                    };
                    let mut matched = vec!['.'; pl.span];
                    for &o in &pl.offsets {
                        matched[o] = oriented[o];
                    }
                    words.push(pl.offsets.iter().map(|&o| oriented[o]).collect());
                    cands.push(HitKey {
                        motif: mi,
                        seq: si,
                        start: p,
                        end: p + pl.span,
                        minus: strand == 1,
                        matched: matched.into_iter().collect(),
                        score: 0.0,
                    });
                }
            }
        }
    }
    if words.is_empty() {
        return Ok(Vec::new());
    }
    let scores = score_match(m, &words, false).ctx("score_match")?;
    Ok(cands
        .into_iter()
        .zip(scores)
        .filter(|(_, s)| *s >= thr)
        .map(|(mut h, s)| {
            h.score = s;
            h
        })
        .collect())
}

fn first_difference(a: &[HitKey], b: &[HitKey]) -> String {
    for i in 0..a.len().max(b.len()) {
        if a.get(i) != b.get(i) {
            return format!(
                "first difference at {i}: scan {:?}, brute force {:?}",
                a.get(i),
                b.get(i)
            );
        }
    }
    String::new()
}

// ---------------------------------------------------------------------------
// Type conversions

#[test]
fn bug_uniform_column_negative_ic() {
    let alph = Alphabet::custom("ABCDE".chars()).unwrap();
    let pcm = Motif::builder(Matrix::from_columns(&[[8.0; 5]]).unwrap(), alph)
        .motif_type(Some(MotifType::Pcm))
        .pseudocount(0.01)
        .build()
        .unwrap();
    // Fails with "ICM values cannot be negative": log2(5) less the
    // column's entropy comes out just below zero.
    convert_type(&pcm, MotifType::Icm, None, false).unwrap();
}

#[test]
fn prop_pcm_to_ppm_matches_formula() {
    run("pcm_to_ppm", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=20);
        let m = gen_motif(rng, &alph, w, MotifType::Pcm, PLAIN);
        let ppm = convert_type(&m, MotifType::Ppm, None, false).ctx("to PPM")?;
        let n = alph.len() as f64;
        let pc = m.pseudocount;
        for (j, (c, p)) in m.matrix.columns().zip(ppm.matrix.columns()).enumerate() {
            let total: f64 = c.iter().sum();
            let s: f64 = p.iter().sum();
            ensure!(
                close(s, 1.0, 1e-12, 0.0),
                "column {j} sums to {s}\n{}",
                desc(&m)
            );
            for (&x, &y) in c.iter().zip(p) {
                let want = (x + pc / n) / (total + pc);
                ensure!(
                    close(y, want, 1e-14, 1e-12),
                    "column {j}: got {y}, want {want}\n{}",
                    desc(&m)
                );
            }
        }
        Ok(())
    });
}

/// Whether an ICM has a column of zeros, which no conversion accepts; see
/// `bug_icm_zero_column_round_trip`.
fn zero_icm_column(m: &Motif) -> bool {
    m.motif_type == MotifType::Icm && m.matrix.columns().any(|c| c.iter().all(|&v| v == 0.0))
}

#[test]
fn bug_icm_zero_column_round_trip() {
    let cols = [[0.25; 4], [0.7, 0.1, 0.1, 0.1]];
    let ppm = Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .build()
        .unwrap();
    let icm = convert_type(&ppm, MotifType::Icm, None, false).unwrap();
    assert_eq!(icm.matrix.col(0), [0.0; 4]);
    // Fails with "PPM columns must sum to 1": 0 / 0 becomes 0 rather than
    // the uniform column that a column without information content is.
    let back = convert_type(&icm, MotifType::Ppm, None, false).unwrap();
    assert_eq!(back.matrix.col(0), [0.25; 4]);
}

#[test]
fn prop_type_round_trips() {
    run("type_round_trips", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=20);
        // Without a pseudocount every conversion can be undone.
        let mut p = gen_motif(rng, &alph, w, MotifType::Ppm, PLAIN).into_parts();
        p.pseudocount = 0.0;
        let ppm = Motif::from_parts(p).unwrap();
        let pwm = convert_type(&ppm, MotifType::Pwm, None, false).ctx("to PWM")?;
        let back = convert_type(&pwm, MotifType::Ppm, None, false).ctx("PWM to PPM")?;
        let icm = convert_type(&pwm, MotifType::Icm, None, false).ctx("PWM to ICM")?;
        let n = alph.len() as f64;
        for (j, ((a, b), h)) in ppm
            .matrix
            .columns()
            .zip(back.matrix.columns())
            .zip(icm.matrix.columns())
            .enumerate()
        {
            for (x, y) in a.iter().zip(b) {
                ensure!(
                    close(*x, *y, 1e-12, 1e-9),
                    "PPM -> PWM -> PPM column {j}: {a:?} became {b:?}\n{}",
                    desc(&ppm)
                );
            }
            // Column heights add up to the column's information content.
            let ent: f64 = a.iter().filter(|&&x| x > 0.0).map(|&x| -x * x.log2()).sum();
            let ic = n.log2() - ent;
            let height: f64 = h.iter().sum();
            ensure!(
                close(height, ic, 1e-9, 1e-9),
                "ICM column {j} is {height} high, IC {ic}\n{}",
                desc(&ppm)
            );
        }
        // ICM back to PPM, except where a column holds no information.
        if !zero_icm_column(&icm) {
            let back = convert_type(&icm, MotifType::Ppm, None, false).ctx("ICM to PPM")?;
            for (a, b) in ppm.matrix.as_slice().iter().zip(back.matrix.as_slice()) {
                ensure!(
                    close(*a, *b, 1e-9, 1e-9),
                    "PPM -> ICM -> PPM: {a} became {b}\n{}",
                    desc(&ppm)
                );
            }
        }
        // Counts through every type and back to probabilities.
        let mut p = gen_motif(rng, &alph, w, MotifType::Pcm, PLAIN).into_parts();
        p.pseudocount = 0.0;
        let pcm = Motif::from_parts(p).unwrap();
        let mut m = convert_type(&pcm, MotifType::Ppm, None, false).ctx("PCM to PPM")?;
        let direct = m.clone();
        let icm = convert_type(&direct, MotifType::Icm, None, false).ctx("to ICM")?;
        if !zero_icm_column(&icm) {
            for t in [MotifType::Pwm, MotifType::Icm, MotifType::Ppm] {
                m = convert_type(&m, t, None, false).ctx("chain")?;
            }
            for (a, b) in direct.matrix.as_slice().iter().zip(m.matrix.as_slice()) {
                ensure!(
                    close(*a, *b, 1e-9, 1e-9),
                    "PCM -> PPM -> PWM -> ICM -> PPM: {a} became {b}\n{}",
                    desc(&pcm)
                );
            }
        }
        Ok(())
    });
}

#[test]
fn prop_conversions_give_valid_motifs() {
    run("conversions_valid", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=20);
        let ty = gen_type(rng, true);
        let m = gen_motif(
            rng,
            &alph,
            w,
            ty,
            Gen {
                raw: true,
                gaps: 0.3,
            },
        );
        ensure!(
            convert_type(&m, ty, None, false).ctx("same type")? == m,
            "converting to the same type changed the motif\n{}",
            desc(&m)
        );
        ensure!(
            ty == MotifType::Cwm || convert_type(&m, MotifType::Cwm, None, false).is_err(),
            "converted to a CWM\n{}",
            desc(&m)
        );
        for to in [
            MotifType::Pcm,
            MotifType::Ppm,
            MotifType::Pwm,
            MotifType::Icm,
        ] {
            let c = match convert_type(&m, to, None, false) {
                // See bug_uniform_column_negative_ic.
                Err(_) if to == MotifType::Icm && uniform_columns(&m) > 0 => continue,
                r => r.ctx(&format!("{ty} to {to}\n{}", desc(&m)))?,
            };
            ensure!(
                c.width() == m.width() && c.gap == m.gap && c.bkg == m.bkg,
                "{ty} to {to} changed the shape\n{}",
                desc(&m)
            );
            match to {
                MotifType::Ppm => {
                    // Weights whose powers of two sum to within 0.01 of 1 are
                    // taken as probabilities as they are.
                    for s in c.matrix.col_sums() {
                        ensure!(
                            close(s, 1.0, 0.01, 0.0),
                            "PPM column sums to {s}\n{}",
                            desc(&m)
                        );
                    }
                }
                MotifType::Pcm => {
                    let want = m.nsites.filter(|&n| n > 1.0).unwrap_or(100.0);
                    let want = if ty == MotifType::Pcm {
                        c.matrix.col_sums()[0]
                    } else {
                        want
                    };
                    for s in c.matrix.col_sums() {
                        ensure!(
                            s == want,
                            "PCM column sums to {s}, not {want}\n{}",
                            desc(&m)
                        );
                    }
                    ensure!(
                        c.matrix
                            .as_slice()
                            .iter()
                            .all(|&v| v >= 0.0 && v == v.round()),
                        "PCM has a value that is not a count: {:?}\n{}",
                        c.matrix.as_slice(),
                        desc(&m)
                    );
                }
                _ => {}
            }
        }
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// Motif operations

const DNA_IUPAC_COMPLEMENT: [(char, char); 15] = [
    ('A', 'T'),
    ('C', 'G'),
    ('G', 'C'),
    ('T', 'A'),
    ('U', 'A'),
    ('M', 'K'),
    ('K', 'M'),
    ('R', 'Y'),
    ('Y', 'R'),
    ('W', 'W'),
    ('S', 'S'),
    ('H', 'D'),
    ('D', 'H'),
    ('B', 'V'),
    ('V', 'B'),
];

fn iupac_complement(c: char, rna: bool) -> char {
    if c == 'A' && rna {
        return 'U';
    }
    DNA_IUPAC_COMPLEMENT
        .iter()
        .find(|(a, _)| *a == c)
        .map_or(c, |(_, b)| *b)
}

fn nucleotide(rng: &mut StdRng) -> Alphabet {
    if rng.random_bool(0.7) {
        Alphabet::Dna
    } else {
        Alphabet::Rna
    }
}

#[test]
fn prop_reverse_complement() {
    run("reverse_complement", 300, |rng| {
        let alph = nucleotide(rng);
        let w = rng.random_range(1..=20);
        let ty = gen_type(rng, true);
        let m = gen_motif(
            rng,
            &alph,
            w,
            ty,
            Gen {
                raw: true,
                gaps: 0.3,
            },
        );
        // Conversions carry the consensus over from the original type, so
        // start from one computed from the matrix, as reverse_complement()
        // computes it.
        let mut p = m.into_parts();
        p.recompute_derived();
        let m = Motif::from_parts(p).unwrap();
        let rc = reverse_complement(&m, false).ctx("rc")?;
        let back = reverse_complement(&rc, false).ctx("rc twice")?;
        let changed: Vec<&str> = [
            ("matrix", back.matrix == m.matrix),
            ("gap", back.gap == m.gap),
            ("bkg", back.bkg == m.bkg),
            ("consensus", back.consensus == m.consensus),
            ("nsites", back.nsites == m.nsites),
            ("pseudocount", back.pseudocount == m.pseudocount),
            ("type", back.motif_type == m.motif_type),
            ("name", back.name == m.name),
        ]
        .iter()
        .filter(|f| !f.1)
        .map(|f| f.0)
        .collect();
        ensure!(
            changed.is_empty(),
            "reverse complementing twice changed {changed:?}: consensus {:?} became {:?}\n{}",
            m.consensus,
            back.consensus,
            desc(&m)
        );
        ensure!(
            close(back.icscore, m.icscore, 1e-9, 1e-9) || ty == MotifType::Cwm,
            "total IC went from {} to {}\n{}",
            m.icscore,
            back.icscore,
            desc(&m)
        );
        // Column j of the reverse complement is column w-1-j with its rows
        // reversed, since A, C, G, T complement to T, G, C, A.
        for j in 0..w {
            let mut c = m.matrix.col(w - 1 - j).to_vec();
            c.reverse();
            ensure!(c == rc.matrix.col(j), "column {j} is wrong\n{}", desc(&m));
        }
        let rna = alph == Alphabet::Rna;
        let want: String = m
            .consensus
            .as_deref()
            .unwrap()
            .chars()
            .rev()
            .map(|c| iupac_complement(c, rna))
            .collect();
        ensure!(
            rc.consensus.as_deref() == Some(want.as_str()),
            "consensus of the reverse complement is {:?}, not {want:?} (from {:?})\n{}",
            rc.consensus,
            m.consensus,
            desc(&m)
        );
        // Gaps are mirrored: a gap after k columns comes after w - k.
        let mut locs: Vec<usize> = m.gap.gaploc.iter().map(|g| w - g).collect();
        locs.reverse();
        ensure!(rc.gap.gaploc == locs, "gaps not mirrored\n{}", desc(&m));
        // Other alphabets are refused unless asked, and then only reversed.
        let custom = custom_alphabet(rng, 2, 6, POOL);
        let other = pick(rng, &[Alphabet::Aa, custom]);
        let o = gen_motif(rng, &other, w, MotifType::Ppm, PLAIN);
        ensure!(
            reverse_complement(&o, false).is_err(),
            "reverse complemented {other}"
        );
        let r = reverse_complement(&o, true).ctx("ignore alphabet")?;
        let mut v = o.matrix.as_slice().to_vec();
        v.reverse();
        ensure!(
            r.matrix.as_slice() == v,
            "{other} motif not simply reversed"
        );
        Ok(())
    });
}

/// A DNA motif equal to its own reverse complement.
fn palindrome(rng: &mut StdRng, w: usize) -> Motif {
    // No zeros and no pseudocount: a pseudocount goes through whole counts,
    // whose rounding gives ties to the earlier letter and so breaks the
    // symmetry between a column and its complement.
    let half: Vec<Vec<f64>> = (0..w.div_ceil(2))
        .map(|_| {
            let c: Vec<f64> = gen_ppm_col(rng, 4).iter().map(|v| v + 0.01).collect();
            let s: f64 = c.iter().sum();
            c.iter().map(|v| v / s).collect()
        })
        .collect();
    let mut cols: Vec<Vec<f64>> = vec![Vec::new(); w];
    for (j, c) in half.iter().enumerate() {
        let mut r = c.clone();
        r.reverse();
        if j == w - 1 - j {
            // The middle column of an odd width must be its own complement.
            let a = (c[0] + c[3]) / 2.0;
            let b = (c[1] + c[2]) / 2.0;
            cols[j] = vec![a, b, b, a];
        } else {
            cols[j] = c.clone();
            cols[w - 1 - j] = r;
        }
    }
    Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .pseudocount(0.0)
        .nsites(gen_nsites(rng))
        .build()
        .unwrap()
}

#[test]
fn prop_palindromes() {
    run("palindromes", 150, |rng| {
        let w = rng.random_range(1..=14);
        let m = palindrome(rng, w);
        let rc = reverse_complement(&m, false).ctx("rc")?;
        ensure!(rc.matrix == m.matrix, "not a palindrome\n{}", desc(&m));
        // The two strands of a palindrome give the same hits.
        let seqs = gen_seqs(rng, &Alphabet::Dna, 4, 80);
        let opts = ScanOptions {
            pvalue: 10f64.powf(rng.random_range(-3.5..-0.5)),
            ..Default::default()
        };
        let t = scan_table(std::slice::from_ref(&m), &seqs, &opts).ctx("scan")?;
        let mut plus: Vec<(u32, i64, i64, u64)> = Vec::new();
        let mut minus = Vec::new();
        for i in 0..t.len() {
            let k = (t.sequence[i], t.start[i], t.end[i], t.score[i].to_bits());
            if t.minus[i] {
                minus.push(k);
            } else {
                plus.push(k);
            }
        }
        plus.sort_unstable();
        minus.sort_unstable();
        ensure!(
            plus == minus,
            "strands differ: {plus:?} vs {minus:?}\n{}\nseqs {seqs:?}",
            desc(&m)
        );
        Ok(())
    });
}

#[test]
fn prop_trim_motif() {
    run("trim_motif", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=20);
        let ty = gen_type(rng, false);
        let m = gen_motif(
            rng,
            &alph,
            w,
            ty,
            Gen {
                raw: true,
                gaps: 0.2,
            },
        );
        let min_ic = rng.random_range(0.0..1.5);
        let from = pick(rng, &[TrimFrom::Both, TrimFrom::Left, TrimFrom::Right]);
        let nsites = m.nsites.unwrap_or(100.0);
        let ic: Vec<f64> = m
            .matrix
            .columns()
            .map(|c| position_icscore(c, &m.bkg, ty, m.pseudocount, nsites, false))
            .collect();
        let got = trim_motif(&m, min_ic, from).ctx("trim")?;
        let left = if from != TrimFrom::Right {
            ic.iter().position(|&v| v >= min_ic).unwrap_or(w)
        } else {
            0
        };
        let right = if from != TrimFrom::Left {
            ic.iter().rposition(|&v| v >= min_ic).map_or(0, |i| i + 1)
        } else {
            w
        };
        match got {
            None => ensure!(
                right <= left,
                "trimmed to nothing, but {left}..{right} should remain\n{}",
                desc(&m)
            ),
            Some(t) => {
                ensure!(
                    left < right && t.matrix == m.matrix.slice_columns(left..right),
                    "kept {} columns, want {left}..{right} (IC {ic:?})\n{}",
                    t.width(),
                    desc(&m)
                );
                for (k, &g) in t.gap.gaploc.iter().enumerate() {
                    ensure!(
                        g > 0 && g < t.width() && m.gap.gaploc.contains(&(g + left)),
                        "gap {k} at {g} does not match the original\n{}",
                        desc(&m)
                    );
                }
            }
        }
        Ok(())
    });
}

#[test]
fn prop_round_and_normalize() {
    run("round_normalize", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=20);
        let m = gen_motif(rng, &alph, w, MotifType::Ppm, PLAIN);
        let n = alph.len() as f64;
        let tol = rng.random_range(0.0..(1.0 / n).min(0.06));
        let r = round_motif(&m, tol).ctx("round")?;
        ensure!(r.motif_type == MotifType::Ppm, "round changed the type");
        for (j, c) in r.matrix.columns().enumerate() {
            let s: f64 = c.iter().sum();
            ensure!(close(s, 1.0, 1e-9, 0.0), "rounded column {j} sums to {s}");
            ensure!(
                c.iter().all(|&v| v == 0.0 || v >= tol),
                "rounded column {j} keeps a value below {tol}: {c:?}\n{}",
                desc(&m)
            );
        }
        if tol == 0.0 {
            for (a, b) in r.matrix.as_slice().iter().zip(m.matrix.as_slice()) {
                ensure!(
                    close(*a, *b, 1e-12, 1e-12),
                    "rounding at 0 changed {b} to {a}"
                );
            }
        }
        // Normalizing leaves no zero probabilities.
        let z = normalize(&m).ctx("normalize")?;
        ensure!(
            z.matrix.as_slice().iter().all(|&v| v > 0.0),
            "normalize left a zero: {:?}\n{}",
            z.matrix.as_slice(),
            desc(&m)
        );
        for s in z.matrix.col_sums() {
            ensure!(close(s, 1.0, 1e-9, 0.0), "normalized column sums to {s}");
        }
        // A PCM keeps its counts and takes a pseudocount of 1 if it had none.
        let pcm = gen_motif(rng, &alph, w, MotifType::Pcm, PLAIN);
        let n = normalize(&pcm).ctx("normalize PCM")?;
        let want = if pcm.pseudocount == 0.0 {
            1.0
        } else {
            pcm.pseudocount
        };
        ensure!(
            n.matrix == pcm.matrix && n.pseudocount == want,
            "normalize changed a PCM's counts or gave it pseudocount {}\n{}",
            n.pseudocount,
            desc(&pcm)
        );
        Ok(())
    });
}

#[test]
fn prop_gaps_and_alphabet_switch() {
    run("gaps_switch", 300, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(2..=20);
        let ty = gen_type(rng, true);
        let m = gen_motif(rng, &alph, w, ty, PLAIN);
        let g = gen_gaps(rng, w);
        let gapped = add_gap(&m, Some(&g.gaploc), &g.mingap, &g.maxgap).ctx("add_gap")?;
        ensure!(gapped.gap == g, "add_gap gave {:?}, not {g:?}", gapped.gap);
        ensure!(gapped.matrix == m.matrix, "add_gap changed the matrix");
        let off = ungap(&gapped, false);
        ensure!(
            !off.gap.isgapped && off.gap.gaploc == g.gaploc,
            "ungap without deleting lost the gaps: {:?}",
            off.gap
        );
        let restored = add_gap(
            &off,
            Some(&off.gap.gaploc),
            &off.gap.mingap,
            &off.gap.maxgap,
        )
        .ctx("again")?;
        ensure!(restored == gapped, "gaps not restored");
        ensure!(
            ungap(&gapped, true) == m,
            "deleting gaps did not give the motif back"
        );

        if alph.is_nucleotide() {
            let s = switch_alph(&m).ctx("switch")?;
            ensure!(
                s.alphabet != m.alphabet && s.matrix == m.matrix,
                "switch_alph"
            );
            let (from, to) = if alph == Alphabet::Dna {
                ('T', 'U')
            } else {
                ('U', 'T')
            };
            let cons = s.consensus.as_deref().unwrap_or("");
            ensure!(
                !cons.contains(from)
                    && cons
                        == m.consensus
                            .as_deref()
                            .unwrap()
                            .replace(from, &to.to_string()),
                "consensus {:?} after switching from {:?}",
                s.consensus,
                m.consensus
            );
            ensure!(switch_alph(&s).ctx("switch back")? == m, "switching twice");
        } else {
            ensure!(switch_alph(&m).is_err(), "switched a {} motif", alph);
        }
        Ok(())
    });
}

/// Whether every value is 0 or at least 1, as in a PPM whose columns each
/// hold a single letter, which the motif builder takes for counts.
fn one_hot(m: &Motif) -> bool {
    m.matrix.as_slice().iter().all(|&v| v == 0.0 || v >= 1.0)
}

fn one_hot_ppm() -> Motif {
    let cols = [[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]];
    Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .pseudocount(1.0)
        .nsites(Some(10.0))
        .build()
        .unwrap()
}

#[test]
fn bug_merge_one_hot_gets_pseudocount() {
    let m = one_hot_ppm();
    let merged = merge_motifs(&[m.clone(), m.clone()], &MergeOptions::default()).unwrap();
    // The merged matrix is all zeros and ones, so from_matrix() takes it for
    // counts and converts it to probabilities with the pseudocount: the
    // first column comes out as [0.932, 0.023, 0.023, 0.023]. The R package
    // returns the motif unchanged.
    assert_eq!(merged.matrix, m.matrix);
}

#[test]
fn bug_shuffle_motifs_changes_type() {
    let m = one_hot_ppm();
    let out = shuffle_motifs(std::slice::from_ref(&m), 1, 1).unwrap();
    // The motif builder guesses the type of the shuffled matrix and takes
    // zeros and ones for counts, which are rescaled to nsites.
    assert_eq!(out[0].motif_type, MotifType::Ppm);
}

#[test]
fn prop_merge_with_copies() {
    run("merge_copies", 150, |rng| {
        let alph = gen_alphabet(rng);
        let w = rng.random_range(1..=14);
        let ty = gen_type(rng, false);
        let m = gen_motif(rng, &alph, w, ty, PLAIN);
        let copies = vec![m.clone(); rng.random_range(2..=4)];
        let opts = MergeOptions {
            min_overlap: rng.random_range(1..=6),
            rc: rng.random_bool(0.5),
            weighted: rng.random_bool(0.5),
            ..Default::default()
        };
        let ppm = convert_type(&m, MotifType::Ppm, None, false).ctx("to PPM")?;
        if one_hot(&ppm) {
            // See bug_merge_one_hot_gets_pseudocount.
            return Ok(());
        }
        let merged = merge_motifs(&copies, &opts).ctx("merge")?;
        ensure!(
            merged.width() == w,
            "merged width {}, not {w}\n{}",
            merged.width(),
            desc(&m)
        );
        for (j, (a, b)) in merged
            .matrix
            .columns()
            .zip(ppm.matrix.columns())
            .enumerate()
        {
            for (x, y) in a.iter().zip(b) {
                ensure!(
                    close(*x, *y, 1e-9, 1e-9),
                    "merged column {j} is {a:?}, the motif's {b:?}\n{}",
                    desc(&m)
                );
            }
        }
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// P-values and scores

/// A motif small enough to list every word: DNA or RNA up to 5 wide, amino
/// acids up to 3 and custom alphabets up to about 4096 words.
fn small_motif(rng: &mut StdRng) -> Motif {
    let alph = gen_alphabet(rng);
    let max_w = match alph.len() {
        4 => 5,
        20 => 3,
        n => ((4096f64).ln() / (n as f64).ln()).floor() as usize,
    };
    let w = rng.random_range(1..=max_w.max(1));
    let ty = gen_type(rng, true);
    let m = gen_motif(rng, &alph, w, ty, PLAIN);
    finite(rng, m)
}

/// The motif with a pseudocount if its PWM would otherwise have infinite
/// entries, where `score_match()` and scanning disagree (see
/// `bug_score_match_differs_from_scan`).
fn finite(rng: &mut StdRng, m: Motif) -> Motif {
    let pwm = convert_type(&m, MotifType::Pwm, None, false).unwrap();
    if pwm.matrix.as_slice().iter().all(|v| v.is_finite()) {
        return m;
    }
    let mut p = m.into_parts();
    p.pseudocount = pick(rng, &[0.01, 0.5, 1.0]);
    p.recompute_derived();
    Motif::from_parts(p).unwrap()
}

#[test]
fn bug_score_match_differs_from_scan() {
    let bkg = vec![
        0.1789181122284846,
        0.21665700955142825,
        0.30256653247854715,
        0.30185834574154,
    ];
    let cols = [[0.25; 4], [1.0, 0.0, 0.0, 0.0]];
    let m = Motif::builder(Matrix::from_columns(&cols).unwrap(), Alphabet::Dna)
        .motif_type(Some(MotifType::Ppm))
        .pseudocount(0.0)
        .nsites(Some(2.0))
        .bkg(Some(bkg))
        .build()
        .unwrap();
    let words = ["AA", "CA", "GA", "TA"];
    let opts = ScanOptions {
        pvalue: 0.999,
        rc: false,
        ..Default::default()
    };
    let hits = scan_table(std::slice::from_ref(&m), &words, &opts).unwrap();
    // score_match() gives 0.965, 0.689, 2.528 and 2.532, the scan 3.286,
    // 3.01, 0.207 and 0.211. score_match() adds the pseudocount to the PWM
    // with normalize(), scanning and motif_pvalue() to the PPM, and the
    // uniform column's counts (half a site per letter) are rounded up for
    // different letters depending on the rounding error of each path.
    assert_eq!(score_match(&m, &words, false).unwrap(), hits.score);
}

#[test]
fn prop_pvalues_match_enumeration() {
    run("pvalue_enumeration", 200, |rng| {
        let m = small_motif(rng);
        let words = all_words(m.alphabet.letters(), m.width());
        let scores = score_match(&m, &words, false).ctx("score_match")?;
        let probs: Vec<f64> = words
            .iter()
            .map(|w| {
                w.chars()
                    .map(|c| m.bkg[letter_index(&m.alphabet, c).unwrap()])
                    .product()
            })
            .collect();
        let total: f64 = probs.iter().sum();
        let mut by_score: BTreeMap<i64, f64> = BTreeMap::new();
        for (s, p) in scores.iter().zip(&probs) {
            *by_score.entry((s * 1000.0).round() as i64).or_default() += p / total;
        }
        // P(score >= s) for each score, highest first.
        let mut tail = Vec::new();
        let mut acc = 0.0;
        for (&s, &p) in by_score.iter().rev() {
            acc += p;
            tail.push((s as f64 / 1000.0, acc));
        }
        let one = std::slice::from_ref(&m);
        let opts = PvalueOptions::default();
        let queries: Vec<(usize, f64)> = tail.iter().map(|&(s, _)| (0, s)).collect();
        let got = motif_pvalue(one, &queries, &opts).ctx(&format!("motif_pvalue\n{}", desc(&m)))?;
        for ((s, want), g) in tail.iter().zip(&got) {
            ensure!(
                close(*g, *want, 1e-12, 1e-9),
                "P-value of {s} is {g}, enumeration gives {want}\n{}",
                desc(&m)
            );
        }
        // Scores between grid points and beyond the range.
        let (lo, hi) = (tail.last().unwrap().0, tail[0].0);
        let range = motif_range(&m, false).ctx("motif_range")?;
        ensure!(
            range == (lo, hi),
            "motif_range {range:?}, enumeration ({lo}, {hi})\n{}",
            desc(&m)
        );
        let edge =
            motif_pvalue(one, &[(0, hi + 0.001), (0, lo - 1.0), (0, lo)], &opts).ctx("edges")?;
        ensure!(
            edge[0] == 0.0 && close(edge[1], 1.0, 1e-12, 0.0) && close(edge[2], 1.0, 1e-12, 0.0),
            "P-values beyond the range: {edge:?}\n{}",
            desc(&m)
        );
        let minp = min_pvalues(one, &opts).ctx("min_pvalues")?[0];
        ensure!(
            close(minp, tail[0].1, 1e-15, 1e-9),
            "min_pvalues {minp}, enumeration {}",
            tail[0].1
        );
        // Thresholds: the lowest score whose P-value is within the cutoff.
        let ps: Vec<f64> = (0..8)
            .map(|_| 10f64.powf(rng.random_range(-6.0..0.0)))
            .chain([1.0, tail[0].1, tail[tail.len() / 2].1])
            .collect();
        let queries: Vec<(usize, f64)> = ps.iter().map(|&p| (0, p)).collect();
        let thr = pvalue_to_score(one, &queries, &opts).ctx("pvalue_to_score")?;
        for (&p, &t) in ps.iter().zip(&thr) {
            // tail is ordered by decreasing score and increasing P-value. The
            // threshold is one grid step above the highest score whose
            // P-value is above the cutoff, which need not be a score any
            // word reaches; it is the lowest score of all if there is none,
            // and infinite if even the best score is too likely.
            let expect = |p: f64| {
                if tail[0].1 > p {
                    return f64::INFINITY;
                }
                match tail.iter().find(|&&(_, q)| q > p) {
                    Some(&(s, _)) => ((s * 1000.0).round() + 1.0) / 1000.0,
                    None => lo,
                }
            };
            let (want, alt) = (expect(p * (1.0 + 1e-9)), expect(p * (1.0 - 1e-9)));
            ensure!(
                t == want || t == alt,
                "threshold for {p} is {t}, enumeration gives {want}\n{}",
                desc(&m)
            );
            if t.is_finite() {
                let back = motif_pvalue(one, &[(0, t)], &opts).ctx("back")?[0];
                ensure!(
                    back <= p * (1.0 + 1e-9),
                    "P-value {back} of threshold {t} above {p}"
                );
            }
        }
        Ok(())
    });
}

#[test]
fn prop_pvalues_monotone() {
    run("pvalue_monotone", 150, |rng| {
        let alph = gen_alphabet(rng);
        let w = if alph == Alphabet::Aa {
            rng.random_range(1..=8)
        } else {
            rng.random_range(1..=20)
        };
        let ty = gen_type(rng, true);
        let m = gen_motif(rng, &alph, w, ty, PLAIN);
        let (lo, hi) = motif_range(&m, false).ctx("range")?;
        let one = std::slice::from_ref(&m);
        let opts = PvalueOptions::default();
        // Thresholds for some cutoffs, then the P-values of random scores, of
        // the best score and of each threshold and the grid point below it,
        // all in one call.
        let ps: Vec<f64> = (0..10)
            .map(|_| 10f64.powf(rng.random_range(-8.0..0.0)))
            .collect();
        let qs: Vec<(usize, f64)> = ps.iter().map(|&x| (0, x)).collect();
        let t = pvalue_to_score(one, &qs, &opts).ctx("thresholds")?;
        let mut scores: Vec<f64> = (0..40)
            .map(|_| rng.random_range(lo - 1.0..hi + 1.0))
            .collect();
        scores.sort_by(f64::total_cmp);
        let mut q: Vec<(usize, f64)> = scores.iter().map(|&s| (0, s)).collect();
        q.push((0, hi));
        for &s in t.iter().filter(|s| s.is_finite()) {
            q.extend([(0, s), (0, s - 0.001)]);
        }
        let p = motif_pvalue(one, &q, &opts).ctx(&format!("motif_pvalue\n{}", desc(&m)))?;
        for k in 1..scores.len() {
            ensure!(
                p[k] <= p[k - 1] && (0.0..=1.0).contains(&p[k]),
                "P-values not decreasing: {} at {} after {} at {}\n{}",
                p[k],
                scores[k],
                p[k - 1],
                scores[k - 1],
                desc(&m)
            );
        }
        // Thresholds and P-values are consistent on the grid.
        let best = p[scores.len()];
        let mut at = p[scores.len() + 1..].chunks(2);
        for (&pv, &s) in ps.iter().zip(&t) {
            if !s.is_finite() {
                ensure!(
                    best > pv,
                    "no threshold for {pv}, but the best score gives {best}"
                );
                continue;
            }
            let at = at.next().unwrap();
            ensure!(
                at[0] <= pv && (s <= lo || at[1] > pv),
                "threshold {s} for {pv}: P-values {at:?} at it and one step below\n{}",
                desc(&m)
            );
        }
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// Scanning

fn gen_scan_motifs(rng: &mut StdRng, alph: &Alphabet, gaps: f64) -> Vec<Motif> {
    let k = rng.random_range(1..=3);
    let max_w = if *alph == Alphabet::Aa { 6 } else { 16 };
    (0..k)
        .map(|_| {
            let w = rng.random_range(1..=max_w);
            let ty = gen_type(rng, true);
            let m = gen_motif(rng, alph, w, ty, Gen { raw: true, gaps });
            finite(rng, m)
        })
        .collect()
}

#[test]
fn prop_scan_matches_brute_force() {
    run("scan_brute_force", 120, |rng| {
        let alph = gen_alphabet(rng);
        let motifs = gen_scan_motifs(rng, &alph, 0.3);
        let count = rng.random_range(1..=6);
        let seqs = gen_seqs(rng, &alph, count, 70);
        let opts = ScanOptions {
            pvalue: 10f64.powf(rng.random_range(-4.0..-0.3)),
            rc: rng.random_bool(0.7),
            use_gaps: rng.random_bool(0.8),
            ..Default::default()
        };
        let ctx = || format!("{}\nseqs {seqs:?}\noptions {opts:?}", descs(&motifs));
        let t = scan_table(&motifs, &seqs, &opts).ctx("scan")?;
        let got = table_keys(&t, &seqs);
        let mut want = Vec::new();
        for (i, m) in motifs.iter().enumerate() {
            want.extend(brute_hits(
                i,
                m,
                &seqs,
                opts.pvalue,
                opts.rc,
                opts.use_gaps,
            )?);
        }
        sort_keys(&mut want);
        ensure!(
            got == want,
            "{} hits from the scan, {} by brute force; {}\n{}",
            got.len(),
            want.len(),
            first_difference(&got, &want),
            ctx()
        );
        // Each hit's P-value is that of its score.
        for (mi, m) in motifs.iter().enumerate() {
            let rows: Vec<usize> = (0..t.len()).filter(|&i| t.motif[i] == mi as u32).collect();
            let q: Vec<(usize, f64)> = rows.iter().map(|&i| (0, t.score[i])).collect();
            let p = motif_pvalue(std::slice::from_ref(m), &q, &PvalueOptions::default())
                .ctx("pvalues")?;
            for (&i, &pv) in rows.iter().zip(&p) {
                ensure!(
                    close(t.pvalue[i], pv, 0.0, 1e-12) && t.pvalue[i] <= opts.pvalue,
                    "hit {i} has P-value {}, its score {} gives {pv}\n{}",
                    t.pvalue[i],
                    t.score[i],
                    ctx()
                );
            }
        }
        // Hits are ordered by motif, sequence and start.
        for i in 1..t.len() {
            let a = (t.motif[i - 1], t.sequence[i - 1], t.start[i - 1]);
            let b = (t.motif[i], t.sequence[i], t.start[i]);
            ensure!(a <= b, "hits out of order at {i}");
        }
        // Thread counts and encoding first make no difference; checked for
        // half the cases to save time.
        if rng.random_bool(0.5) {
            return Ok(());
        }
        let threaded = scan_table(
            &motifs,
            &seqs,
            &ScanOptions {
                nthreads: 3,
                ..opts.clone()
            },
        )
        .ctx("threaded scan")?;
        ensure!(threaded == t, "3 threads gave different hits\n{}", ctx());
        let enc = EncodedSequences::new(&seqs, &alph, 1).ctx("encode")?;
        let encoded = scan_encoded(&motifs, &enc, &opts).ctx("encoded scan")?;
        ensure!(
            encoded == t,
            "scan_encoded differs from scan_table\n{}",
            ctx()
        );
        for i in 0..t.len() {
            ensure!(
                t.matched_encoded(i, &enc) == t.matched(i, &seqs)
                    || matches!(alph, Alphabet::Custom(_)),
                "matched letters of hit {i} differ when decoded"
            );
        }
        Ok(())
    });
}

#[test]
fn prop_scan_infinite_pwms() {
    // Motifs with zero probabilities and no pseudocount, whose PWMs have
    // infinite entries until a pseudocount is added. The scan and
    // motif_pvalue() must agree on the scores and P-values.
    run("scan_infinite_pwm", 100, |rng| {
        let alph = if rng.random_bool(0.6) {
            nucleotide(rng)
        } else {
            gen_alphabet(rng)
        };
        let w = rng.random_range(1..=8);
        let ty = pick(rng, &[MotifType::Ppm, MotifType::Pcm, MotifType::Icm]);
        let mut p = gen_motif(rng, &alph, w, ty, PLAIN).into_parts();
        p.pseudocount = 0.0;
        let m = Motif::from_parts(p).unwrap();
        let pwm = convert_type(&m, MotifType::Pwm, None, false).ctx("PWM")?;
        if pwm.matrix.as_slice().iter().all(|v| v.is_finite()) {
            return Ok(());
        }
        let seqs = gen_seqs(rng, &alph, 3, 60);
        let opts = ScanOptions {
            pvalue: 10f64.powf(rng.random_range(-3.0..-0.3)),
            ..Default::default()
        };
        let one = std::slice::from_ref(&m);
        let t = scan_table(one, &seqs, &opts).ctx("scan")?;
        let q: Vec<(usize, f64)> = t.score.iter().map(|&s| (0, s)).collect();
        let pv = motif_pvalue(one, &q, &PvalueOptions::default()).ctx("pvalue")?;
        for i in 0..t.len() {
            ensure!(
                close(t.pvalue[i], pv[i], 0.0, 1e-12),
                "hit {i} has P-value {}, motif_pvalue() gives {} for its score {}\n{}\nseqs {seqs:?}",
                t.pvalue[i],
                pv[i],
                t.score[i],
                desc(&m)
            );
        }
        // Scores are grid sums of the PWM motif_pvalue() uses.
        let sm = universalmotif::pvalue::scoring_matrix(&m, None).ctx("scoring_matrix")?;
        for (i, h) in table_keys(&t, &seqs).iter().enumerate() {
            let grid: i64 = h
                .matched
                .chars()
                .enumerate()
                .map(|(j, c)| {
                    let r = letter_index(&alph, c).unwrap();
                    (sm.pwm.get(r, j) * 1000.0) as i64
                })
                .sum();
            ensure!(
                grid as f64 / 1000.0 == h.score,
                "hit {i} scores {}, the PWM of motif_pvalue() gives {}\n{}\nseqs {seqs:?}",
                h.score,
                grid as f64 / 1000.0,
                desc(&m)
            );
        }
        Ok(())
    });
}

/// Reverse complement a sequence, keeping letters outside the alphabet as
/// they are (so `U` stays `U` in DNA).
fn rc_string(s: &str, rna: bool) -> String {
    let alph = if rna { Alphabet::Rna } else { Alphabet::Dna };
    s.chars()
        .rev()
        .map(|c| {
            if letter_index(&alph, c).is_some() {
                complement(c, rna)
            } else {
                c
            }
        })
        .collect()
}

#[test]
fn prop_scan_reverse_complement_mirrors() {
    run("scan_rc_mirror", 100, |rng| {
        let alph = nucleotide(rng);
        let rna = alph == Alphabet::Rna;
        let motifs = gen_scan_motifs(rng, &alph, 0.3);
        let count = rng.random_range(1..=5);
        let seqs = gen_seqs(rng, &alph, count, 80);
        let rcs: Vec<String> = seqs.iter().map(|s| rc_string(s, rna)).collect();
        let opts = ScanOptions {
            pvalue: 10f64.powf(rng.random_range(-4.0..-0.5)),
            ..Default::default()
        };
        let a = scan_table(&motifs, &seqs, &opts).ctx("scan")?;
        let b = scan_table(&motifs, &rcs, &opts).ctx("scan rc")?;
        let flip = |t: &HitTable, s: &[String]| {
            let mut v: Vec<HitKey> = table_keys(t, s);
            for h in v.iter_mut() {
                let len = s[h.seq].chars().count();
                (h.start, h.end) = (len - h.end, len - h.start);
                h.minus = !h.minus;
            }
            sort_keys(&mut v);
            v
        };
        let want = flip(&b, &rcs);
        let got = table_keys(&a, &seqs);
        ensure!(
            got == want,
            "hits on the reverse complement do not mirror: {}\n{}\nseqs {seqs:?}",
            first_difference(&got, &want),
            descs(&motifs)
        );
        Ok(())
    });
}

#[test]
fn prop_scan_overlap_removal() {
    run("scan_no_overlaps", 100, |rng| {
        let alph = gen_alphabet(rng);
        let motifs = gen_scan_motifs(rng, &alph, 0.2);
        let count = rng.random_range(1..=4);
        let seqs = gen_seqs(rng, &alph, count, 90);
        let no = NoOverlaps {
            by: pick(rng, &[OverlapPriority::Pvalue, OverlapPriority::Score]),
            ignore_strand: rng.random_bool(0.5),
            ignore_motif: rng.random_bool(0.5),
        };
        let opts = ScanOptions {
            pvalue: 10f64.powf(rng.random_range(-3.0..-0.3)),
            ..Default::default()
        };
        let full = scan_table(&motifs, &seqs, &opts).ctx("scan")?;
        let kept = scan_table(
            &motifs,
            &seqs,
            &ScanOptions {
                no_overlaps: Some(no),
                ..opts.clone()
            },
        )
        .ctx("scan without overlaps")?;
        let all = table_keys(&full, &seqs);
        let keep = table_keys(&kept, &seqs);
        let ctx = || format!("{no:?}\n{}\nseqs {seqs:?}", descs(&motifs));
        ensure!(
            keep.iter().all(|k| all.contains(k)),
            "kept a hit the full scan does not have\n{}",
            ctx()
        );
        let pv: HashMap<(usize, usize, usize, bool, String), f64> = (0..full.len())
            .map(|i| {
                (
                    (
                        full.motif[i] as usize,
                        full.sequence[i] as usize,
                        full.start[i] as usize,
                        full.minus[i],
                        full.matched(i, &seqs),
                    ),
                    full.pvalue[i],
                )
            })
            .collect();
        let group = |h: &HitKey| {
            (
                h.seq,
                if no.ignore_motif { 0 } else { h.motif },
                !no.ignore_strand && h.minus,
            )
        };
        let prio = |h: &HitKey| match no.by {
            OverlapPriority::Pvalue => -pv[&(h.motif, h.seq, h.start, h.minus, h.matched.clone())],
            OverlapPriority::Score => h.score,
        };
        let overlap = |a: &HitKey, b: &HitKey| a.start < b.end && b.start < a.end;
        for (i, a) in keep.iter().enumerate() {
            for b in &keep[i + 1..] {
                ensure!(
                    group(a) != group(b) || !overlap(a, b),
                    "kept overlapping hits {a:?} and {b:?}\n{}",
                    ctx()
                );
            }
        }
        for h in all.iter().filter(|h| !keep.contains(h)) {
            ensure!(
                keep.iter()
                    .any(|k| group(k) == group(h) && overlap(k, h) && prio(k) >= prio(h)),
                "dropped {h:?} without a better overlapping hit\n{}",
                ctx()
            );
        }
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// Comparison and clustering

fn uniform_columns(m: &Motif) -> usize {
    let ppm = convert_type(m, MotifType::Ppm, None, false).unwrap();
    let n = m.alphabet.len() as f64;
    ppm.matrix
        .columns()
        .filter(|c| c.iter().map(|v| (v - 1.0 / n).powi(2)).sum::<f64>() < 1e-10)
        .count()
}

#[test]
fn prop_comparisons() {
    run("compare", 100, |rng| {
        let alph = gen_alphabet(rng);
        let k = rng.random_range(2..=5);
        let motifs: Vec<Motif> = (0..k)
            .map(|_| {
                let w = rng.random_range(1..=12);
                let ty = gen_type(rng, false);
                gen_motif(rng, &alph, w, ty, PLAIN)
            })
            .collect();
        let opts = CompareOptions {
            min_overlap: rng.random_range(1..=6),
            rc: rng.random_bool(0.6),
            ..Default::default()
        };
        let ctx = || format!("{opts:?}\n{}", descs(&motifs));
        let r = compare_matrix(&motifs, &opts, true).ctx("compare")?;
        for i in 0..k {
            for j in 0..k {
                let (s, p) = (r.scores[i][j], r.pvalues[i][j]);
                ensure!(
                    close(s, r.scores[j][i], 1e-9, 0.0),
                    "scores not symmetric: {s} vs {} at {i},{j}\n{}",
                    r.scores[j][i],
                    ctx()
                );
                ensure!(
                    (-1.0 - 1e-9..=1.0 + 1e-9).contains(&s) && (0.0..=1.0).contains(&p),
                    "score {s} or P-value {p} out of range\n{}",
                    ctx()
                );
            }
            if uniform_columns(&motifs[i]) == 0 {
                ensure!(
                    close(r.scores[i][i], 1.0, 1e-9, 0.0),
                    "self-score {}\n{}",
                    r.scores[i][i],
                    ctx()
                );
            }
        }
        // Permuting the motifs permutes the results.
        let mut perm: Vec<usize> = (0..k).collect();
        for i in (1..k).rev() {
            perm.swap(i, rng.random_range(0..=i));
        }
        let shuffled: Vec<Motif> = perm.iter().map(|&i| motifs[i].clone()).collect();
        let r2 = compare_matrix(&shuffled, &opts, true).ctx("compare permuted")?;
        for a in 0..k {
            for b in 0..k {
                let (i, j) = (perm[a], perm[b]);
                ensure!(
                    close(r2.scores[a][b], r.scores[i][j], 1e-12, 0.0)
                        && close(r2.pvalues[a][b], r.pvalues[i][j], 1e-15, 1e-6),
                    "permuting changed {i},{j}: score {} vs {}, P-value {} vs {}\n{}",
                    r2.scores[a][b],
                    r.scores[i][j],
                    r2.pvalues[a][b],
                    r.pvalues[i][j],
                    ctx()
                );
            }
        }
        let threaded = compare_matrix(
            &motifs,
            &CompareOptions {
                nthreads: 3,
                ..opts.clone()
            },
            true,
        )
        .ctx("threaded")?;
        ensure!(
            format!("{threaded:?}") == format!("{r:?}"),
            "3 threads gave different results\n{}",
            ctx()
        );
        // A motif and its reverse complement score as the motif with itself.
        // The motif is taken as probabilities first, since the reverse
        // complement keeps the background as it is, which changes the
        // probabilities a PWM stands for.
        if alph.is_nucleotide() && opts.rc {
            let m = &convert_type(&motifs[0], MotifType::Ppm, None, false).ctx("PPM")?;
            let pair = [m.clone(), reverse_complement(m, false).ctx("rc")?];
            let s = compare_matrix(&pair, &opts, false)
                .ctx("compare rc")?
                .scores;
            ensure!(
                close(s[0][1], s[0][0], 1e-9, 0.0),
                "against its reverse complement {}, against itself {}\n{}",
                s[0][1],
                s[0][0],
                desc(m)
            );
        }
        Ok(())
    });
}

/// Parse a Newick tree, returning its leaf labels and checking that every
/// node is a pair with branch lengths.
fn newick_leaves(s: &str) -> Result<Vec<String>, String> {
    fn node(b: &[u8], i: &mut usize, leaves: &mut Vec<String>) -> Result<(), String> {
        if b.get(*i) == Some(&b'(') {
            *i += 1;
            node(b, i, leaves)?;
            if b.get(*i) != Some(&b',') {
                return Err(format!("expected ',' at {i}"));
            }
            *i += 1;
            node(b, i, leaves)?;
            if b.get(*i) != Some(&b')') {
                return Err(format!("expected ')' at {i}"));
            }
            *i += 1;
        } else {
            let start = *i;
            while *i < b.len() && !b":,();".contains(&b[*i]) {
                *i += 1;
            }
            leaves.push(String::from_utf8_lossy(&b[start..*i]).into_owned());
        }
        Ok(())
    }
    fn length(b: &[u8], i: &mut usize) -> Result<f64, String> {
        if b.get(*i) != Some(&b':') {
            return Err(format!("expected ':' at {i}"));
        }
        *i += 1;
        let start = *i;
        while *i < b.len() && !b",();".contains(&b[*i]) {
            *i += 1;
        }
        std::str::from_utf8(&b[start..*i])
            .unwrap()
            .parse()
            .map_err(|e| format!("branch length: {e}"))
    }
    // Lengths follow every child; parse them by stripping them in a second
    // pass over a copy with lengths checked separately.
    let b = s.as_bytes();
    let mut stripped = Vec::with_capacity(b.len());
    let mut i = 0;
    while i < b.len() {
        if b[i] == b':' {
            let v = length(b, &mut i)?;
            if !v.is_finite() {
                return Err(format!("branch length {v}"));
            }
        } else {
            stripped.push(b[i]);
            i += 1;
        }
    }
    if stripped.last() != Some(&b';') {
        return Err("no closing ';'".into());
    }
    stripped.pop();
    let mut leaves = Vec::new();
    let mut i = 0;
    node(&stripped, &mut i, &mut leaves)?;
    if i != stripped.len() {
        return Err(format!("trailing text at {i}"));
    }
    Ok(leaves)
}

/// Heights of single linkage merges: the edges of a minimum spanning tree.
fn mst_weights(d: &[Vec<f64>]) -> Vec<f64> {
    let n = d.len();
    let mut in_tree = vec![false; n];
    let mut best = vec![f64::INFINITY; n];
    best[0] = 0.0;
    let mut out = Vec::new();
    for step in 0..n {
        let u = (0..n)
            .filter(|&i| !in_tree[i])
            .min_by(|&a, &b| best[a].total_cmp(&best[b]))
            .unwrap();
        in_tree[u] = true;
        if step > 0 {
            out.push(best[u]);
        }
        for v in 0..n {
            if !in_tree[v] {
                best[v] = best[v].min(d[u.max(v)][u.min(v)]);
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

#[test]
fn prop_hclust() {
    run("hclust", 300, |rng| {
        let n = rng.random_range(2..=14);
        let ties = rng.random_bool(0.3);
        let mut d = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..i {
                let v: f64 = rng.random_range(0.0..1.0);
                let v = if ties { (v * 5.0).round() / 5.0 } else { v };
                d[i][j] = v;
                d[j][i] = v;
            }
        }
        let methods = [
            Linkage::WardD,
            Linkage::Single,
            Linkage::Complete,
            Linkage::Average,
            Linkage::Mcquitty,
            Linkage::Median,
            Linkage::Centroid,
            Linkage::WardD2,
        ];
        for method in methods {
            let ctx = || format!("{method} on {d:?}");
            let t = hclust(&d, method).ctx("hclust")?;
            ensure!(
                t.merge.len() == n - 1,
                "{} merges\n{}",
                t.merge.len(),
                ctx()
            );
            let mut order = t.order.clone();
            order.sort_unstable();
            ensure!(order == (0..n).collect::<Vec<_>>(), "order {:?}", t.order);
            let monotone = !matches!(method, Linkage::Median | Linkage::Centroid);
            if monotone {
                ensure!(
                    t.height.windows(2).all(|h| h[0] <= h[1] + 1e-12),
                    "heights not increasing: {:?}\n{}",
                    t.height,
                    ctx()
                );
            }
            let mut prev: Option<Vec<usize>> = None;
            for k in 1..=n {
                let g = t.cutree(k).ctx("cutree")?;
                let mut next = 0;
                for &x in &g {
                    ensure!(x <= next, "groups not numbered by first appearance: {g:?}");
                    if x == next {
                        next += 1;
                    }
                }
                ensure!(
                    next == k,
                    "cutree({k}) gave {next} groups: {g:?}\n{}",
                    ctx()
                );
                // Cutting into more groups splits groups, never joins them.
                if let Some(p) = &prev {
                    for a in 0..n {
                        for b in 0..n {
                            ensure!(
                                g[a] != g[b] || p[a] == p[b],
                                "cutree({k}) joins {a} and {b}, split at {}\n{}",
                                k - 1,
                                ctx()
                            );
                        }
                    }
                }
                prev = Some(g);
            }
            if monotone {
                let h = pick(rng, &t.height) + if rng.random_bool(0.5) { 0.0 } else { 1e-6 };
                let g = t.cutree_height(h).ctx("cutree_height")?;
                let groups = g.iter().max().unwrap() + 1;
                let want = n - t.height.iter().filter(|&&x| x <= h).count();
                ensure!(groups == want, "cut at {h}: {groups} groups, want {want}");
            }
            if method == Linkage::Single {
                let mut hs = t.height.clone();
                hs.sort_by(f64::total_cmp);
                let mst = mst_weights(&d);
                ensure!(
                    hs.iter().zip(&mst).all(|(a, b)| close(*a, *b, 1e-12, 0.0)),
                    "single linkage heights {hs:?}, spanning tree {mst:?}"
                );
            }
            let labels: Vec<String> = (0..n).map(|i| format!("leaf{i}")).collect();
            let mut leaves = newick_leaves(&t.newick(&labels))
                .map_err(|e| format!("{e}: {}\n{}", t.newick(&labels), ctx()))?;
            leaves.sort();
            let mut want = labels.clone();
            want.sort();
            ensure!(leaves == want, "Newick leaves {leaves:?}");
            let link = t.linkage();
            ensure!(
                link.last().unwrap()[3] == n as f64,
                "last merge holds {} leaves",
                link.last().unwrap()[3]
            );
        }
        Ok(())
    });
}

#[test]
fn prop_motif_tree() {
    run("motif_tree", 40, |rng| {
        let alph = nucleotide(rng);
        let n = rng.random_range(2..=8);
        let motifs: Vec<Motif> = (0..n)
            .map(|_| {
                let w = rng.random_range(3..=12);
                gen_motif(rng, &alph, w, MotifType::Ppm, PLAIN)
            })
            .collect();
        let t = motif_tree(
            &motifs,
            &TreeOptions {
                min_overlap: rng.random_range(1..=6),
                ..Default::default()
            },
        )
        .ctx("motif_tree")?;
        let labels: Vec<String> = motifs.iter().map(|m| m.name.clone()).collect();
        let leaves = newick_leaves(&t.newick(&labels))?;
        ensure!(leaves.len() == n, "{} leaves for {n} motifs", leaves.len());
        ensure!(
            t.height.iter().all(|&h| (0.0..=1.0 + 1e-9).contains(&h)),
            "heights {:?}",
            t.height
        );
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// Sequences

fn sorted_chars(s: &str) -> Vec<char> {
    let mut v: Vec<char> = s.chars().collect();
    v.sort_unstable();
    v
}

#[test]
fn prop_shuffles() {
    run("shuffles", 200, |rng| {
        let alph = gen_alphabet(rng);
        let k = rng.random_range(1..=4);
        let seqs: Vec<String> = (0..rng.random_range(1..=4))
            .map(|_| {
                let len = rng.random_range(k..=120);
                gen_seq(rng, &alph, len, 0.05)
            })
            .collect();
        let seed = rng.random::<u64>();
        let ctx = || format!("k {k} seed {seed} seqs {seqs:?}");
        let euler = shuffle_sequences(
            &seqs,
            &ShuffleOptions {
                k,
                seed,
                ..Default::default()
            },
        )
        .ctx("euler")?;
        for (a, b) in seqs.iter().zip(&euler) {
            let letters = sequence_letters(&[a]);
            ensure!(
                count_kmers(a, k, &letters) == count_kmers(b, k, &letters),
                "Euler shuffle changed the {k}-mers: {a} -> {b}\n{}",
                ctx()
            );
            let (ac, bc): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
            ensure!(
                ac[..k - 1] == bc[..k - 1] && ac[ac.len() + 1 - k..] == bc[bc.len() + 1 - k..],
                "Euler shuffle moved the ends: {a} -> {b}"
            );
        }
        for method in [ShuffleMethod::Linear, ShuffleMethod::Markov] {
            let out = shuffle_sequences(
                &seqs,
                &ShuffleOptions {
                    k,
                    seed,
                    method,
                    ..Default::default()
                },
            )
            .ctx("shuffle")?;
            for (a, b) in seqs.iter().zip(&out) {
                let ok = if method == ShuffleMethod::Linear {
                    sorted_chars(a) == sorted_chars(b)
                } else {
                    let l = sequence_letters(&[a]);
                    b.chars().count() == a.chars().count() && b.chars().all(|c| l.contains(&c))
                };
                ensure!(ok, "{method:?} shuffle: {a} -> {b}\n{}", ctx());
            }
        }
        // Windows shuffle each part, keeping the letters overall.
        let size = rng.random_range(k.max(2)..=30) as f64;
        let overlap = rng.random_range(0..size as usize) as f64;
        let windowed = shuffle_sequences(
            &seqs,
            &ShuffleOptions {
                k,
                seed,
                window: Some(universalmotif::kmer::WindowOptions { size, overlap }),
                nthreads: 2,
                ..Default::default()
            },
        );
        if let Ok(w) = windowed {
            for (a, b) in seqs.iter().zip(&w) {
                ensure!(
                    sorted_chars(a) == sorted_chars(b),
                    "windowed shuffle changed the letters: {a} -> {b}"
                );
            }
        }
        // Results do not depend on the thread count.
        let threaded = shuffle_sequences(
            &seqs,
            &ShuffleOptions {
                k,
                seed,
                nthreads: 3,
                ..Default::default()
            },
        )
        .ctx("threaded")?;
        ensure!(
            threaded == euler,
            "3 threads shuffled differently\n{}",
            ctx()
        );
        Ok(())
    });
}

#[test]
fn prop_shuffle_motifs_keeps_columns() {
    run("shuffle_motifs", 150, |rng| {
        let alph = gen_alphabet(rng);
        let motifs: Vec<Motif> = (0..rng.random_range(1..=4))
            .map(|_| {
                let w = rng.random_range(1..=10);
                // Columns holding one letter can make a shuffled motif look
                // like counts; see bug_shuffle_motifs_changes_type.
                loop {
                    let m = gen_motif(rng, &alph, w, MotifType::Ppm, PLAIN);
                    if m.matrix
                        .columns()
                        .all(|c| c.iter().any(|&v| v > 0.0 && v < 1.0))
                    {
                        break m;
                    }
                }
            })
            .collect();
        let k = rng.random_range(1..=3);
        let out = shuffle_motifs(&motifs, k, rng.random()).ctx("shuffle_motifs")?;
        let cols = |ms: &[Motif]| {
            let mut v: Vec<Vec<u64>> = ms
                .iter()
                .flat_map(|m| {
                    m.matrix
                        .columns()
                        .map(|c| c.iter().map(|x| x.to_bits()).collect::<Vec<_>>())
                        .collect::<Vec<_>>()
                })
                .collect();
            v.sort();
            v
        };
        ensure!(
            cols(&out) == cols(&motifs),
            "shuffling motifs changed the columns\n{}",
            descs(&motifs)
        );
        for (a, b) in motifs.iter().zip(&out) {
            ensure!(
                a.width() == b.width() && a.motif_type == b.motif_type,
                "widths or types changed"
            );
        }
        Ok(())
    });
}

#[test]
fn prop_kmer_counts() {
    run("kmer_counts", 300, |rng| {
        let alph = gen_alphabet(rng);
        let len = rng.random_range(0..150);
        let seq = gen_seq(rng, &alph, len, 0.08);
        let k = rng.random_range(1..=3);
        // Count over a subset of the letters, so some k-mers are left out.
        let mut letters: Vec<char> = alph
            .letters()
            .iter()
            .copied()
            .filter(|_| rng.random_bool(0.8))
            .collect();
        if letters.is_empty() {
            letters.push(alph.letters()[0]);
        }
        let words = kmers(&letters, k);
        let index: HashMap<&str, usize> = words
            .iter()
            .enumerate()
            .map(|(i, w)| (w.as_str(), i))
            .collect();
        let mut want = vec![0u64; words.len()];
        let chars: Vec<char> = seq.chars().collect();
        for w in chars.windows(k) {
            let s: String = w.iter().collect();
            if let Some(&i) = index.get(s.as_str()) {
                want[i] += 1;
            }
        }
        let got = count_kmers(&seq, k, &letters);
        ensure!(
            got == want,
            "counts {got:?}, brute force {want:?} for {seq:?}"
        );
        Ok(())
    });
}

#[test]
fn prop_get_bkg_counts_windows() {
    run("get_bkg", 200, |rng| {
        let dna = rng.random_bool(0.5);
        let alph = if dna {
            Alphabet::Dna
        } else {
            gen_alphabet(rng)
        };
        let kmax = rng.random_range(1..=3);
        let seqs: Vec<String> = (0..rng.random_range(1..=4))
            .map(|_| {
                let len = rng.random_range(kmax..=80);
                gen_seq(rng, &alph, len, 0.03)
            })
            .collect();
        let rc = dna && rng.random_bool(0.5);
        let merge = rng.random_bool(0.6);
        let explicit = rng.random_bool(0.3);
        let opts = BkgOptions {
            k: (1..=kmax).collect(),
            rc,
            merge,
            alphabet: explicit.then(|| alph.letters().to_vec()),
            ..Default::default()
        };
        let rows = get_bkg(&seqs, &[], &opts).ctx("get_bkg")?;
        let mut all: Vec<String> = seqs.clone();
        if rc {
            // Every letter is complemented here, U to A included.
            all.extend(seqs.iter().map(|s| {
                s.chars()
                    .rev()
                    .map(|c| complement(c, false))
                    .collect::<String>()
            }));
        }
        let letters: Vec<char> = match &opts.alphabet {
            Some(a) => a.clone(),
            None => sequence_letters(&all),
        };
        let windows = |s: &str, k: usize| {
            let c: Vec<char> = s.chars().collect();
            c.windows(k)
                .filter(|w| w.iter().all(|x| letters.contains(x)))
                .count() as u64
        };
        for k in 1..=kmax {
            let of_k: Vec<_> = rows
                .iter()
                .filter(|r| r.kmer.chars().count() == k)
                .collect();
            if merge {
                let total: u64 = of_k.iter().map(|r| r.count).sum();
                let want: u64 = all.iter().map(|s| windows(s, k)).sum();
                ensure!(
                    total == want,
                    "k {k}: counts sum to {total}, {want} windows\nseqs {seqs:?} {opts:?}"
                );
                let p: f64 = of_k.iter().map(|r| r.probability).sum();
                ensure!(
                    want == 0 || close(p, 1.0, 1e-9, 0.0),
                    "probabilities sum to {p}"
                );
            } else {
                for (i, s) in all.iter().enumerate().take(seqs.len()) {
                    let total: u64 = of_k
                        .iter()
                        .filter(|r| r.sequence == Some(i))
                        .map(|r| r.count)
                        .sum();
                    ensure!(
                        total == windows(s, k),
                        "k {k} sequence {i}: {total} counted\nseqs {seqs:?} {opts:?}"
                    );
                }
            }
        }
        Ok(())
    });
}

#[test]
fn prop_encoding_round_trips() {
    run("encoding", 300, |rng| {
        let alph = gen_alphabet(rng);
        let count = rng.random_range(1..=5);
        let seqs = gen_seqs(rng, &alph, count, 200);
        let enc = EncodedSequences::new(&seqs, &alph, rng.random_range(1..=2)).ctx("encode")?;
        ensure!(enc.len() == seqs.len(), "{} sequences encoded", enc.len());
        let custom = matches!(alph, Alphabet::Custom(_));
        for (i, s) in seqs.iter().enumerate() {
            let want: String = s
                .chars()
                .map(|c| match letter_index(&alph, c) {
                    Some(r) if !custom => alph.letters()[r],
                    _ => c,
                })
                .collect();
            ensure!(
                enc.decode(i) == want,
                "decoded {:?}, want {want:?}",
                enc.decode(i)
            );
            ensure!(enc.seq_len(i) == s.len(), "length {}", enc.seq_len(i));
            // Pieces of any size encode the same as the whole.
            let mut piecewise = EncodedSequences::empty(&alph).ctx("empty")?;
            let mut r = piecewise.record();
            let b = s.as_bytes();
            let mut at = 0;
            while at < b.len() {
                let step = rng.random_range(1..=9).min(b.len() - at);
                r.push(&b[at..at + step]);
                at += step;
            }
            piecewise.push(r);
            ensure!(piecewise.decode(0) == want, "piecewise encoding of {s:?}");
            if s.len() > 1 {
                let (a, z) = (rng.random_range(0..s.len()), s.len());
                let mut out = Vec::new();
                enc.decode_into(i, a, z, &mut out);
                ensure!(out == want.as_bytes()[a..z], "decoding {a}..{z} of {s:?}");
            }
        }
        ensure!(
            enc.total_len() == seqs.iter().map(String::len).sum::<usize>(),
            "total length"
        );
        Ok(())
    });
}

#[test]
fn prop_masking() {
    run("masking", 300, |rng| {
        let alph = if rng.random_bool(0.7) {
            nucleotide(rng)
        } else {
            Alphabet::Aa
        };
        let len = rng.random_range(0..120);
        let seq = gen_seq(rng, &alph, len, 0.05);
        let n = seq.len() as i64;
        let ranges: Vec<(i64, i64)> = (0..rng.random_range(0..5))
            .map(|_| {
                let a = rng.random_range(0..=n);
                let b = rng.random_range(0..=n);
                (a.min(b), a.max(b))
            })
            .collect();
        let masked = mask_ranges(&seq, &ranges, '-').ctx("mask_ranges")?;
        let want: String = seq
            .chars()
            .enumerate()
            .map(|(i, c)| {
                if ranges.iter().any(|&(a, b)| (a..b).contains(&(i as i64))) {
                    '-'
                } else {
                    c
                }
            })
            .collect();
        ensure!(
            masked == want,
            "mask_ranges {ranges:?} of {seq:?} gave {masked:?}"
        );
        if n > 0 {
            ensure!(
                mask_ranges(&seq, &[(0, n + 1)], '-').is_err(),
                "masking past the end"
            );
        }
        // Pattern masking, then the reverse complement on what is left.
        let plen = rng.random_range(1..=3);
        let pattern = gen_seq(rng, &alph, plen, 0.0).to_ascii_uppercase();
        let rc = alph.is_nucleotide() && rng.random_bool(0.5);
        let got = mask_seqs(&[&seq], &pattern, rc, 'N', &alph).ctx("mask_seqs")?;
        let mask = |s: &[u8], p: &[u8]| -> Vec<u8> {
            let mut out = s.to_vec();
            if s.len() >= p.len() {
                for i in 0..=s.len() - p.len() {
                    if s[i..i + p.len()].eq_ignore_ascii_case(p) {
                        out[i..i + p.len()].fill(b'N');
                    }
                }
            }
            out
        };
        let mut want = mask(seq.as_bytes(), pattern.as_bytes());
        if rc {
            let rp = rc_string(&pattern, alph == Alphabet::Rna);
            want = mask(&want, rp.as_bytes());
        }
        ensure!(
            got[0].as_bytes() == want,
            "masking {pattern} (rc {rc}) in {seq:?} gave {:?}",
            got[0]
        );
        Ok(())
    });
}

#[test]
fn prop_complexity_bounds() {
    run("complexity", 300, |rng| {
        let alph = gen_alphabet(rng);
        let len = rng.random_range(2..200);
        let seq = gen_seq(rng, &alph, len, 0.05);
        let own: String = sequence_letters(&[&seq]).iter().collect();
        let distinct = own.chars().count();
        let max_word = rng.random_range(1..=8);
        for method in [
            ComplexityMethod::WoottonFederhen,
            ComplexityMethod::WoottonFederhenFast,
            ComplexityMethod::Trifonov,
            ComplexityMethod::TrifonovFast,
        ] {
            // Against its own letters, as the default alphabet of DNA
            // would not cover other letters in the sequence.
            let c = complexity(&seq, method, Some(&own), max_word);
            if distinct < 2
                && matches!(
                    method,
                    ComplexityMethod::WoottonFederhen | ComplexityMethod::WoottonFederhenFast
                )
            {
                continue;
            }
            ensure!(
                (-1e-12..=1.0 + 1e-12).contains(&c),
                "{method} complexity {c} of {seq:?}"
            );
        }
        if seq.len() >= 4 {
            let d = complexity(&seq, ComplexityMethod::Dust, None, 3);
            ensure!(d >= 0.0, "DUST {d} of {seq:?}");
        }
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// File formats

fn names_and_metadata(rng: &mut StdRng, m: Motif, i: usize) -> Motif {
    let mut p = m.into_parts();
    p.name = format!("motif_{i}_{}", rng.random_range(0..1000));
    p.altname = rng.random_bool(0.5).then(|| format!("alt{i}"));
    p.family = rng.random_bool(0.5).then(|| format!("fam{i}"));
    p.organism = rng.random_bool(0.5).then(|| format!("org{i}"));
    Motif::from_parts(p).unwrap()
}

fn gen_io_motifs(rng: &mut StdRng, alph: &Alphabet) -> Vec<Motif> {
    (0..rng.random_range(1..=4))
        .map(|i| {
            let w = rng.random_range(1..=20);
            let ty = gen_type(rng, false);
            let m = gen_motif(rng, alph, w, ty, PLAIN);
            names_and_metadata(rng, m, i)
        })
        .collect()
}

fn same_matrix(a: &Motif, b: &Motif, tol: f64) -> bool {
    a.matrix.nrow() == b.matrix.nrow()
        && a.matrix.ncol() == b.matrix.ncol()
        && a.matrix
            .as_slice()
            .iter()
            .zip(b.matrix.as_slice())
            .all(|(x, y)| (x - y).abs() <= tol)
}

fn as_type(ms: &[Motif], t: MotifType) -> Vec<Motif> {
    ms.iter()
        .map(|m| convert_type(m, t, None, false).unwrap())
        .collect()
}

#[test]
fn bug_transfac_single_site_counts() {
    let text = "ID one\nP0\tA\tC\tG\tT\n01\t1\t0\t0\t0\n02\t0\t1\t0\t0\nXX\n//\n";
    let m = &transfac::read(text).unwrap()[0];
    // The counts become probabilities and then counts again, and an nsites
    // of 1 or less is taken as 100 on the way back, so the motif has
    // columns of 100 with nsites 1. The R package does the same.
    assert_eq!(m.nsites, Some(1.0));
    assert_eq!(m.matrix.col_sums(), [1.0, 1.0]);
}

#[test]
fn prop_formats_round_trip() {
    run("formats", 150, |rng| {
        let alph = gen_alphabet(rng);
        let motifs = gen_io_motifs(rng, &alph);
        let ctx = || descs(&motifs);
        let nuc = alph.is_nucleotide();
        let custom = matches!(alph, Alphabet::Custom(_));
        let ppms = as_type(&motifs, MotifType::Ppm);
        let pcms = as_type(&motifs, MotifType::Pcm);
        let check = |fmt: &str,
                     back: Vec<Motif>,
                     want: &[Motif],
                     tol: f64,
                     alt: bool,
                     meta: bool|
         -> Check {
            ensure!(
                back.len() == want.len(),
                "{fmt}: read {} motifs of {}\n{}",
                back.len(),
                want.len(),
                ctx()
            );
            for (b, w) in back.iter().zip(want) {
                ensure!(
                    b.name == w.name,
                    "{fmt}: name {:?} became {:?}",
                    w.name,
                    b.name
                );
                ensure!(b.alphabet == w.alphabet, "{fmt}: alphabet {}", b.alphabet);
                ensure!(
                    same_matrix(b, w, tol),
                    "{fmt}: matrix of {} changed\nread {}\n{}",
                    w.name,
                    desc(b),
                    ctx()
                );
                if alt {
                    ensure!(b.altname == w.altname, "{fmt}: altname {:?}", b.altname);
                }
                if meta {
                    ensure!(
                        b.family == w.family && b.organism == w.organism,
                        "{fmt}: family or organism"
                    );
                }
            }
            Ok(())
        };

        // MEME: probabilities to six decimals.
        if !custom || alph.letters().iter().all(|c| c.is_ascii_alphanumeric()) {
            let text = meme::write(&motifs, &meme::WriteOptions::default()).ctx("meme write")?;
            let back = meme::read(&text, &meme::ReadOptions::default()).ctx("meme read")?;
            check("MEME", back.clone(), &ppms, 1e-6, true, false)?;
            for (b, m) in back.iter().zip(&motifs) {
                ensure!(
                    b.nsites == Some(m.nsites.unwrap_or(100.0)),
                    "MEME: nsites {:?} became {:?}",
                    m.nsites,
                    b.nsites
                );
            }
        }
        // JASPAR: counts.
        if !custom {
            let text = jaspar::write(&motifs).ctx("jaspar write")?;
            check(
                "JASPAR",
                jaspar::read(&text).ctx("jaspar read")?,
                &pcms,
                0.0,
                true,
                false,
            )?;
        } else {
            ensure!(
                jaspar::write(&motifs).is_err(),
                "JASPAR wrote a custom alphabet"
            );
        }
        if nuc {
            let text = transfac::write(&motifs, &transfac::WriteOptions::default())
                .ctx("transfac write")?;
            let back = transfac::read(&text).ctx("transfac read")?;
            // TRANSFAC files hold DNA.
            let dna: Vec<Motif> = pcms
                .iter()
                .map(|m| {
                    if m.alphabet == Alphabet::Rna {
                        switch_alph(m).unwrap()
                    } else {
                        m.clone()
                    }
                })
                .collect();
            // Counts of a single site come back as 100; see
            // bug_transfac_single_site_counts.
            if pcms.iter().all(|m| m.matrix.col_sums()[0] > 1.0) {
                check("TRANSFAC", back, &dna, 0.0, true, true)?;
            }
            let text = cisbp::write(&motifs).ctx("cisbp write")?;
            check(
                "CIS-BP",
                cisbp::read(&text).ctx("cisbp read")?,
                &ppms,
                0.0,
                true,
                true,
            )?;
        }
        if alph == Alphabet::Dna {
            let text = homer::write(&motifs, &homer::Threshold::default()).ctx("homer write")?;
            check(
                "HOMER",
                homer::read(&text).ctx("homer read")?,
                &ppms,
                1.5e-3,
                false,
                false,
            )?;
            let text = uniprobe::write(&motifs).ctx("uniprobe write")?;
            check(
                "UniPROBE",
                uniprobe::read(&text).ctx("uniprobe read")?,
                &ppms,
                0.0,
                false,
                false,
            )?;
        }
        // Plain matrices, as rows or columns.
        let positions = pick(
            rng,
            &[matrix_io::Positions::Columns, matrix_io::Positions::Rows],
        );
        let text = matrix_io::write(
            &motifs,
            &matrix_io::WriteOptions {
                positions,
                motif_type: Some(MotifType::Ppm),
                digits: 8,
                ..Default::default()
            },
        )
        .ctx("matrix write")?;
        let back = matrix_io::read(
            &text,
            &matrix_io::ReadOptions {
                positions,
                alphabet: alph.clone(),
                motif_type: Some(MotifType::Ppm),
                ..Default::default()
            },
        )
        .ctx("matrix read")?;
        check("matrix", back, &ppms, 1e-7, false, false)?;
        Ok(())
    });
}

// ---------------------------------------------------------------------------
// Statistics

fn lchoose(n: u64, k: u64) -> f64 {
    libm::lgamma(n as f64 + 1.0) - libm::lgamma(k as f64 + 1.0) - libm::lgamma((n - k) as f64 + 1.0)
}

#[test]
fn prop_binomial() {
    run("binomial", 300, |rng| {
        let n = rng.random_range(0..=300u64);
        let r = rng.random_range(0.0..1.0);
        let p = pick(rng, &[r, 1e-4, 0.5, 0.999]);
        let d: Vec<f64> = (0..=n).map(|x| stats::dbinom(x, n, p)).collect();
        let total: f64 = d.iter().sum();
        ensure!(
            close(total, 1.0, 1e-10, 0.0),
            "dbinom sums to {total} (n {n}, p {p})"
        );
        let mut prev = 1.0;
        for k in 0..=n + 1 {
            let up = stats::pbinom_upper(k, n, p);
            let want: f64 = d.iter().skip(k as usize).sum();
            ensure!(
                (0.0..=1.0).contains(&up) && up <= prev + 1e-15,
                "pbinom_upper({k}, {n}, {p}) = {up} after {prev}"
            );
            ensure!(
                close(up, want, 1e-12, 1e-8),
                "pbinom_upper({k}, {n}, {p}) = {up}, sum {want}"
            );
            prev = up;
        }
        Ok(())
    });
}

#[test]
fn prop_hypergeometric_and_fisher() {
    run("hypergeometric", 300, |rng| {
        let nr = rng.random_range(0..=120u64);
        let nb = rng.random_range(0..=120u64);
        let n = rng.random_range(0..=nr + nb);
        let lo = n.saturating_sub(nb);
        let hi = n.min(nr);
        let dens: Vec<f64> = (lo..=hi).map(|x| stats::dhyper(x, nr, nb, n)).collect();
        let total: f64 = dens.iter().sum();
        ensure!(
            close(total, 1.0, 1e-10, 0.0),
            "dhyper sums to {total} ({nr}, {nb}, {n})"
        );
        let exact = |x: u64| (lchoose(nr, x) + lchoose(nb, n - x) - lchoose(nr + nb, n)).exp();
        for (i, x) in (lo..=hi).enumerate() {
            ensure!(
                close(dens[i], exact(x), 1e-300, 1e-8),
                "dhyper({x}, {nr}, {nb}, {n}) = {}, want {}",
                dens[i],
                exact(x)
            );
        }
        let mut prev_lower = 0.0;
        for x in lo as i64 - 2..=hi as i64 + 1 {
            let lower = stats::phyper(x, nr, nb, n, true);
            let upper = stats::phyper(x, nr, nb, n, false);
            let want: f64 = (lo as i64..=x.min(hi as i64))
                .map(|y| dens[(y - lo as i64) as usize])
                .sum();
            ensure!(
                close(lower + upper, 1.0, 1e-12, 0.0),
                "phyper tails sum to {} at {x} ({nr}, {nb}, {n})",
                lower + upper
            );
            ensure!(
                close(lower, want, 1e-12, 1e-8) && lower >= prev_lower - 1e-15,
                "phyper({x}, {nr}, {nb}, {n}) = {lower}, sum {want}, previous {prev_lower}"
            );
            prev_lower = lower;
        }
        // Fisher's test over the tables with these margins: a in the corner,
        // n in the first row and nr in the first column.
        let mut prev = 1.0 + 1e-15;
        for a in lo..=hi {
            let (b, c) = (n - a, nr - a);
            let d = nb - b;
            let f = stats::fisher_greater(a, b, c, d);
            let want: f64 = dens[(a - lo) as usize..].iter().sum();
            ensure!(
                close(f, want, 1e-12, 1e-8) && f <= prev && (0.0..=1.0).contains(&f),
                "fisher_greater({a}, {b}, {c}, {d}) = {f}, sum {want}, previous {prev}"
            );
            let lf = stats::ln_fisher_greater(a, b, c, d);
            ensure!(
                close(lf, f.ln(), 1e-9, 1e-9),
                "ln_fisher_greater {lf}, log of {f}"
            );
            prev = f;
        }
        Ok(())
    });
}

#[test]
fn prop_fisher_tiny_pvalues() {
    run("fisher_tiny", 100, |rng| {
        // Tables so lopsided that the P-value underflows a double.
        let a = rng.random_range(300..=2000u64);
        let d = rng.random_range(300..=2000u64);
        let b = rng.random_range(0..=5u64);
        let c = rng.random_range(0..=5u64);
        let lf = stats::ln_fisher_greater(a, b, c, d);
        let (r, n, k) = (a + c, a + b + c + d, a + b);
        let terms: Vec<f64> = (a..=r.min(k))
            .map(|x| lchoose(r, x) + lchoose(n - r, k - x) - lchoose(n, k))
            .collect();
        let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let want = max + terms.iter().map(|t| (t - max).exp()).sum::<f64>().ln();
        ensure!(
            close(lf, want, 1e-6, 1e-8),
            "ln_fisher_greater({a}, {b}, {c}, {d}) = {lf}, want {want}"
        );
        Ok(())
    });
}

#[test]
fn prop_p_adjust_bh() {
    run("p_adjust_bh", 300, |rng| {
        let n = rng.random_range(1..=60);
        let ties = rng.random_bool(0.3);
        let p: Vec<f64> = (0..n)
            .map(|_| {
                let v = 10f64.powf(rng.random_range(-12.0..0.0));
                if ties { (v * 10.0).ceil() / 10.0 } else { v }
            })
            .collect();
        let q = stats::p_adjust_bh(&p);
        let mut sorted = p.clone();
        sorted.sort_by(f64::total_cmp);
        for i in 0..n {
            ensure!(q[i] >= p[i] && q[i] <= 1.0, "q {} for p {}", q[i], p[i]);
            for j in 0..n {
                ensure!(
                    p[i] > p[j] || q[i] <= q[j],
                    "not monotone: p {} -> q {}, p {} -> q {}",
                    p[i],
                    q[i],
                    p[j],
                    q[j]
                );
            }
            // The smallest n p(j) / j over the ranks j at or above p[i].
            let want = (0..n)
                .filter(|&j| sorted[j] >= p[i])
                .map(|j| n as f64 * sorted[j] / (j + 1) as f64)
                .fold(1.0f64, f64::min);
            ensure!(
                close(q[i], want, 0.0, 1e-12),
                "q {} for p {}, want {want}",
                q[i],
                p[i]
            );
        }
        let ln: Vec<f64> = p.iter().map(|x| x.ln()).collect();
        let lq = stats::ln_p_adjust_bh(&ln, n as f64);
        for (a, b) in lq.iter().zip(&q) {
            ensure!(
                close(a.exp().min(1.0), *b, 1e-15, 1e-9),
                "ln_p_adjust_bh gives {}, p_adjust_bh {b}",
                a.exp()
            );
        }
        Ok(())
    });
}
