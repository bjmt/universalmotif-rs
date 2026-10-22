//! Benchmarks for the main functions at realistic sizes.
//!
//! ```sh
//! cargo build --release --example bench
//! /usr/bin/time -l target/release/examples/bench scan_many [threads]
//! ```
//!
//! Each run does one task and prints the setup and task times; measure peak
//! memory with `/usr/bin/time` around it. Motifs come from the MEME file in
//! `UM_BENCH_MOTIFS`, for which the HOCOMOCO v11 human set is a good choice, or
//! are otherwise 400 random motifs of widths 8 to 25. With
//! `UM_BENCH_ALPHABET=AA`, the random motifs and sequences are proteins
//! instead.

use std::time::Instant;

use rand::SeedableRng;
use rand::rngs::StdRng;
use universalmotif::background::{MatchOptions, match_bkg};
use universalmotif::cluster::{Method, TreeOptions, hclust, motif_tree};
use universalmotif::compare::{CompareOptions, compare_matrix};
use universalmotif::complexity::{ComplexityOptions, Method as Cx, sequence_complexity};
use universalmotif::convert::convert_type;
use universalmotif::coocc::{CooccHit, CooccOptions, motif_coocc};
use universalmotif::create::{CreateOptions, random};
use universalmotif::enrich::{EnrichOptions, enrich_motifs};
use universalmotif::finder::{FinderOptions, motif_finder};
use universalmotif::io::meme;
use universalmotif::kmer::{BkgOptions, WindowOptions, get_bkg};
use universalmotif::merge::{MergeSimilarOptions, merge_similar};
use universalmotif::ops::add_gap;
use universalmotif::pvalue::{PvalueOptions, motif_pvalue};
use universalmotif::scan::{ScanOptions, scan_sequences, scan_table};
use universalmotif::shuffle::{ShuffleMethod, ShuffleOptions, create_sequences, shuffle_sequences};
use universalmotif::{Alphabet, Matrix, Motif, MotifType};

fn bench_alphabet() -> Alphabet {
    match std::env::var("UM_BENCH_ALPHABET").as_deref() {
        Ok("AA") => Alphabet::Aa,
        _ => Alphabet::Dna,
    }
}

fn bench_motifs() -> (String, Vec<Motif>) {
    let motifs = match std::env::var("UM_BENCH_MOTIFS") {
        Ok(path) => {
            let text = std::fs::read_to_string(path).unwrap();
            meme::read(&text, &meme::ReadOptions::default()).unwrap()
        }
        Err(_) => {
            let mut rng = StdRng::seed_from_u64(1);
            let opts = CreateOptions {
                nsites: Some(50.0),
                pseudocount: 1.0,
                ..Default::default()
            };
            (0..400)
                .map(|i| random(8 + i % 18, bench_alphabet(), &opts, &mut rng).unwrap())
                .collect()
        }
    };
    let text = meme::write(&motifs, &meme::WriteOptions::default()).unwrap();
    (text, motifs)
}

fn seqs(n: usize, len: usize) -> Vec<String> {
    create_sequences(bench_alphabet().letters(), n, len, None, 1, 1).unwrap()
}

fn main() {
    let task = std::env::args().nth(1).expect("task");
    let threads: usize = std::env::args()
        .nth(2)
        .map(|s| s.parse().unwrap())
        .unwrap_or(1);
    let t0 = Instant::now();
    let (text, motifs) = bench_motifs();
    let mut setup = t0.elapsed().as_secs_f64();
    let timed = |label: &str, setup: f64, f: &mut dyn FnMut() -> String| {
        let t = Instant::now();
        let out = f();
        println!(
            "{label}\tsetup {setup:.2}s\ttask {:.3}s\t{out}",
            t.elapsed().as_secs_f64()
        );
    };
    match task.as_str() {
        "read_meme" => timed("read_meme x20", setup, &mut || {
            let mut n = 0;
            for _ in 0..20 {
                n += meme::read(&text, &meme::ReadOptions::default())
                    .unwrap()
                    .len();
            }
            format!("{n} motifs")
        }),
        "write_meme" => timed("write_meme x20", setup, &mut || {
            let mut n = 0;
            for _ in 0..20 {
                n += meme::write(&motifs, &meme::WriteOptions::default())
                    .unwrap()
                    .len();
            }
            format!("{n} bytes")
        }),
        "pvalue" => timed("motif_pvalue, 10 scores per motif", setup, &mut || {
            let q: Vec<(usize, f64)> = (0..motifs.len())
                .flat_map(|i| (0..10).map(move |k| (i, k as f64)))
                .collect();
            let opts = PvalueOptions {
                nthreads: threads,
                ..Default::default()
            };
            let p = motif_pvalue(&motifs, &q, &opts).unwrap();
            format!("{} p-values", p.len())
        }),
        "scan_many" | "scan_long" => {
            let t = Instant::now();
            let (ms, s) = if task == "scan_many" {
                (motifs.clone(), seqs(1000, 1000))
            } else {
                (motifs[..20].to_vec(), seqs(50, 1_000_000))
            };
            setup += t.elapsed().as_secs_f64();
            let label = if task == "scan_many" {
                "scan all motifs x 1 Mb (1000 x 1 kb)"
            } else {
                "scan 20 motifs x 50 Mb (50 x 1 Mb)"
            };
            timed(label, setup, &mut || {
                let opts = ScanOptions {
                    nthreads: threads,
                    ..Default::default()
                };
                format!("{} hits", scan_table(&ms, &s, &opts).unwrap().len())
            })
        }
        "scan_gapped" | "scan_gapped_expanded" => {
            // 20 motifs with a gap of 0 to 5 in the middle (six
            // combinations each), scanned with their gaps or as the 120
            // ungapped motifs they expand into.
            let t = Instant::now();
            let s = seqs(10, 1_000_000);
            let gapped: Vec<Motif> = motifs[..20]
                .iter()
                .map(|m| add_gap(m, None, &[0], &[5]).unwrap())
                .collect();
            let expanded: Vec<Motif> = gapped
                .iter()
                .flat_map(|m| {
                    let pwm = convert_type(m, MotifType::Pwm, None, false).unwrap();
                    let cols: Vec<Vec<f64>> = pwm.matrix.columns().map(<[f64]>::to_vec).collect();
                    let at = m.gap.gaploc[0];
                    (0..=5).map(move |g| {
                        let mut c = cols[..at].to_vec();
                        c.extend(std::iter::repeat_n(vec![0.0; 4], g));
                        c.extend_from_slice(&cols[at..]);
                        Motif::builder(Matrix::from_columns(&c).unwrap(), Alphabet::Dna)
                            .motif_type(Some(MotifType::Pwm))
                            .build()
                            .unwrap()
                    })
                })
                .collect();
            setup += t.elapsed().as_secs_f64();
            let (label, ms) = if task == "scan_gapped" {
                ("scan 20 gapped motifs (6 gap sizes) x 10 Mb", &gapped)
            } else {
                ("scan their 120 expansions x 10 Mb", &expanded)
            };
            timed(label, setup, &mut || {
                let opts = ScanOptions {
                    nthreads: threads,
                    ..Default::default()
                };
                format!("{} hits", scan_table(ms, &s, &opts).unwrap().len())
            })
        }
        "compare_scores" | "compare_pvalues" => timed(
            &format!("compare_matrix all x all ({task})"),
            setup,
            &mut || {
                let opts = CompareOptions {
                    nthreads: threads,
                    ..Default::default()
                };
                let m = compare_matrix(&motifs, &opts, task == "compare_pvalues").unwrap();
                format!("{} rows", m.scores.len())
            },
        ),
        "motif_tree" => timed("motif_tree, all motifs", setup, &mut || {
            let opts = TreeOptions {
                min_overlap: 6,
                rc: true,
                method: Method::Complete,
                nthreads: threads,
            };
            format!("{} merges", motif_tree(&motifs, &opts).unwrap().merge.len())
        }),
        "merge_similar" => timed("merge_similar, all motifs", setup, &mut || {
            let opts = MergeSimilarOptions {
                nthreads: threads,
                ..Default::default()
            };
            format!("{} merged", merge_similar(&motifs, &opts).unwrap().1.len())
        }),
        "hclust" => {
            let n = 3000;
            let d: Vec<Vec<f64>> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| ((i * 7919 + j * 104729) % 1000) as f64)
                        .collect()
                })
                .collect();
            timed("hclust 3000 complete", setup, &mut || {
                format!(
                    "{} merges",
                    hclust(&d, Method::Complete).unwrap().merge.len()
                )
            })
        }
        "finder" => {
            let t = Instant::now();
            let mut s = seqs(500, 500);
            for (i, x) in s.iter_mut().enumerate().take(300) {
                let p = 100 + (i * 37) % 300;
                x.replace_range(p..p + 10, "TTGACGTCAA");
            }
            setup += t.elapsed().as_secs_f64();
            timed("motif_finder 500 x 500 bp, widths 6-15", setup, &mut || {
                let opts = FinderOptions {
                    nthreads: threads,
                    ..Default::default()
                };
                format!("{} motifs", motif_finder(&s, None, &opts).unwrap().len())
            })
        }
        "shuffle_euler" | "shuffle_markov" | "shuffle_linear" => {
            let t = Instant::now();
            let s = seqs(100, 100_000);
            setup += t.elapsed().as_secs_f64();
            let method = match task.as_str() {
                "shuffle_euler" => ShuffleMethod::Euler,
                "shuffle_markov" => ShuffleMethod::Markov,
                _ => ShuffleMethod::Linear,
            };
            timed(&format!("{task} k=2, 10 Mb"), setup, &mut || {
                let opts = ShuffleOptions {
                    k: 2,
                    method,
                    nthreads: threads,
                    ..Default::default()
                };
                format!("{} seqs", shuffle_sequences(&s, &opts).unwrap().len())
            })
        }
        "create" => timed("create_sequences 10 Mb", setup, &mut || {
            format!(
                "{} seqs",
                create_sequences(Alphabet::Dna.letters(), 100, 100_000, None, 1, threads)
                    .unwrap()
                    .len()
            )
        }),
        "bkg" | "bkg_window" => {
            let t = Instant::now();
            let s = if task == "bkg" {
                seqs(100, 100_000)
            } else {
                seqs(100, 10_000)
            };
            setup += t.elapsed().as_secs_f64();
            let names: Vec<String> = (0..s.len()).map(|i| format!("s{i}")).collect();
            let label = if task == "bkg" {
                "get_bkg k=1:3 merged, 10 Mb"
            } else {
                "get_bkg k=1:3 per sequence, 100 bp windows, 1 Mb"
            };
            timed(label, setup, &mut || {
                let opts = BkgOptions {
                    merge: task == "bkg",
                    window: (task == "bkg_window").then_some(WindowOptions {
                        size: 100.0,
                        overlap: 0.0,
                    }),
                    nthreads: threads,
                    ..Default::default()
                };
                format!("{} rows", get_bkg(&s, &names, &opts).unwrap().len())
            })
        }
        "complexity_wf" | "complexity_dust" | "complexity_trifonov" => {
            let t = Instant::now();
            let s = seqs(100, 10_000);
            setup += t.elapsed().as_secs_f64();
            let method = match task.as_str() {
                "complexity_wf" => Cx::WoottonFederhen,
                "complexity_dust" => Cx::Dust,
                _ => Cx::Trifonov,
            };
            timed(
                &format!("sequence_complexity {task}, 1 Mb, window 20"),
                setup,
                &mut || {
                    let opts = ComplexityOptions {
                        method,
                        nthreads: threads,
                        ..Default::default()
                    };
                    format!("{} windows", sequence_complexity(&s, &opts).unwrap().len())
                },
            )
        }
        "enrich" => {
            let t = Instant::now();
            let s = seqs(1000, 500);
            setup += t.elapsed().as_secs_f64();
            timed(
                "enrich_motifs, 1000 x 500 bp, shuffled background",
                setup,
                &mut || {
                    let opts = EnrichOptions {
                        nthreads: threads,
                        ..Default::default()
                    };
                    format!(
                        "{} rows",
                        enrich_motifs(&motifs, &s, None, &opts).unwrap().len()
                    )
                },
            )
        }
        "coocc" | "coocc_nodist" => {
            let t = Instant::now();
            let s = seqs(2000, 500);
            let hits: Vec<CooccHit> = scan_sequences(
                &motifs,
                &s,
                &ScanOptions {
                    pvalue: 1e-3,
                    nthreads: threads,
                    ..Default::default()
                },
            )
            .unwrap()
            .iter()
            .map(|h| CooccHit {
                motif: h.motif,
                sequence: h.sequence,
                start: h.start as i64,
            })
            .collect();
            setup += t.elapsed().as_secs_f64();
            let label = format!(
                "motif_coocc, {} hits, {}",
                hits.len(),
                if task == "coocc" {
                    "max_distance 50"
                } else {
                    "no distances"
                }
            );
            timed(&label, setup, &mut || {
                let opts = CooccOptions {
                    max_distance: (task == "coocc").then_some(50),
                    nthreads: threads,
                    ..Default::default()
                };
                format!(
                    "{} pairs",
                    motif_coocc(&hits, motifs.len(), s.len(), &opts)
                        .unwrap()
                        .len()
                )
            })
        }
        "match_bkg" => {
            let t = Instant::now();
            let targets = seqs(2000, 300);
            let mut universe = Vec::new();
            for (i, x) in seqs(50_000, 400).into_iter().enumerate() {
                universe.push(x[..150 + (i * 13) % 250].to_string());
            }
            setup += t.elapsed().as_secs_f64();
            timed(
                "match_bkg 2000 targets x 50,000 universe, 3 each",
                setup,
                &mut || {
                    let opts = MatchOptions {
                        n_per_target: 3,
                        ..Default::default()
                    };
                    format!(
                        "{} draws",
                        match_bkg(&targets, &universe, &[], &opts)
                            .unwrap()
                            .universe_index
                            .len()
                    )
                },
            )
        }
        _ => panic!("unknown task {task}"),
    }
}
