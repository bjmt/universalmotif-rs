//! Sequence motifs: creating, converting, reading and writing them, and
//! using them to scan, compare, cluster and test sequences.
//!
//! This crate is a port of the [universalmotif] R/Bioconductor package and
//! the engine behind its Python version; besides fixes for bugs in the R
//! package, [DIVERGENCES.md] lists where the two differ.
//!
//! [universalmotif]: https://bioconductor.org/packages/universalmotif
//! [DIVERGENCES.md]: https://github.com/bjmt/universalmotif-rs/blob/main/DIVERGENCES.md
//!
//! # The motif type
//!
//! A [`Motif`] is a position matrix over an [`Alphabet`] of DNA, RNA, amino
//! acids or a custom set of letters, together with its metadata: name,
//! background frequencies, pseudocount, site count, gaps and so on. Its
//! [`MotifType`] says whether the matrix holds counts, probabilities, weights,
//! information content or contributions. Motifs are immutable once built and
//! always valid; [`MotifBuilder`] and [`MotifParts`] build them.
//!
//! # A first example
//!
//! ```
//! use universalmotif::create::{CreateOptions, from_strings};
//! use universalmotif::compare::{CompareOptions, compare_motifs};
//! use universalmotif::io::meme;
//! use universalmotif::scan::{ScanOptions, scan_table};
//!
//! // Two motifs from aligned sites.
//! let opts = CreateOptions { nsites: Some(20.0), pseudocount: 1.0, ..Default::default() };
//! let ap1 = from_strings(&["TGACTCA", "TGAGTCA", "TGACTAA"], None, &opts)?;
//! let creb = from_strings(&["TGACGTCA", "TGACGTCA", "TTACGTAA"], None, &opts)?;
//!
//! // Write them as MEME text and read them back.
//! let text = meme::write(&[ap1, creb], &meme::WriteOptions::default())?;
//! let motifs = meme::read(&text, &meme::ReadOptions::default())?;
//!
//! // Scan both strands of a sequence at a P-value cutoff of 1e-3.
//! let seqs = ["CCCCTGACTCACCCCCATGACGTCATCCC"];
//! let hits = scan_table(&motifs, &seqs, &ScanOptions { pvalue: 1e-3, ..Default::default() })?;
//! assert!(hits.len() >= 2);
//! println!("{} at {}..{}", hits.matched(0, &seqs), hits.start[0], hits.end[0]);
//!
//! // Compare the first motif with both.
//! let comparisons = compare_motifs(&motifs, &[0], &CompareOptions::default())?;
//! assert_eq!(comparisons.len(), 2);
//! # Ok::<(), universalmotif::Error>(())
//! ```
//!
//! # Where things are
//!
//! | Task | Modules |
//! |---|---|
//! | Build and convert motifs | [`create`], [`convert`], [`ops`], [`merge`] |
//! | Read and write motif and FASTA files | [`io`] |
//! | Read genomes a piece at a time from FASTA, bgzip and 2bit files | [`genome`] |
//! | P-values and match scores | [`pvalue`], [`scoring`], [`cdf`] |
//! | Scan sequences, including gapped motifs and whole genomes | [`scan`], [`packed`] |
//! | Compare and cluster motifs | [`compare`], [`cluster`], [`metrics`] |
//! | Enrichment, co-occurrence and positional enrichment | [`enrich`], [`coocc`], [`peaks`] |
//! | Genomic anchors and matched backgrounds | [`proximity`], [`background`], [`ranges`] |
//! | Sequences: k-mers, complexity, shuffling, sampling, masking | [`kmer`], [`complexity`], [`shuffle`], [`sample`], [`mask`] |
//! | De novo motif discovery | [`finder`] |
//!
//! Functions that do a lot of work take a thread count, `nthreads`, in their
//! options, where 0 means all cores.
//!
//! Options are structs with public fields and a `Default`, so a call sets only
//! the fields it needs and takes the rest from `..Default::default()`. Until
//! version 1.0, options may gain fields in new minor versions. Results, such
//! as [`scan::Hit`] or [`compare::Comparison`], are marked `#[non_exhaustive]`,
//! so they can gain fields without breaking code that reads them.

#![warn(missing_docs)]

pub mod alphabet;
pub mod convert;
pub mod create;
pub mod error;
pub mod matrix;
pub mod merge;
pub mod motif;
pub mod ops;
pub mod packed;

pub mod background;
pub mod cdf;
pub mod cluster;
pub mod compare;
pub mod complexity;
pub mod coocc;
pub mod dedup;
pub mod enrich;
pub mod finder;
pub mod genome;
pub mod io;
pub mod kmer;
pub mod mask;
pub mod metrics;
pub mod peaks;
pub mod proximity;
pub mod pvalue;
pub mod ranges;
pub mod sample;
pub mod scan;
pub mod scoring;
pub mod shuffle;
#[doc(hidden)]
pub mod stats;

pub use alphabet::Alphabet;
pub use error::{Error, Result};
pub use matrix::Matrix;
pub use motif::{GapInfo, Motif, MotifBuilder, MotifParts, MotifType, Strand};

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
