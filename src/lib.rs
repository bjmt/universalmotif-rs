//! universalmotif: sequence motif types, conversions and analysis.
//!
//! The central type is [`Motif`], a position matrix over an [`Alphabet`]
//! together with its metadata. Motifs are created with the functions in
//! [`create`] or with [`MotifBuilder`], converted between count, probability,
//! weight, information content and contribution matrices with
//! [`convert::convert_type`], and transformed with the functions in [`ops`].
//!
//! The behaviour follows the universalmotif R/Bioconductor package; see
//! `DIVERGENCES.md` for the places where it intentionally differs.

pub mod alphabet;
pub mod convert;
pub mod create;
pub mod error;
pub mod matrix;
pub mod merge;
pub mod motif;
pub mod ops;

pub mod cdf;
pub mod cluster;
pub mod compare;
pub mod complexity;
pub mod dedup;
pub mod finder;
pub mod io;
pub mod kmer;
pub mod pvalue;
pub mod sample;
pub mod scan;
pub mod shuffle;
pub mod stats;

pub use alphabet::Alphabet;
pub use error::{Error, Result};
pub use matrix::Matrix;
pub use motif::{GapInfo, Motif, MotifBuilder, MotifParts, MotifType, Strand};

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
