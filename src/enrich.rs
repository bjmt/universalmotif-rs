//! Motif enrichment (`enrich_motifs_lite()`).
//!
//! Target and background sequences are scanned for hits, and each motif's
//! hits are compared between the two sets with a one-sided Fisher's exact
//! test, either on the number of sequences with a hit or on the number of
//! hits per scannable position.

use crate::error::{Error, Result, arg_err};
use crate::motif::Motif;
use crate::scan::{ScanOptions, scan_table};
use crate::shuffle::{ShuffleMethod, ShuffleOptions, shuffle_sequences};
use crate::stats::{fisher_greater, p_adjust_bh};

/// What [`enrich_motifs`] counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnrichTest {
    /// Sequences containing at least one hit.
    Seqs,
    /// Hits, out of all the positions a motif could occupy.
    Sites,
}

impl EnrichTest {
    pub fn parse(s: &str) -> Result<EnrichTest> {
        match s {
            "seqs" => Ok(EnrichTest::Seqs),
            "sites" => Ok(EnrichTest::Sites),
            _ => arg_err(format!("unknown test '{s}'; expected 'seqs' or 'sites'")),
        }
    }
}

/// Options for [`enrich_motifs`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnrichOptions {
    /// P-value cutoff for hits.
    pub pvalue: f64,
    pub test: EnrichTest,
    /// Scan both strands.
    pub rc: bool,
    /// k-lets kept when shuffling the targets into a background.
    pub shuffle_k: usize,
    /// Seed for that shuffle.
    pub seed: u64,
    /// Added to every cell of the 2x2 table.
    pub pseudocount: u64,
    pub nthreads: usize,
}

impl Default for EnrichOptions {
    fn default() -> Self {
        EnrichOptions {
            pvalue: 1e-4,
            test: EnrichTest::Seqs,
            rc: true,
            shuffle_k: 2,
            seed: 1,
            pseudocount: 1,
            nthreads: 1,
        }
    }
}

/// Enrichment of one motif.
#[derive(Debug, Clone, PartialEq)]
pub struct Enrichment {
    /// Index of the motif.
    pub motif: usize,
    pub target_seq_n: usize,
    /// Target sequences with at least one hit.
    pub target_seq_hits: usize,
    /// Hits in the target sequences.
    pub target_site_hits: usize,
    pub bkg_seq_n: usize,
    pub bkg_seq_hits: usize,
    pub bkg_site_hits: usize,
    /// Ratio of the hit rates, with a floor of half a hit on the background
    /// rate to avoid dividing by zero.
    pub enrichment: f64,
    pub log2_enrichment: f64,
    /// One-sided Fisher's exact P-value.
    pub pvalue: f64,
    /// Benjamini-Hochberg q-value across all motifs.
    pub qvalue: f64,
}

struct Counts {
    seqs: Vec<usize>,
    sites: Vec<usize>,
}

fn count_hits<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    seqs: &[S],
    opts: &ScanOptions,
) -> Result<Counts> {
    let hits = scan_table(motifs, seqs, opts)?;
    let mut sites = vec![0; motifs.len()];
    let mut nseqs = vec![0; motifs.len()];
    // Hits come ordered by motif and sequence, so a new sequence for a motif
    // shows up as a change from the previous row.
    let mut last: Option<(i64, i64)> = None;
    for (&m, &s) in hits.motif.iter().zip(&hits.sequence) {
        sites[m as usize] += 1;
        if last != Some((m, s)) {
            nseqs[m as usize] += 1;
            last = Some((m, s));
        }
    }
    Ok(Counts { seqs: nseqs, sites })
}

/// Test DNA or RNA motifs for enrichment in `sequences` relative to
/// `background` (by default the targets shuffled with the Euler method,
/// keeping `shuffle_k`-lets), as universalmotif's `enrich_motifs_lite()`.
/// Returns one result per motif, in input order and unfiltered.
pub fn enrich_motifs<S: AsRef<str> + Sync>(
    motifs: &[Motif],
    sequences: &[S],
    background: Option<&[S]>,
    opts: &EnrichOptions,
) -> Result<Vec<Enrichment>> {
    let Some(first) = motifs.first() else {
        return arg_err("no motifs given");
    };
    if motifs.iter().any(|m| m.alphabet != first.alphabet) {
        return arg_err("all motifs must share the same alphabet");
    }
    if !first.alphabet.is_nucleotide() {
        return Err(Error::Unsupported(
            "enrichment is only supported for DNA and RNA motifs".into(),
        ));
    }
    if !(opts.pvalue > 0.0 && opts.pvalue < 1.0) {
        return arg_err("pvalue must be in (0, 1)");
    }
    if sequences.is_empty() {
        return arg_err("no sequences given");
    }
    let shuffled: Vec<String>;
    let bkg: Vec<&str> = match background {
        Some(b) => b.iter().map(|s| s.as_ref()).collect(),
        None => {
            shuffled = shuffle_sequences(
                sequences,
                &ShuffleOptions {
                    k: opts.shuffle_k,
                    method: ShuffleMethod::Euler,
                    seed: opts.seed,
                    window: None,
                    nthreads: opts.nthreads,
                },
            )?;
            shuffled.iter().map(|s| s.as_str()).collect()
        }
    };
    if bkg.is_empty() {
        return arg_err("no background sequences");
    }
    let scan = ScanOptions {
        pvalue: opts.pvalue,
        rc: opts.rc,
        no_overlaps: None,
        nthreads: opts.nthreads,
    };
    let tgt = count_hits(motifs, sequences, &scan)?;
    let neg = count_hits(motifs, &bkg, &scan)?;
    let (n_tgt, n_bkg) = (sequences.len(), bkg.len());
    let tgt_lens: Vec<usize> = sequences.iter().map(|s| s.as_ref().len()).collect();
    let bkg_lens: Vec<usize> = bkg.iter().map(|s| s.len()).collect();
    let strands = if opts.rc { 2 } else { 1 };
    let positions = |lens: &[usize], w: usize| -> usize {
        lens.iter()
            .map(|&l| (l + 1).saturating_sub(w))
            .sum::<usize>()
            * strands
    };

    let pc = opts.pseudocount;
    let mut out: Vec<Enrichment> = Vec::with_capacity(motifs.len());
    for (i, m) in motifs.iter().enumerate() {
        let (a, b, c, d, eff) = match opts.test {
            EnrichTest::Seqs => {
                let (a, b) = (tgt.seqs[i], neg.seqs[i]);
                let pos_rate = a as f64 / n_tgt as f64;
                let neg_rate = b as f64 / n_bkg as f64;
                let neg_den = if neg_rate > 1e-9 {
                    neg_rate
                } else {
                    0.5 / n_bkg as f64
                };
                (a, b, n_tgt - a, n_bkg - b, pos_rate / neg_den)
            }
            EnrichTest::Sites => {
                let w = m.width();
                let (ppos, npos) = (positions(&tgt_lens, w), positions(&bkg_lens, w));
                let (a, b) = (tgt.sites[i], neg.sites[i]);
                let pos_rate = if ppos > 0 {
                    a as f64 / ppos as f64
                } else {
                    0.0
                };
                let neg_den = if npos > 0 && b > 0 {
                    b as f64 / npos as f64
                } else {
                    0.5 / npos.max(1) as f64
                };
                let eff = if pos_rate > 0.0 {
                    pos_rate / neg_den
                } else {
                    0.0
                };
                (a, b, ppos.saturating_sub(a), npos.saturating_sub(b), eff)
            }
        };
        let pvalue = fisher_greater(a as u64 + pc, b as u64 + pc, c as u64 + pc, d as u64 + pc);
        out.push(Enrichment {
            motif: i,
            target_seq_n: n_tgt,
            target_seq_hits: tgt.seqs[i],
            target_site_hits: tgt.sites[i],
            bkg_seq_n: n_bkg,
            bkg_seq_hits: neg.seqs[i],
            bkg_site_hits: neg.sites[i],
            enrichment: eff,
            log2_enrichment: if eff > 0.0 {
                eff.log2()
            } else {
                f64::NEG_INFINITY
            },
            pvalue,
            qvalue: 0.0,
        });
    }
    let q = p_adjust_bh(&out.iter().map(|e| e.pvalue).collect::<Vec<_>>());
    for (e, q) in out.iter_mut().zip(q) {
        e.qvalue = q;
    }
    Ok(out)
}
