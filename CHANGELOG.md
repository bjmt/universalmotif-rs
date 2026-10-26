# Changelog

## Unreleased

- **Breaking:** `HitTable::motif` and `HitTable::sequence` are `Vec<u32>`
  rather than `Vec<i64>`, which takes 8 bytes off every hit; positions stay
  `i64`.
- Scanning keeps hit positions as 64-bit numbers throughout, so hits in a
  single sequence longer than 4 Gb no longer get wrapped positions.
- Large scans gather each thread's hits into a list of the final size and
  no longer copy the first batch, which lowers peak memory by about a
  quarter for scans with millions of hits.
- `scoring::get_matches()` rounds a score between grid points up, as
  scanning and P-values do, rather than returning words scoring below it.
- `complexity()`, `sequence_complexity()` and `fasta::write()` work by
  letter rather than by byte, so sequences with letters outside ASCII no
  longer panic or come out garbled.
- `get_bkg()` with reverse complements counts letters found only on the
  reverse strands.
- Markov shuffles leave a window shorter than k as it is.
- The median of no values is NaN rather than a panic.
- Genomes reject `.fai` records with no letters per line and BGZF blocks
  too small to be valid, rather than panicking.
- The error for letters outside a custom alphabet no longer says
  "alphabet alphabet".
- `score_match()`, `motif_range()`, `get_matches()` and P-values add a
  missing pseudocount the way scanning does, so a word always scores what
  a scan reports, and `normalize()` gives a PCM without one a pseudocount
  of 1. A CWM goes through its probabilities.
- Converting with a new pseudocount updates the consensus.
- A uniform column converts to an IC column of zeros, rather than one just
  below zero that fails validation, and back to a uniform column.
- Merging and `shuffle_motifs()` keep probability matrices whose columns
  each hold one letter as probabilities, rather than taking them for
  counts.
- TRANSFAC counts from a single site are read as they are.
- `compare_motifs()` and `compare_matrix()` with P-values gather each
  query's q-values in one pass, which makes thousands of motifs several
  times faster, and scanning prepares one motif per thread at a time, so
  that long P-value tables stay within its memory budget.
- Comparisons centre each column once and build their consensus strings in
  parallel, which makes them up to 3 times faster on one thread and 6
  times faster on four for amino acid motifs, and P-value tables are built
  in cache-sized blocks, which makes amino acid P-values about 1.5 times
  faster. Results are unchanged.
- Enrichment builds each motif's P-value table once for the targets and
  the background, DNA and RNA scans score both strands in one pass, and
  co-occurrence distances only visit the sequences both motifs hit, which
  makes enrichment about 1.6 times faster, DNA scans about 1.35 times
  faster and co-occurrence up to 1.9 times faster. Gapped scans skip
  starts whose next piece cannot reach the threshold over any gap size,
  which makes them about 1.25 times faster. Results are unchanged.

## 0.1.1 (2026-10-24)

- `create::from_typed_matrix()` builds a motif from a matrix of a known
  type rather than detecting it, such as log-odds filters with no negative
  values.
- `scoring::score_match()` adds up the rounded matrix entries as integers,
  as scanning does, so that a word scores exactly what a scan reports;
  adding them as floats could be 0.001 off.
- Genomes open FASTA files whose last line has no line ending.
- `create_sequences()` rejects frequencies that are all zero, and
  `motif_proximity()` an empty list of distances.
- Count matrices must have at least one count in every column, and a
  matrix read as counts because it holds whole numbers explains when its
  columns have different sums.
- Merging keeps the first of two motifs with the same IC as the anchor,
  such as a motif and its reverse complement, whatever the rounding.

## 0.1.0 (2026-10-23)

The first release, a port of the universalmotif R/Bioconductor package
covering most of what it does:

- The `Motif` type, with count, probability, weight, information content
  and contribution matrices over DNA, RNA, amino acid or custom alphabets,
  and conversions between them.
- Creating motifs from sites, consensus strings or matrices, and trimming,
  rounding, merging, reverse complementing, gapping and switching the
  alphabet of motifs.
- Readers and writers for MEME, JASPAR, HOMER, TRANSFAC, CIS-BP, UniPROBE
  and plain matrix files, and for FASTA files, gzipped or not.
- Exact P-values and score thresholds, and scores and probabilities of
  single matches.
- Scanning sequences on both strands for any alphabet, including gapped
  motifs, sequences packed at two bits per letter, FASTA files read a batch
  at a time, and sparse matrices of motif matches for single-cell ATAC
  data.
- Motif comparison with P-values, scored on the overlap or on the whole
  query, merging of similar motifs, and hierarchical clustering into
  trees.
- Enrichment, co-occurrence, positional enrichment, proximity to genomic
  anchors, composition-matched backgrounds and *de novo* motif discovery.
- Genomes read a piece at a time from plain or bgzipped FASTA files and 2bit
  files.
- Sequence utilities: k-mer counts and backgrounds, shuffling, random
  sequences, complexity, masking, sliding windows, and sampling and
  implanting sites.

Positions count from 0 and leave the end out, as in BED files and Rust
ranges. Where the crate differs from the R package on purpose,
[DIVERGENCES.md](DIVERGENCES.md) says how.
