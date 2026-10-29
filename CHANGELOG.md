# Changelog

## Unreleased

Breaking changes:

- `ScanOptions` has a new `first_hit` field, and `CompareOptions` new
  `method` and `overhang` fields, so struct literals need them or
  `..Default::default()`.

Other changes:

- `ScanOptions::first_hit` keeps only the leftmost hit of each motif in
  each sequence and stops scanning that sequence for the motif soon after
  it.
- `motif_matches()` with `MatchValue::Matches` stops scanning a sequence
  for a motif after its first hit, which is up to 1.7 times faster when
  most sequences contain the motifs.
- `CompareOptions::method` compares motifs by Euclidean distance
  (`ColumnMetric::Eucl`) as well as Pearson correlation, with exact
  P-values from either null.
- `CompareOptions::overhang` with `Overhang::Background` compares columns
  hanging over either motif's ends with the background, so that every
  alignment counts the same number of column scores.
- Comparison P-values round column correlations to 201 bins instead of
  51, which makes them more accurate, especially for short, near-perfect
  alignments, at about three times the cost of the P-values. P-values and
  q-values change, while scores and alignments do not.
- Comparisons are otherwise about 10% faster, with identical results.
- `complexity()` without an alphabet scores DNA against `ACGT` and RNA
  against `ACGU`, as the Python package and the R package do, rather than
  against the letters of the sequence.
- `compare::anchor_layout()` lowers `min_overlap` to the width of the
  narrower motif, as comparisons do, so motifs narrower than it are
  aligned rather than left at the anchor's edge.

## 0.2.0 (2026-10-26)

Breaking changes:

- `HitTable::motif` and `HitTable::sequence` are `Vec<u32>` rather than
  `Vec<i64>`, which takes 8 bytes off every hit; positions stay `i64`.
- `ops::normalize()` gives a PCM without a pseudocount a pseudocount of 1,
  rather than returning it unchanged.

Fixes:

- `scoring::score_match()`, `motif_range()`, `get_matches()` and P-values
  add a missing pseudocount the way scanning does, so a word always scores
  what a scan reports; a CWM goes through its probabilities.
- `scoring::get_matches()` rounds a score between grid points up, as
  scanning and P-values do, rather than returning words scoring below it.
- Scanning keeps hit positions as 64-bit numbers, so hits in a single
  sequence longer than 4 Gb no longer get wrapped positions.
- Converting with a new pseudocount updates the consensus.
- A uniform column converts to an IC column of zeros, rather than one just
  below zero that fails validation, and back to a uniform column.
- Merging and `shuffle_motifs()` keep probability matrices whose columns
  each hold one letter as probabilities, rather than taking them for
  counts.
- TRANSFAC counts from a single site are read as they are.
- `complexity()`, `sequence_complexity()` and `fasta::write()` work by
  letter rather than by byte, so sequences with letters outside ASCII no
  longer panic or come out garbled.
- `get_bkg()` with reverse complements counts letters found only on the
  reverse strands.
- Markov shuffles leave a window shorter than k as it is.
- Genomes reject `.fai` records with no letters per line and BGZF blocks
  too small to be valid, rather than panicking, and the median of no
  values is NaN rather than a panic.
- The error for letters outside a custom alphabet no longer says
  "alphabet alphabet".

Speed and memory, with results unchanged:

- Comparisons centre each column once, build their consensus strings in
  parallel and gather q-values in one pass: up to 3 times faster on one
  thread and 6 times on four for amino acid motifs, and far faster for
  thousands of motifs.
- DNA and RNA scans score both strands in one pass (about 1.35 times
  faster), gapped scans skip starts that cannot reach the threshold
  (about 1.25 times), and enrichment builds each P-value table once for
  targets and background (about 1.6 times).
- P-value tables are built in cache-sized blocks (amino acid P-values
  about 1.5 times faster), and `pvalue::min_pvalues()` is computed in
  closed form, in microseconds; its values agree with the old ones to
  about 1e-13.
- Large scans gather hits into lists of their final size, lowering peak
  memory by about a quarter for millions of hits, and keep long P-value
  tables within their memory budget.
- Co-occurrence distances only visit the sequences both motifs hit (up to
  1.9 times faster), and `io::meme::read()` takes linear rather than
  quadratic time in the number of motifs.

Additions:

- `kmer::get_bkg_table()` gives `get_bkg()`'s output as columns, with
  k-mers as indices into one list, in about a third of the memory.

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
