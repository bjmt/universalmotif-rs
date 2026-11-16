# Changelog

## Unreleased

New:

- `progress::with_progress()` runs a call with a `progress::Progress`
  handle, to which scanning, `motif_finder()`, `gapped_motif_finder()`,
  `shuffle_sequences()` and the functions built on them report their
  progress, and through which they can be cancelled
  (`Error::Cancelled`).

## 0.4.0 (2026-11-15)

Breaking changes:

- `scan::scan_fasta()` takes a new argument, `trim_names`.
- `MergeSimilarOptions` and `MergedClusters` have new fields, so struct
  literals need them or `..Default::default()`.
- The enums `compare::Overhang`, `compare::NullModel`,
  `similarity::CwmMetric`, `merge::MergeLinkage`, `enrich::EnrichTest`,
  `proximity::ProximityMethod` and `scan::MatchValue`, and the results
  `pvalue::ScoringMatrix` and `sample::Implant`, are `#[non_exhaustive]`, so
  a `match` on one of the enums needs a wildcard arm.

Changes in behaviour:

- `ScanOptions::default()` Bonferroni-corrects gapped motifs' P-values for
  the gap sizes tried (`correct_gaps: true`), which let several times too
  many chance hits through before.
- `NoOverlaps::default()` compares hits on both strands, so a palindromic
  site is kept once.
- Scans respect each motif's strand: motifs for the `+` or `-` strand are
  searched on that strand only, whatever `rc` says, and enrichment's site
  test counts positions on the strands searched.
- `ShuffleOptions::default()` keeps dinucleotides (`k: 2`).
- MEME files without a `strands:` line give DNA and RNA motifs for both
  strands.
- `motif_finder()` with `rc: false` gives motifs for the `+` strand.
- `sample::sample_sites()` and `sample::implant_motifs()` give gapped
  motifs' sites a spacer in each gap, of a size drawn within its range,
  rather than leaving the gaps out; implants keep the sequence's letters
  in their spacers.
- `create::from_strings()` refuses dots and dashes among DNA, RNA or amino
  acid letters when guessing the alphabet, rather than making a custom
  alphabet of them.
- `merge::merge_motifs()` merges CWMs as CWMs, averaging their signed
  contributions over alignments found by cosine similarity, and
  `merge::merge_similar()` groups CWMs by cosine similarity reaching
  `MergeSimilarOptions::min_score` rather than by q-value.
- `merge::merge_motifs()` aligns motifs narrower than `min_overlap` over
  their whole width, as comparisons do, rather than leaving them out.
- `compare::compare_motifs()` and `similarity::compare_cwms()` leave a
  query's comparison with itself out of the Benjamini-Hochberg correction
  of its q-values, and give it its P-value as its q-value; its tiny
  P-value made the best matches' q-values about half what they should be.
- `compare::anchor_layout()` aligns CWMs by the cosine of their
  contributions, through the new `similarity::cwm_anchor_layout()`.
- `Genome` and `scan_fasta(trim_names = true)` name sequences by the first
  word of the header, and a repeated name is an error.
- Gap locations must be in increasing order; `ops::add_gap()` refuses
  others, which made printing the motif panic.
- `scoring::get_scores()` lists at most 16,777,216 scores, as for a DNA
  motif of 12 columns, where wider motifs ran out of memory.
- `shuffle::shuffle_sequences()` refuses a k whose k-mer counts would not
  fit in memory.

Other changes:

- `gapped::gapped_motif_finder()`, experimental for now, discovers motifs
  made of two parts with a spacer of varying length, testing each in
  sequences held out from the search, from pairs of k-mers counted at every
  spacing by `finder::spaced_pairs()`.
- `MergeSimilarOptions::linkage` with `MergeLinkage::Complete` joins only
  motifs that are all similar to each other, and `MergedClusters::warnings`
  names a cluster holding more than a quarter of the motifs.
- `fasta::first_word()` and `fasta::trim_names()`.
- `mask::mask_complexity()` and `mask::dust_intervals()` find low-complexity
  regions with symmetric DUST, giving the same regions as dustmasker, and
  mask them, along with lowercase letters if asked.
- `similarity::compare_cwms()` and `similarity::cwm_similarity_matrix()`
  compare CWMs by cosine similarity, with P-values, or by TF-MoDISco's
  continuous Jaccard similarity, keeping their signs unless
  `ignore_sign` is set.
- `contributions::scan_contributions()` finds motif hits in the
  contribution scores of a sequence model, by the similarity and importance
  of each window against CWMs, gapped CWMs included.
- `background::match_bkg_values()` matches on GC fractions and lengths
  rather than sequences, `background::composition()` gives the GC and N
  fractions of a sequence, and `genome::Genome::composition_many()` measures
  many windows without holding them.
- `shuffle::shuffle_codes()` makes several k-mer-preserving shuffles of
  each of a set of sequences given as letter codes, keeping unknown
  letters in place, with each shuffle seeded by its sequence's place in
  the whole set.
- Faster: comparisons with `complete: true` (from quadratic in the number
  of motifs to linear for a few queries), removing overlapping hits in
  dense clusters, positional enrichment, proximity permutations and
  deduplication in `finder::motif_finder()`.
- `scan::motif_matches()` scans motifs in batches whose P-value tails fit
  in about 32 MB when P-values choose among overlapping hits, where it held
  every motif's tails at once.

Fixes:

- `finder::motif_finder()` names its motifs in order of rank, where one
  dropped by its q-value left a number out.
- `finder::motif_finder()` and `enrich::enrich_motifs()` leave sequences
  too short to shuffle as they are in the backgrounds they shuffle, where
  one empty sequence made them fail.
- Windowed `shuffle::shuffle_sequences()` shuffles sequences shorter than
  the window whole, where it returned them unchanged.
- `shuffle::create_sequences()` refuses more than one frequency for a
  one-letter alphabet, where it looped forever.
- `genome::Genome` refuses FASTA files with a blank line inside a record,
  which it read with the bases after it shifted.

## 0.3.0 (2026-10-29)

Breaking changes:

- `ScanOptions` has a new `first_hit` field, and `CompareOptions` new
  `method` and `overhang` fields, so struct literals need them or
  `..Default::default()`.

Changes in behaviour:

- Comparison P-values round column correlations to 201 bins instead of
  51, which makes them more accurate, especially for short, near-perfect
  alignments, at about three times the cost of the P-values. P-values and
  q-values change, while scores and alignments do not.
- `complexity()` without an alphabet scores DNA against `ACGT` and RNA
  against `ACGU`, as the Python package and the R package do, rather than
  against the letters of the sequence.
- `compare::anchor_layout()` lowers `min_overlap` to the width of the
  narrower motif, as comparisons do, so motifs narrower than it are
  aligned rather than left at the anchor's edge.

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
- Comparisons are about 10% faster, with identical results.

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
