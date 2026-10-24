# Changelog

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
