# Differences from the R package

This crate, and the Python package built on it, follow the behaviour of the
universalmotif R package. The test suite checks the results against
reference values exported from R (see `scripts/make_reference.R`). Bugs in
the R package found while porting it are fixed in the port rather than
copied, and are not listed here. The other differences are below.

## Additions

- **FASTA files.** The port reads and writes FASTA files (reading gzipped
  ones too), since scanning needs sequences; the R package relies on
  Biostrings for this.
- **Scanning without holding sequences as text.** Scanning can read a FASTA
  file a batch at a time, or take sequences packed at two bits per base
  (a byte per letter for other alphabets), which hold a human genome in
  about 0.8 GB. Packed sequences keep every letter outside the alphabet but
  not the case of the alphabet's own letters.
- **CIS-BP and UniPROBE writers.** The R package can read these formats but
  not write them. The port writes both (UniPROBE in its minimal form), and its
  files read correctly with the R package's `read_cisbp()` and
  `read_uniprobe()`.
- **Other alphabets.** The R package's lite functions for scanning,
  comparison, merging, trees, enrichment and `view_motifs_lite()` only take
  DNA and RNA motifs. The port also takes amino acid and custom-alphabet
  motifs, with results that match the R package's full `scan_sequences()`
  and `compare_motifs()` (PCC). Only DNA and RNA have a reverse strand, so
  other alphabets are scanned and aligned on the forward strand alone.
  Protein letters are matched in either case and custom letters exactly.
  The parametric null for comparison P-values remains DNA and RNA only,
  since its Dirichlet grid assumes four letters; other alphabets use the
  empirical null.
- **Gapped scanning.** `scan_sequences_lite()` ignores gaps. The port scans
  gapped motifs as the full `scan_sequences()` does: every combination of
  gap sizes is tried, matched letters show the gaps as dots, and P-values
  are those of the motif without its gaps, so they are not adjusted for the
  number of combinations. The port can Bonferroni-correct them for it
  (`correct_gaps`), which the R package cannot. Each piece between gaps is
  scored once per
  position and the combinations are built from those scores, which is
  about twice as fast as scanning every combination separately.
- **Complete comparison scores.** With `CompareOptions::complete`,
  comparisons score whole queries rather than overlaps, as TOMTOM does by
  default: each query column is scored against the median of its
  correlations with the null columns, so that columns hanging over the end
  of the target count as typical ones, and the P-value is that of the best
  of all the alignments tried rather than a Bonferroni correction. The
  Python package's `compare_motifs(complete=True)` agrees with TOMTOM's
  P-values to a correlation of 0.999.
- **Motif match matrices.** `scan::motif_matches` summarizes the hits of
  each motif in each sequence as a sparse matrix of matches, hit counts or
  best scores without keeping the hits, for annotating the peaks of
  single-cell ATAC data; the Python package returns it as a SciPy matrix
  or AnnData object with `match_motifs()` and `add_motif_matches()`.
- **Genomes on disk.** `genome::Genome` reads any part of a genome from a
  FASTA file, plain or compressed with bgzip, or a 2bit file, through the
  `.fai` and `.gzi` indexes `samtools faidx` writes (or ones built when the
  file is opened). The Python package wraps it as `Genome`, a mapping that
  `scan_sequences()`, `motif_proximity()` and `match_bkg(genome=...)` read
  a chromosome or a window at a time, with `extract_ranges()` for the
  sequences of peaks; the R package needs Biostrings and BSgenome for this.
- **Matching on signal.** The Python package's `signal_stats()` summarizes
  bigWig signal over ranges, and `match_bkg()` can match the windows it
  draws from a genome to the targets on that signal, such as chromatin
  accessibility, as well as on GC content and length.
- **Genomic intervals.** The Python package's functions that take ranges
  also accept PyRanges objects and polars-bio frames, converting their
  0-based positions where needed, and `to_pyranges()` and
  `to_polars_bio()` hand tables back. `genomic_hits()` moves hits found in
  sequences extracted from a genome back to genomic coordinates, flipping
  the strand of hits in reverse complemented ranges; the R package works
  with GenomicRanges instead.
- **A command line.** The Python package installs `universalmotif`, also
  available as `umotif`, with subcommands that run its main functions on
  motif, FASTA, genome and BED files and write motifs, FASTA or
  tab-separated tables.
- **P-values on a log scale.** `motif_finder()` and `enrich_motifs()` also
  report `log10_pvalue` and `log10_qvalue`, which stay finite where the
  P-values themselves are too small for a double and come out as 0.
- **Warnings for results that come out empty or odd.** The Python
  package warns when a motif cannot reach the P-value cutoff of a scan or
  enrichment test and so can have no hits (naming the smallest P-value it
  can give), when DNA motifs are scanned against sequences that look like
  RNA or the other way round, when a string naming a file is scanned as a
  sequence, when a score is above a motif's highest possible score, when
  `merge_motifs()` leaves out motifs that cannot overlap the anchor, and
  when FASTA names repeat. The R package returns these results silently.
- **Lower case letters.** `score_match()`, `prob_match()` and `mask_seqs()`
  accept lower case DNA, RNA and amino acid letters (Biostrings sequences
  are always upper case); masking keeps the case of the letters it leaves.

## Names

- **Spelling.** Names use underscores where the R package uses dots, and
  Canadian spelling throughout, so `summarise_motifs()` is
  `summarize_motifs()`, and `motif_peaks()` reports `best_centre` and
  `centres` where the R package has `best.center` and `centers`.
- **Positions.** Positions on sequences and genomes count from 0 and
  leave the end out, as in BED files, Python slices and Rust ranges, where
  the R package counts from 1 and includes the end. A hit the R package
  reports at 11 to 20 starts at 10 and ends at 20 here, and indices such as
  `motif_i` count from 0 as well. This applies to hits, ranges, windows,
  implant positions, masked ranges, positional enrichment centres and MEME
  site positions, which are converted when read. Window tables call their
  last column `end` rather than `stop`.
- **k-mers.** The R package calls words of k letters k-lets, and the port
  calls them k-mers, the more common term, so `count_klets()` and
  `get_klets()` are `count_kmers()` and `get_kmers()`, and `get_bkg()`
  names its first column `kmer` rather than `klet`.

## Conventions

- **Total IC after a new pseudocount.** Converting a motif to another type
  with a new pseudocount gives the total IC of the input with that
  pseudocount, so that converting never changes it. The R package
  recomputes it from the new matrix as its new type, which differs for
  ICM outputs and for PWMs converted to PPMs.

## Not reproducible

- **Random sequences, shuffles and samples.** `shuffle_sequences()`,
  `create_sequences()`, `shuffle_motifs()`, `sample_sites()`,
  `implant_motifs()`, `match_bkg()`, `motif_finder()` and
  `enrich_motifs_lite()` without background sequences, and
  `motif_proximity(method = "permutation")` use
  the port's own random number generator, so a given seed gives different
  results from the R package. Results do not depend on the number of
  threads, which is not true of the R package's shuffles.
- **Random motifs.** `create_motif()` with a width draws from R's random
  number generator, which the port does not reproduce. The port also gives a
  random PPM a random `nsites` whenever none is supplied, whereas
  the R package only does so when the type was left at its default.

## Not supported by design

- The R package's own YAML format (`read_motifs()` and `write_motifs()`).
- Higher-order backgrounds and multifreq matrices.
- Schneider's small-sample correction for information content.
- The non-lite versions of functions that have a lite counterpart in
  the R package.
- Bioconductor containers. The port takes sequences as strings (or
  mappings from names), genomes as mappings from chromosome names to
  sequences instead of BSgenome objects, and genomic ranges as tables
  instead of GRanges.
- ggplot2 and ggtree. The Python package draws with matplotlib instead:
  `view_motifs()` and `view_logo()` reproduce the R package's letter
  polygons exactly (they are tested against them), apart from drawing
  letters below zero upside down by default (`flip_neg=True`) and stacking
  them as a mirror image of those above zero, with the largest furthest
  from zero, while trees and peak plots follow the R package's layouts
  without matching them point for point. Tree tips are ordered as `hclust()` orders them rather than
  ladderized as ggtree does, and the unrooted `radial`, `equal_angle` and
  `daylight` layouts are not available. Plots are not part of the Rust
  crate.
- `motif_pvalue(method = "exhaustive")`: the branch-and-bound P-value
  search and the random score search. Only the exact dynamic method is
  ported.
