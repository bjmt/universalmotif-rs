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
- **Contribution matrices.** CWMs are compared by the cosine similarity
  of their signed contributions by default, or by TF-MoDISco's continuous
  Jaccard similarity, so that a negative pattern does not match the
  positive one with the same letters, and they are merged and clustered
  as CWMs.
- **Hits in contribution scores.** `scan_contributions()` finds where a
  sequence model relied on each CWM, matching the CWMs against the model's
  contribution scores rather than against the sequence, with thresholds
  from each TF-MoDISco pattern's seqlets in the Python package. With
  hypothetical scores it can also match database motifs, through their
  log-odds, against the letters the model wants at each position.
- **Seqlets.** `call_seqlets()` calls the stretches of a model's
  contribution scores that stand out from a null, either the scores
  shuffled or the model's scores for shuffled sequences, at a false
  discovery rate. Batches called with the same null give the seqlets of
  one call, so scores need not be kept once called.
- **Choosing and refining motifs.** `select_motifs()` keeps the motifs
  that reach enough sequences or sites, such as seqlets, not reached by
  those before them, leaving out motifs whose CWMs repeat one kept.
  `refine_motifs()` rebuilds motifs from their hits, widened and trimmed.
- **Motifs from sequence models.** The R package has CWMs and
  `trim_cwm()`. The Python package also reads the contribution scores of
  chromBPNet, BPNet and other models (`read_contributions()`) and
  TF-MoDISco's results (`read_modisco()`, `read_seqlets()`,
  `modisco_thresholds()`), and draws TF-MoDISco patterns as a report with
  their closest database motifs (`plot_modisco()`, also as an HTML page).
  `contribution_motifs()`, experimental for now, finds a model's motifs
  from its scores in one call, from seqlets to patterns and hits, and
  `combine_patterns()` gathers the patterns of several runs into one set.
- **Experiments on a model.** The Python package makes the inputs of
  common experiments on sequence models: one-hot batches, sites inserted
  into background sequences, pairs of sites at a range of spacings,
  mutagenesis of single letters and of pairs, shuffled or deleted regions,
  and the shuffled references DeepLIFT and DeepSHAP take, with plots of
  the results. `to_filters()` and `from_filters()` turn motifs into
  convolution filters and back.
- **Variants.** The Python package scores how variants change motif sites,
  insertions and deletions included, along a sequence or across a genome
  (`variant_motif_effects()`), reads VCF files (`read_vcf()`), and makes
  batches of both alleles for a model (`variant_batches()`).
- **Discovery.** `motif_finder()` can rerun its search on random halves of
  the sequences (`stability`) and report how often each motif is found
  again. `gapped_motif_finder()`, experimental for now, finds motifs made
  of two parts with a spacer of varying length, testing each in sequences
  held out from the search.
- **First hits only.** Scans can keep only the first hit of each motif in
  each sequence and stop scanning that sequence for the motif after it.
- **Package-wide threads.** The Python package's `set_nthreads()`, or the
  `UNIVERSALMOTIF_NTHREADS` environment variable, sets the threads used by
  every function whose `nthreads` is left out, where the R package takes
  `nthreads` per call only.
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
  with GenomicRanges instead. `read_bed()` and `write_bed()` read and
  write BED and narrowPeak files, `resize_ranges()` gives ranges one width
  around their centres or summits, and `drop_overlaps()` keeps only those
  that do not overlap one kept before.
- **Conversions.** The R package converts motifs to and from other
  Bioconductor packages. The Python package converts them to and from
  Biopython's motifs and to logomaker's matrices.
- **A command line.** The Python package installs `universalmotif`, also
  available as `umotif`, with subcommands that run its main functions on
  motif, FASTA, genome, BED, VCF and contribution score files, and write
  motifs, FASTA, tab-separated tables, figures or HTML reports.
- **Progress bars.** Long calls of the Python package, such as scanning a
  genome or discovering motifs, show a progress bar once they have run for
  a second (`set_progress()`), and Ctrl-C stops them.
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

- **Underscores.** Names use underscores where the R package uses dots,
  so `min.ic` becomes `min_ic`.
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
- **The R package's format.** `read_motifs()` and `write_motifs()`, which
  read and write the R package's own YAML-based format, are
  `read_universalmotif()` and `write_universalmotif()`. Higher-order
  matrices, backgrounds of longer words and the extra P-values of gaps are
  left out when reading, as the port does not hold them.

## Defaults

These defaults differ from the R package's, to suit most uses better.

- **FASTA names.** Sequences read from FASTA files are named by the first
  word of the header, as genomes on disk are and as most tools do, and a
  name that repeats is an error; `trim_names = False` keeps whole headers.
- **Overlapping hits.** Removing overlapping hits, with `no_overlaps` or
  the Python package's `dedup_hits()`, compares hits on both strands by
  default, so that a palindromic site counts once, where the R package's
  `no.overlaps` works on each strand apart.
- **Strands.** Scans respect each motif's strand, so a motif for the `+`
  or `-` strand is only searched on that strand, where the R package
  searches both strands whenever `RC = TRUE`. The Python package scans RNA
  on one strand by default, and MEME files without a `strands:` line give
  DNA and RNA motifs for both strands, where the R package assumes `+`.
  Motifs found by `motif_finder()` on one strand are for that strand.
- **Shuffling.** `shuffle_sequences()` keeps dinucleotides by default
  (`k = 2`), as backgrounds usually should, where the R package keeps
  single letters.
- **Q-value filtering.** The Python package's `enrich_motifs()`,
  `motif_peaks()` and `motif_proximity()` return every motif unless
  `qvalue` is given, where the R package keeps only those with q-values up
  to 0.1 by default.
- **Discovery backgrounds.** Without background sequences,
  `motif_finder()` and `gapped_motif_finder()` compare the sequences with
  three shuffles of themselves, where the R package's `motif_finder()` uses
  one, so that motifs in a few percent of the sequences are found whatever
  the seed.
- **Faint letters.** Logos draw every letter, however short, where the R
  package leaves out those under 1% of the logo's height.
- **Logos of contribution matrices.** The Python package's `view_motifs()`
  draws CWMs as contributions by default, and other motifs as information
  content, where the R package draws every motif as information content
  unless told otherwise.

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
  without matching them point for point. Tree tips are ordered as
  `hclust()` orders them rather than ladderized as ggtree does, and the
  unrooted `radial`, `equal_angle` and `daylight` layouts are not
  available. Plots are not part of the Rust crate.
- `motif_pvalue(method = "exhaustive")`: the branch-and-bound P-value
  search and the random score search. Only the exact dynamic method is
  ported.
