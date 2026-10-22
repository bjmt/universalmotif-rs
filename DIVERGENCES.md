# Differences from the R package

This crate, and the Python package built on it, follow the behaviour of the
universalmotif R package. The test suite checks the results against
reference values exported from R (see `scripts/make_reference.R`). Where the
R package has a bug, the port fixes it rather than copying it, and the
change is listed here so the R package can be fixed to match.

## Fixed in the port

### Motif objects

- **`trim_motifs()` invents `nsites`.** The R package writes `nsites = 100`
  into any trimmed motif that had no `nsites`, because it fills the slot in
  before computing the per-column IC. The port uses 100 only for that
  computation and leaves `nsites` empty.
- **`create_motif()` fails on single-column matrices with row names.** The
  rows are sorted with `matrix[order(rownames(matrix)), ]` without
  `drop = FALSE`, so a one-column matrix collapses to a vector, loses its row
  names, and the call errors. The port handles these matrices.
- **`motif_rc()` drops gap information.** The motif is rebuilt without its
  `gapinfo` slot. The port mirrors the gap positions onto the reverse
  complement.
- **Trimming leaves gaps in the wrong place.** `trim_motifs()` and
  `trim_cwm()` remove columns but keep the old gap positions. The port shifts
  them and drops any gap that no longer falls inside the motif.
- **Gap positions are only checked when there are several gaps.**
  `validObject_universalmotif()` skips the range check on `gaploc` when it has
  length 1. The port always checks it, and also enforces the
  `mingap <= maxgap` and `maxgap > 0` rules that the R package only applies
  in `add_gap()`.
- **Some invalid values pass validation.** The port also rejects matrices
  containing NaN and pseudocounts that are negative or not finite.
- **Letters outside a custom alphabet are ignored.** When `create_motif()` is
  given several strings and an explicit custom alphabet, letters not in the
  alphabet are silently dropped from the counts. The port reports an error.
- **Any 20-row matrix is taken to be amino acids.** When guessing the
  alphabet of a labelled matrix, the R package treats every 20-row matrix as
  AA regardless of its row names. The port requires the row names to be the
  twenty amino acids and otherwise makes a custom alphabet.

### Motif files

- **TRANSFAC files with a header block cannot be read.** `read_transfac()`
  treats the `VV ... //` header at the top of TRANSFAC releases as a motif
  and fails because it has no matrix. The port skips records without a matrix.
- **Only the first word of a field is kept.** `read_transfac()` keeps only
  the first word of `HC` and `OS` (so `OS Homo sapiens` gives `Homo`), and
  `read_cisbp()` does the same for `TF Name`, `Family` and `Species`. The port
  keeps the whole value.
- **TRANSFAC names depend on field order.** When the `ID` line comes before
  the `AC` line, `read_transfac()` uses the identifier as both the name and
  the alternate name. The port always uses the accession as the name.
- **TRANSFAC and JASPAR rows are assumed to be in alphabetical order.**
  `read_transfac()` ignores the letters on the `P0` line and `read_jaspar()`
  ignores the order of the row labels. The port reorders the values to match
  the labels.
- **Unusual JASPAR alphabets.** `read_jaspar()` treats any labelled matrix
  that is not RNA or protein as DNA, and any unlabelled matrix as DNA.
  The port makes a custom alphabet from unknown labels and reads unlabelled
  20-row matrices as protein.
- **HOMER names with spaces.** `read_homer()` splits header lines on any
  whitespace. The port splits tab-separated headers on tabs only.
- **MEME custom alphabet backgrounds.** `read_meme()` assigns background
  frequencies by position, which puts them on the wrong letters when MEME's
  letter order (letters before digits) differs from character order. The port
  matches them by letter.
- **MEME sites at positions of 10,000 or more.** In BLOCKS sections the
  position then fills its parentheses (`(12345)`), the columns shift, and
  `read_meme()` returns the wrong field as the site. The port finds the site
  after the closing parenthesis.
- **MEME block diagrams over three lines.** `read_meme()` joins only the
  first two lines of a combined block diagram. The port joins all of them.
- **Numbers are parsed with correct rounding.** R's parser can be one unit
  in the last place off (`as.numeric("5.3e-032")` is not the closest double
  to 5.3e-32), so values read by the two packages can differ in the last
  bit.

### Scanning and P-values

- **HOMER thresholds are on the wrong scale.** `write_homer()` computes the
  log-odds detection threshold from the motif's log2 PWM (with its own
  background and pseudocount). HOMER scores sites with natural logs against
  a fixed 0.25 background, flooring probabilities at 0.001
  (`Motif::logxform()` in HOMER's `cpp/Motif2.cpp`), so the thresholds the R
  package writes are about 44% too high on HOMER's scale. The port computes
  them with HOMER's own scoring of the probabilities as written. P-value
  thresholds use a uniform background, as HOMER assumes, and are raised by
  0.001 per position so that integer-grid rounding cannot let in sites above
  the P-value.
- **Short sequences stop the scan.** `scan_sequences_lite()` fails if any
  sequence is shorter than a motif. The port reports no hits for that motif
  and sequence and carries on.
- **Lower-case sequences.** The port scans lower-case letters like
  upper-case ones (Biostrings sequences are always upper case).

### Comparison

- **Tied alignments are decided by rounding.** When two alignments score
  the same, as both strands of a palindromic motif do,
  `compare_motifs_lite()` keeps whichever sums up larger in the last bit,
  which depends on the compiler (clang on arm64 fuses the multiply-adds in
  the column correlations, for example). Since the P-value depends on the
  alignment chosen, results can differ between machines. The port keeps
  the first alignment found (forward strand, leftmost offset) unless a
  later one is better by more than 1e-10. This also affects which motifs
  `view_motifs_lite()` shows as reverse complements.
- **RNA consensus strings use U.** `compare_motifs()` writes T in the
  consensus strings of aligned RNA motifs, unlike every other RNA consensus
  in the package. The port writes U.

### Plotting

- **Tip logos fail for motifs with the same name.** `motif_tree_lite()`
  labels the tree's tips with deduplicated names (`creb`, `creb.1`) but
  names the tip logos with the original ones, so drawing them stops with
  "`motif` value(s) not found in the supplied list". The Python port matches
  logos to tips by position.

### Sequences

- **Markov shuffling and sequence generation skip a letter of context.**
  `shuffle_sequences(method = "markov")` and `create_sequences()` with
  k-let frequencies draw each new letter given the k-1 letters before the
  previous one, rather than the k-1 letters just generated, so the output
  does not follow the k-let frequencies it was built from (for k = 2 the
  letter just written is ignored). The port conditions on the last k-1
  letters.
- **Wootton-Federhen complexity fails for windows of 171 letters or more.**
  `calc_complexity()` and `sequence_complexity()` compute the factorials in
  the formula directly, which overflow a double at 171! and give NaN. The
  port works with log-factorials.
- **k-lets with unknown letters can be counted.** `get_bkg()` and
  `count_klets()` mark letters outside the alphabet with a large negative
  number and skip k-lets whose index comes out negative. With large
  alphabets or long k-lets the valid letters can outweigh the marker (for
  amino acids from k = 5), and from k = 7 for DNA the index overflows, so
  such k-lets are added to the wrong counts. The port skips any k-let
  containing an unknown letter.
- **`calc_windows(return.incomp = FALSE)` drops a complete last window.**
  The last window is removed when `stop - start` is below the window size,
  which is also true of a full window (`calc_windows(9, 3, 0, FALSE)` loses
  7-9). The port removes it only when it is shorter than the window.
- **`motif_finder()` breaks ties between seeds in hash order.** Candidate
  k-mers with the same Fisher P-value are ranked in the iteration order of
  an `unordered_map`, so which seed is tried first depends on the standard
  library. The port ranks ties by k-mer. Results agree with the R package
  whenever the best seeds have distinct P-values.

### Enrichment and co-occurrence

- **`motif_coocc()` odds ratios overflow.** The 2x2 cells are R integers,
  so `both * neither` and `a_only * b_only` overflow once they pass 2^31,
  and the odds ratio becomes NA with a warning (for example with 60,000
  shared sequences out of 200,000). The port computes it in floating point.
- **`motif_proximity()` counts wide hits too often.** A hit is counted as
  near an anchor when any of it overlaps the near zone, but the null
  probability is the zone's share of the universe, as if hits were single
  positions. A 10 bp hit overlaps a 21 bp zone from 30 start positions, so
  near-zone counts run about 40% high at small distances and the binomial
  test reports enrichment for hits placed at random (p = 0.007 for 20,000
  random 10 bp hits and 1,000 random anchors, against p = 0.87 for the same
  hits as single positions). The distances behind the KS test have the
  same problem. The port places every hit at its middle position, which
  makes the counts, distances and nulls agree; distances are the number of
  positions from that middle to the anchor (0 inside it).
- **`motif_proximity()` never clips zones at chromosome ends.** The
  anchors are meant to take the universe's chromosome lengths so that
  `trim()` clips widened zones, but the universe built from
  `seqlengths(hits)` is created without lengths, so nothing is clipped.
  When counting anchors, zones near the ends then include positions off
  the chromosome, which inflates the expected count. The port clips
  wherever the lengths are known.
- **`match_bkg(genome = ...)` with one target draws random window
  widths.** The window widths come from `sample(target.widths, ...)`, and
  R's `sample()` treats a single number `n` as `1:n`, so with one target
  sequence of width 150 the windows get any width from 1 to 150. The port
  always uses the target widths.
- **`motif_peaks()` orders motifs by locale.** Motifs are numbered (and
  ties in P-value broken) by `sort()` on their names, which follows the
  locale's collation rules, so `motif.i` can differ between machines. The
  port sorts names by code point.

## Additions

- **FASTA files.** The port reads and writes FASTA files, since scanning
  needs sequences; the R package relies on Biostrings for this.
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

## Known issues shared with the R package

- **`ppm_to_pcm()` can produce negative counts.** The rounding error is added
  to the largest count, which goes negative when many small probabilities all
  round up (for example `create_motif("BZJX", nsites = 50)`). Both packages
  then reject the motif. Fixing this needs a different rounding scheme, which
  would change counts for ordinary motifs too, so it is left for now.

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
  polygons exactly (they are tested against them), while trees and peak
  plots follow the R package's layouts without matching them point for
  point. Tree tips are ordered as `hclust()` orders them rather than
  ladderized as ggtree does, and the unrooted `radial`, `equal_angle` and
  `daylight` layouts are not available. Plots are not part of the Rust
  crate.
- `motif_pvalue(method = "exhaustive")`: the branch-and-bound P-value
  search and the random score search. Only the exact dynamic method is
  ported.
