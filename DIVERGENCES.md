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
- **`ppm_to_pcm()` can make uneven or negative counts.** Each value is
  rounded and the whole rounding error goes to the largest count, so a
  uniform column at 50 sites becomes 11, 13, 13 and 13, and when many small
  probabilities round up the largest count goes negative and the motif is
  rejected (for example `create_motif("BZJX", nsites = 50)`). The port keeps
  this rule when the rounded counts fall short or overshoot by one, which
  covers ordinary DNA columns, and otherwise rounds down and gives the
  remaining sites to the largest remainders, so the uniform column becomes
  13, 13, 12 and 12. Amino acid columns overshoot more often, so their counts,
  and the PWMs made from them, differ from the R package's more often too.
- **The total IC and consensus go stale.** `round_motif()` and `normalize()`
  change the matrix or pseudocount but keep the old total IC and consensus,
  and `convert_type()` with a new pseudocount keeps the old total IC. The
  port recomputes them.
- **Relative entropy is NaN for zero probabilities.**
  `ppm_to_icm(relative_entropy = TRUE)` computes `0 * log2(0)` for letters
  that never occur, which is NaN, so the column, and a logo drawn with
  `view_motifs(relative_entropy = TRUE)`, comes out empty. The port takes
  those letters as 0.
- **CWMs get a total IC.** Contribution weights are not probabilities, so
  the IC the R package reports for a CWM has no meaning. The port gives CWMs
  none (NaN) and leaves it out when printing them, and computes one when a
  CWM is converted to another type.
- **A matrix with negative values becomes a near-uniform PPM.**
  `create_motif()` takes a matrix with negative values to be a PWM and,
  with the default `type = "PPM"`, converts it, so a matrix of contribution
  weights silently becomes a PPM of nearly uniform columns. The Python
  package asks for the type of such a matrix, `"PWM"` or `"CWM"`, instead.
- **Some invalid values pass validation.** The port also rejects matrices
  containing NaN and pseudocounts that are negative or not finite.
- **Letters outside a custom alphabet are ignored.** When `create_motif()` is
  given several strings and an explicit custom alphabet, letters not in the
  alphabet are silently dropped from the counts. The port reports an error.
- **Any 20-row matrix is taken to be amino acids.** When guessing the
  alphabet of a labelled matrix, the R package treats every 20-row matrix as
  AA regardless of its row names. The port requires the row names to be the
  twenty amino acids and otherwise makes a custom alphabet.

- **`motif_range(allow.nonfinite = TRUE)` returns NA.** With `-Inf`
  entries in the PWM, both ends of the range come back as NA. The port gives
  `-Inf` for the lowest score and the true highest score.
- **`log_string_pval()` fails for the values it exists for.** It parses the
  string as a `long double`, which is no wider than a double on some
  platforms (arm64 macOS among them), so P-values below about 1e-308, such
  as `"1e-400"`, raise an error there. The port works from the mantissa and
  exponent and handles any value.

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
- **MEME files get a placeholder E-value.** `write_meme()` writes `E= 0`
  for motifs without an E-value, which reads back as an E-value of 0, and
  leaves trailing spaces on the MOTIF and matrix lines. The port leaves
  out the optional E-value when there is none and the trailing spaces.
- **`meme_alph()` declares each complementary pair twice.** Every core
  letter gets its own `X ~ Y` line, so A ~ T is followed later by T ~ A.
  MEME's own alphabet files declare each pair once, as the port does.
- **Rounded probabilities are taken for information content.** A matrix
  of probabilities written to two decimals can have columns summing to
  0.98 or 1.02, and guessing its type, as `read_matrix()` and
  `create_motif()` do, only accepts sums between 0.99 and 1.01, so such a
  matrix is read as an ICM. The port accepts sums within 0.02 of 1 when
  guessing and renormalises the columns.
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
- **Gapped motifs mostly fail to scan.** The full `scan_sequences()` fails
  for gapped motifs with `RC = TRUE` (its default), for a motif with a
  single gap that can have size zero, and when an ungapped motif comes
  before a gapped one in the list. The port scans both strands and any
  gaps.

### Comparison

- **Tied alignments are decided by rounding.** When two alignments score
  the same, as both strands of a palindromic motif do,
  `compare_motifs_lite()` keeps whichever sums up larger in the last bit,
  which depends on the compiler (clang on arm64 fuses the multiply-adds in
  the column correlations, for example). Since the P-value depends on the
  alignment chosen, results can differ between machines. The port keeps
  the first alignment found (forward strand, leftmost offset) unless a
  later one is better by more than 1e-10, or ties with it over a longer
  overlap. Preferring the longer overlap keeps the mean score the same
  whichever way round two motifs are compared; the R package's choice
  depends on the order in which offsets are tried, so its score matrices
  can be asymmetric. This also affects which motifs `view_motifs_lite()`
  shows as reverse complements.
- **RNA consensus strings use U.** `compare_motifs()` writes T in the
  consensus strings of aligned RNA motifs, unlike every other RNA consensus
  in the package. The port writes U.
- **The parametric null gives P-values of exactly 0.** The parametric
  null describes random columns with a coarse grid of 56 columns, which
  some real columns never match perfectly, so a near-perfect alignment
  such as a motif with itself can score beyond anything the null reaches
  and get a P-value of 0. The port gives such scores the probability of
  the null's highest score, a small but honest upper bound.
- **Aligned consensus strings follow their own rule.** The consensus
  strings in `compare_motifs()` results take the most probable letters of
  each column until they reach 75%, which can give a different letter from
  the motif's own consensus. The port uses the motif consensus rule.
- **WPCC weights are not normalised.** `compare_columns(method = "WPCC")`
  weights each letter by the sum of the two backgrounds, which adds up to
  2, and uses those weights directly as if they summed to 1, so the
  weighted means are doubled and WPCC with a uniform background is not
  PCC. The port divides the weights by their sum.
- **KL and ALLR give NaN for zero probabilities.** `0 * log(0)` is NaN in
  R, so a letter with probability zero in both columns makes the whole
  score NaN. The port takes those terms as zero; a zero facing a nonzero
  probability still gives an infinite score, as it should.

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
- **`calc_complexity()` scores each string against its own letters.**
  Without `alph`, each string's alphabet is the set of letters it contains,
  so ATATATATAT counts as a two-letter sequence and scores twice as complex
  as it does in `sequence_complexity()`, which uses ACGT for DNA, and a
  homopolymer gives NaN. The port uses ACGT or ACGU for DNA and RNA in both
  functions, and otherwise the letters of all the strings together.
- **`get_bkg(to.meme = ...)` leaves trailing spaces.** Each line of the
  background file ends in a space. The port leaves it out.
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

- **`mask_ranges()` fails for common inputs.** The replacement strings are
  built with `sapply()`, which collapses them into a character vector when
  every sequence has exactly one range, and `replaceAt()` then fails; a
  sequence with no ranges at all also raises an error. The port masks the
  ranges given and leaves other sequences unchanged.

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
- **`motif_finder()` P-values ignore the search.** Each motif is the best
  of the thousands of words tried at each width, but its P-value is that
  of a single Fisher's exact test, and the Benjamini-Hochberg correction
  only counts the motifs that were kept. Random sequences compared with
  random sequences therefore give "motifs" with q-values around 1e-4. The
  port multiplies each P-value by the number of words tried in the search
  that found it and counts every word tried at every width in the
  q-values, which in tests on random sequences leaves nothing below 0.05;
  with that, the default q-value cutoff is 0.05 rather than 0.001. The
  motifs found are the same.
- **`motif_coocc(self.pairs = TRUE)` tests a motif against itself.** The
  2x2 table of a motif with itself has nothing in the off-diagonal cells,
  so the test always gives an infinite odds ratio and a P-value near zero,
  and those P-values also enter the multiple testing correction of the
  other pairs. A site on both strands at one position, as palindromes
  give, also counts as two hits at distance zero. The port does not test
  self pairs: `both` counts the sequences with hits at two or more
  positions and `a_only` those with a hit at one, and hits at the same
  position count once.
- **`motif_peaks()` orders motifs by locale.** Motifs are numbered (and
  ties in P-value broken) by `sort()` on their names, which follows the
  locale's collation rules, so `motif.i` can differ between machines. The
  port sorts names by code point.
- **`motif_proximity()` numbers motifs by name.** `motif.i` is a motif's
  rank among the sorted names rather than its place in the input, and ties
  in P-value are broken in that order. The port keeps the index that
  `scan_sequences()` gave each motif, or numbers motifs in order of first
  appearance.

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
- **Motif match matrices.** `scan::motif_matches` summarises the hits of
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
- **Matching on signal.** The Python package's `signal_stats()` summarises
  bigWig signal over ranges, and `match_bkg()` can match the windows it
  draws from a genome to the targets on that signal, such as chromatin
  accessibility, as well as on GC content and length.
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
  polygons exactly (they are tested against them), while trees and peak
  plots follow the R package's layouts without matching them point for
  point. Tree tips are ordered as `hclust()` orders them rather than
  ladderized as ggtree does, and the unrooted `radial`, `equal_angle` and
  `daylight` layouts are not available. Plots are not part of the Rust
  crate.
- `motif_pvalue(method = "exhaustive")`: the branch-and-bound P-value
  search and the random score search. Only the exact dynamic method is
  ported.
