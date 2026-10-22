# Differences from universalmotif

umotif follows the behaviour of the universalmotif R package, and its test
suite checks the results against reference values exported from R (see
`scripts/make_reference.R`). Where universalmotif has a bug, umotif fixes it
rather than copying it, and the change is listed here so the R package can be
fixed to match.

## Fixed in umotif

### Motif objects

- **`trim_motifs()` invents `nsites`.** universalmotif writes `nsites = 100`
  into any trimmed motif that had no `nsites`, because it fills the slot in
  before computing the per-column IC. umotif uses 100 only for that
  computation and leaves `nsites` empty.
- **`create_motif()` fails on single-column matrices with row names.** The
  rows are sorted with `matrix[order(rownames(matrix)), ]` without
  `drop = FALSE`, so a one-column matrix collapses to a vector, loses its row
  names, and the call errors. umotif handles these matrices.
- **`motif_rc()` drops gap information.** The motif is rebuilt without its
  `gapinfo` slot. umotif mirrors the gap positions onto the reverse
  complement.
- **Trimming leaves gaps in the wrong place.** `trim_motifs()` and
  `trim_cwm()` remove columns but keep the old gap positions. umotif shifts
  them and drops any gap that no longer falls inside the motif.
- **Gap positions are only checked when there are several gaps.**
  `validObject_universalmotif()` skips the range check on `gaploc` when it has
  length 1. umotif always checks it, and also enforces the
  `mingap <= maxgap` and `maxgap > 0` rules that universalmotif only applies
  in `add_gap()`.
- **Some invalid values pass validation.** umotif also rejects matrices
  containing NaN and pseudocounts that are negative or not finite.
- **Letters outside a custom alphabet are ignored.** When `create_motif()` is
  given several strings and an explicit custom alphabet, letters not in the
  alphabet are silently dropped from the counts. umotif reports an error.
- **Any 20-row matrix is taken to be amino acids.** When guessing the
  alphabet of a labelled matrix, universalmotif treats every 20-row matrix as
  AA regardless of its row names. umotif requires the row names to be the
  twenty amino acids and otherwise makes a custom alphabet.

### Motif files

- **TRANSFAC files with a header block cannot be read.** `read_transfac()`
  treats the `VV ... //` header at the top of TRANSFAC releases as a motif
  and fails because it has no matrix. umotif skips records without a matrix.
- **Only the first word of a field is kept.** `read_transfac()` keeps only
  the first word of `HC` and `OS` (so `OS Homo sapiens` gives `Homo`), and
  `read_cisbp()` does the same for `TF Name`, `Family` and `Species`. umotif
  keeps the whole value.
- **TRANSFAC names depend on field order.** When the `ID` line comes before
  the `AC` line, `read_transfac()` uses the identifier as both the name and
  the alternate name. umotif always uses the accession as the name.
- **TRANSFAC and JASPAR rows are assumed to be in alphabetical order.**
  `read_transfac()` ignores the letters on the `P0` line and `read_jaspar()`
  ignores the order of the row labels. umotif reorders the values to match
  the labels.
- **Unusual JASPAR alphabets.** `read_jaspar()` treats any labelled matrix
  that is not RNA or protein as DNA, and any unlabelled matrix as DNA.
  umotif makes a custom alphabet from unknown labels and reads unlabelled
  20-row matrices as protein.
- **HOMER names with spaces.** `read_homer()` splits header lines on any
  whitespace. umotif splits tab-separated headers on tabs only.
- **MEME custom alphabet backgrounds.** `read_meme()` assigns background
  frequencies by position, which puts them on the wrong letters when MEME's
  letter order (letters before digits) differs from character order. umotif
  matches them by letter.
- **MEME sites at positions of 10,000 or more.** In BLOCKS sections the
  position then fills its parentheses (`(12345)`), the columns shift, and
  `read_meme()` returns the wrong field as the site. umotif finds the site
  after the closing parenthesis.
- **MEME block diagrams over three lines.** `read_meme()` joins only the
  first two lines of a combined block diagram. umotif joins all of them.
- **Numbers are parsed with correct rounding.** R's parser can be one unit
  in the last place off (`as.numeric("5.3e-032")` is not the closest double
  to 5.3e-32), so values read by the two packages can differ in the last
  bit.

## Additions

- **CIS-BP and UniPROBE writers.** universalmotif can read these formats but
  not write them. umotif writes both (UniPROBE in its minimal form), and its
  files read correctly with universalmotif's `read_cisbp()` and
  `read_uniprobe()`.

## Not reproducible

- **Random motifs.** `create_motif()` with a width draws from R's random
  number generator, which umotif does not reproduce. umotif also gives a
  random PPM a random `nsites` whenever none is supplied, whereas
  universalmotif only does so when the type was left at its default.

## Known issues shared with universalmotif

- **`ppm_to_pcm()` can produce negative counts.** The rounding error is added
  to the largest count, which goes negative when many small probabilities all
  round up (for example `create_motif("BZJX", nsites = 50)`). Both packages
  then reject the motif. Fixing this needs a different rounding scheme, which
  would change counts for ordinary motifs too, so it is left for now.

## Not supported by design

- universalmotif's own YAML format (`read_motifs()` and `write_motifs()`).
- Higher-order backgrounds and multifreq matrices.
- Schneider's small-sample correction for information content.
- The non-lite versions of functions that have a lite counterpart in
  universalmotif.
