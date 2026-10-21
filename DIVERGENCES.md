# Differences from universalmotif

umotif follows the behaviour of the universalmotif R package, and its test
suite checks the results against reference values exported from R (see
`scripts/make_reference.R`). Where universalmotif has a bug, umotif fixes it
rather than copying it, and the change is listed here so the R package can be
fixed to match.

## Fixed in umotif

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

- Higher-order backgrounds and multifreq matrices.
- Schneider's small-sample correction for information content.
- The non-lite versions of functions that have a lite counterpart in
  universalmotif.
