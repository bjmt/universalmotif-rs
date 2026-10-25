#!/usr/bin/env Rscript
#
# Export reference values from the universalmotif R package for randomized
# inputs, to cross-check the crate beyond the hand-picked cases of
# make_reference.R.
#
# Usage (from the root of this repository):
#   Rscript scripts/make_reference_random.R [path/to/universalmotif]
#
# The motifs, sequences and arguments are drawn from R's random number
# generator with a fixed seed and recorded alongside the results, so the
# crate's tests read their inputs from the JSON files rather than
# regenerating them. Output goes to tests/reference/random_*.json.

args <- commandArgs(trailingOnly = TRUE)
um_path <- if (length(args)) args[1] else "~/universalmotif"
um_path <- normalizePath(um_path)

suppressMessages(pkgload::load_all(um_path, quiet = TRUE, export_all = TRUE))
library(jsonlite)

out_dir <- file.path("tests", "reference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

commit <- system2("git", c("-C", um_path, "rev-parse", "HEAD"), stdout = TRUE)
um_version <- as.character(utils::packageVersion("universalmotif", lib.loc = NULL))

opt <- function(x) if (length(x)) x else NULL
vec <- function(x) I(unname(as.numeric(x)))

motif_record <- function(m) {
  ei <- m@extrainfo
  list(
    name = m@name,
    altname = opt(m@altname),
    family = opt(m@family),
    organism = opt(m@organism),
    alphabet = m@alphabet,
    type = m@type,
    icscore = m@icscore,
    nsites = opt(m@nsites),
    pseudocount = m@pseudocount,
    bkg = vec(m@bkg[rownames(m@motif)]),
    bkgsites = opt(m@bkgsites),
    consensus = opt(m@consensus),
    strand = m@strand,
    pval = opt(m@pval),
    qval = opt(m@qval),
    eval = opt(m@eval),
    extrainfo = lapply(seq_along(ei), function(i) {
      nm <- names(ei)[i]
      I(c(if (is.null(nm) || is.na(nm)) "" else nm, unname(ei[i])))
    }),
    gap = list(
      isgapped = m@gapinfo@isgapped,
      gaploc = vec(m@gapinfo@gaploc),
      mingap = vec(m@gapinfo@mingap),
      maxgap = vec(m@gapinfo@maxgap)
    ),
    matrix = lapply(seq_len(ncol(m@motif)), function(j) vec(m@motif[, j]))
  )
}

attempt <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = suppressMessages(suppressWarnings(expr))),
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
}

as_motif_result <- function(r) {
  if (r$ok) list(motif = motif_record(r$value)) else list(error = r$error)
}

num_or_error <- function(expr) {
  r <- attempt(expr)
  if (r$ok) list(value = vec(r$value)) else list(error = r$error)
}

df_cols <- function(d) {
  d <- as.data.frame(d)
  lapply(d, function(x) if (is.numeric(x)) vec(x) else I(as.character(x)))
}

rows <- function(m) lapply(seq_len(nrow(m)), function(i) vec(m[i, ]))

write_ref <- function(x, file) {
  x <- c(list(meta = list(universalmotif_commit = commit,
                          universalmotif_version = um_version,
                          generated = format(Sys.time(), "%Y-%m-%d"),
                          seed = 20261025)),
         x)
  writeLines(toJSON(x, auto_unbox = TRUE, null = "null", na = "string",
                    digits = I(17), pretty = FALSE),
             file.path(out_dir, file))
  cat("wrote", file.path(out_dir, file),
      round(file.size(file.path(out_dir, file)) / 1024), "KB\n")
}

set.seed(20261025)

# ---------------------------------------------------------------------------
# Random motifs
#
# Columns are drawn from gamma variates of varying shape, so some are nearly
# uniform and some nearly pure, and a fifth of them have zeros. Widths cover
# one and two columns as well as 15 to 25.

alph_rows <- list(DNA = DNA_BASES, RNA = RNA_BASES, AA = AA_STANDARD2,
                  QRSTU = c("Q", "R", "S", "T", "U"))

rand_col <- function(n) {
  x <- rgamma(n, shape = sample(c(0.15, 0.4, 1, 4), 1))
  if (runif(1) < 0.2) x[sample.int(n, max(1, n %/% 4))] <- 0
  if (sum(x) == 0) x[sample.int(n, 1)] <- 1
  x / sum(x)
}
rand_bkg <- function(n) { b <- runif(n, 0.3, 1); b / sum(b) }

n_motifs <- 44
widths <- sample(c(1, 1, 2, 2, sample(3:12, 20, TRUE), sample(15:25, 20, TRUE)))
alphs <- sample(c(rep("DNA", 22), rep("RNA", 8), rep("AA", 7), rep("QRSTU", 7)))
types <- sample(rep(c("PCM", "PPM", "PWM", "ICM"), 11))
# Wide amino acid motifs make large records; keep them to a few columns.
widths[alphs == "AA"] <- pmin(widths[alphs == "AA"], sample(c(1, 2, 4, 6, 9, 15, 18), 7))

make_motif <- function(i) {
  alph <- alphs[i]
  rws <- alph_rows[[alph]]
  n <- length(rws)
  w <- widths[i]
  probs <- matrix(vapply(seq_len(w), function(j) rand_col(n), numeric(n)),
                  nrow = n, dimnames = list(rws, NULL))
  nsites <- sample(c(5, 7, 13, 21, 37, 51, 99, 101, 250), 1)
  pseudo <- sample(c(0, 0.5, 1, 1, 2, 3), 1)
  bkg <- if (runif(1) < 0.5) NULL else setNames(rand_bkg(n), rws)
  name <- sprintf("m%02d", i)
  # create_motif(bkg = NULL) aborts R for custom alphabets, so the argument
  # is only passed when there is a background.
  extra <- list(alphabet = alph, pseudocount = pseudo, name = name)
  if (!is.null(bkg)) extra$bkg <- bkg
  if (types[i] == "PCM") {
    counts <- apply(probs, 2, function(p) rmultinom(1, nsites, p))
    counts <- matrix(counts, nrow = n, dimnames = list(rws, NULL))
    m <- do.call(create_motif, c(list(input = counts, type = "PCM"), extra))
  } else {
    if (runif(1) >= 0.25) extra$nsites <- nsites
    m <- do.call(create_motif, c(list(input = probs, type = "PPM"), extra))
    if (types[i] != "PPM") m <- suppressMessages(convert_type(m, types[i]))
  }
  m@multifreq <- list()
  m
}
motifs <- lapply(seq_len(n_motifs), make_motif)
is_nt <- vapply(motifs, function(m) m@alphabet %in% c("DNA", "RNA"), logical(1))

letters_of <- function(m) rownames(m@motif)
random_words <- function(m, n) {
  vapply(seq_len(n), function(i) paste(sample(letters_of(m), ncol(m@motif), TRUE),
                                        collapse = ""), character(1))
}
# A word drawn from the motif's own probabilities, so that it tends to score
# high.
motif_word <- function(m) {
  p <- suppressMessages(convert_type(m, "PPM"))@motif
  paste(apply(p, 2, function(col) sample(rownames(p), 1, prob = col + 1e-9)),
        collapse = "")
}

# convert_type() between all four types, with a new pseudocount for some.
convert_cases <- list()
for (i in seq_along(motifs)) {
  m <- motifs[[i]]
  for (t in c("PCM", "PPM", "PWM", "ICM")) {
    variants <- list(list())
    if (t != "PCM" && i %% 3 == 0) variants <- c(variants, list(list(pseudocount = 0.25)))
    if (t == "ICM" && i %% 4 == 0) variants <- c(variants, list(list(relative_entropy = TRUE)))
    for (v in variants) {
      r <- attempt(do.call(convert_type, c(list(motifs = m, type = t), v)))
      convert_cases[[length(convert_cases) + 1]] <- c(
        list(input = i - 1L, args = c(list(type = t), v)), as_motif_result(r))
    }
  }
}

# Reverse complements, trimming and rounding.
ops_cases <- list()
add_op <- function(op, i, args, expr) {
  r <- attempt(expr)
  ops_cases[[length(ops_cases) + 1]] <<- c(list(op = op, input = i - 1L, args = args),
                                          as_motif_result(r))
}
for (i in seq_along(motifs)) {
  m <- motifs[[i]]
  if (is_nt[i]) add_op("motif_rc", i, list(), motif_rc(m))
  ic <- round(runif(1, 0.05, 1.2), 3)
  tf <- sample(c("both", "left", "right"), 1)
  add_op("trim_motifs", i, list(min_ic = ic, trim_from = tf),
         trim_motifs(m, min.ic = ic, trim.from = tf))
  tol <- sample(c(0.02, 0.05, 0.1, 0.2), 1)
  add_op("round_motif", i, list(pct_tolerance = tol), round_motif(m, tol))
}

avg_cases <- lapply(c("a.mean", "g.mean", "median", "fzt"), function(a)
  list(average = a, value = vec(average_ic(motifs, a))))

# Scoring helpers.
scoring <- list()
add_scoring <- function(fn, i, args, out) {
  scoring[[length(scoring) + 1]] <<- c(list(fn = fn, motif = i - 1L, args = args), out)
}
for (i in seq_along(motifs)) {
  m <- motifs[[i]]
  add_scoring("motif_range", i, list(), num_or_error(motif_range(m)))
  th <- sort(round(c(0, runif(3), 1), 4))
  for (tt in c("total", "fromzero"))
    add_scoring("motif_score", i, list(threshold = vec(th), threshold_type = tt),
                num_or_error(motif_score(m, th, threshold.type = tt)))
  words <- c(random_words(m, 6), vapply(1:3, function(k) motif_word(m), character(1)))
  add_scoring("score_match", i, list(matches = I(words)),
              num_or_error(score_match(m, words)))
  for (az in c(TRUE, FALSE))
    add_scoring("prob_match", i, list(matches = I(words), allow_zero = az),
                num_or_error(prob_match(m, words, allow.zero = az)))
  if (nrow(m@motif)^ncol(m@motif) <= 4096) {
    all <- suppressMessages(get_scores(m))
    add_scoring("get_scores", i, list(), list(value = vec(all)))
    # Scores taken from the grid itself, since the R package truncates
    # scores between grid points where the port rounds them up.
    for (s in unique(all[c(1, max(1, length(all) %/% 50), max(1, length(all) %/% 5))])) {
      r <- attempt(get_matches(m, s))
      add_scoring("get_matches", i, list(score = s),
                  if (r$ok) list(value = I(r$value)) else list(error = r$error))
    }
  }
}

# motif_pvalue(method = "dynamic") in both directions, with the motif's own
# background and with another one.
pv_cases <- list()
pvals_in <- c(0.5, 0.1, 0.01, 1e-3, 1e-4, 1e-5, 1e-7, 1e-10)
for (i in seq_along(motifs)) {
  m <- motifs[[i]]
  rng <- suppressMessages(motif_range(m))
  scores <- c(round(seq(rng[1], rng[2], length.out = 9), 3),
              round(runif(3, rng[1], rng[2]), 3), rng[2] + 0.5)
  bkgs <- list(NULL)
  if (i %% 2 == 0) bkgs <- c(bkgs, list(setNames(rand_bkg(nrow(m@motif)), letters_of(m))))
  for (b in bkgs) {
    r1 <- attempt(if (is.null(b)) motif_pvalue(m, score = scores)
                  else motif_pvalue(m, score = scores, bkg.probs = b))
    r2 <- attempt(if (is.null(b)) motif_pvalue(m, pvalue = pvals_in)
                  else motif_pvalue(m, pvalue = pvals_in, bkg.probs = b))
    pv_cases[[length(pv_cases) + 1]] <- list(
      motif = i - 1L, bkg = if (is.null(b)) NULL else vec(b),
      scores = vec(scores), pvalues = if (r1$ok) vec(r1$value) else NULL,
      pvalue_in = vec(pvals_in), score_out = if (r2$ok) vec(r2$value) else NULL,
      error = c(if (!r1$ok) r1$error, if (!r2$ok) r2$error))
  }
}

write_ref(list(motifs = lapply(motifs, motif_record), convert = convert_cases,
               ops = ops_cases, average_ic = avg_cases, scoring = scoring,
               pvalue = pv_cases), "random_motifs.json")

# ---------------------------------------------------------------------------
# Scanning and enrichment
#
# DNA and RNA motifs go through scan_sequences_lite() and the others through
# the full scan_sequences() on the forward strand, as the port's alphabet
# tests do. Sequences vary in composition and carry sites drawn from the
# motifs.

rand_seq <- function(letters, len, prob = NULL)
  paste(sample(letters, len, TRUE, prob = prob), collapse = "")
plant_words <- function(seqs, ms, n) {
  for (k in seq_len(n)) {
    j <- sample.int(length(seqs), 1)
    w <- motif_word(ms[[sample.int(length(ms), 1)]])
    if (nchar(w) >= nchar(seqs[j])) next
    p <- sample.int(nchar(seqs[j]) - nchar(w) + 1, 1)
    substr(seqs[j], p, p + nchar(w) - 1) <- w
  }
  seqs
}
seq_set <- function(letters, ms, n, lens, nplant) {
  s <- vapply(seq_len(n), function(k)
    rand_seq(letters, sample(lens, 1), prob = rand_bkg(length(letters))), character(1))
  s <- plant_words(s, ms, nplant)
  names(s) <- sprintf("s%02d", seq_len(n))
  s
}

# scan_sequences_lite() scores a PCM that has zero counts and no pseudocount
# with infinite entries, because normalize() leaves PCMs unchanged, and then
# reports hits at nearly every position. Such motifs are scanned here as
# PPMs, which the R package handles correctly and the port treats the same
# way as the PCMs themselves.
scan_motifs <- lapply(motifs, function(m) {
  if (m@type == "PCM" && m@pseudocount == 0 && any(m@motif == 0))
    suppressMessages(convert_type(m, "PPM"))
  else m
})

dna_i <- which(vapply(motifs, function(m) m@alphabet == "DNA", logical(1)))
rna_i <- which(vapply(motifs, function(m) m@alphabet == "RNA", logical(1)))
aa_i <- which(vapply(motifs, function(m) m@alphabet == "AA", logical(1)))
cu_i <- which(vapply(motifs, function(m) m@alphabet == "QRSTU", logical(1)))

dna_seqs <- seq_set(DNA_BASES, motifs[dna_i], 16, 30:300, 40)
rna_seqs <- seq_set(RNA_BASES, motifs[rna_i], 8, 30:200, 16)
aa_seqs <- seq_set(AA_STANDARD2, motifs[aa_i], 8, 30:200, 16)
cu_seqs <- seq_set(alph_rows$QRSTU, motifs[cu_i], 8, 30:200, 16)

scan_cases <- list()
add_scan <- function(set, idx, seqs, args) {
  r <- attempt(do.call(scan_sequences_lite,
                       c(list(motifs = scan_motifs[idx], sequences = seqs,
                              return.granges = FALSE), args)))
  case <- list(set = set, motifs = I(idx - 1L), args = args)
  if (r$ok) {
    d <- r$value
    case$hits <- list(
      motif_i = I(d$motif.i), sequence_i = I(d$sequence.i), start = I(d$start),
      end = I(d$end), strand = I(d$strand), score = vec(d$score),
      score_pct = vec(d$score.pct), match = I(d$match), pvalue = vec(d$pvalue))
  } else case$error <- r$error
  scan_cases[[length(scan_cases) + 1]] <<- case
}
dna_set <- Biostrings::DNAStringSet(dna_seqs)
rna_set <- Biostrings::RNAStringSet(rna_seqs)
add_scan("dna", dna_i, dna_set, list(pvalue = 1e-2))
add_scan("dna", dna_i, dna_set, list(pvalue = 1e-3, RC = FALSE))
add_scan("dna", dna_i, dna_set, list(pvalue = 1e-4))
add_scan("dna", dna_i, dna_set, list(pvalue = 1e-3, no.overlaps = TRUE))
add_scan("dna", dna_i, dna_set, list(pvalue = 1e-2, no.overlaps = TRUE,
                                     no.overlaps.by = "score",
                                     no.overlaps.by.strand = TRUE))
add_scan("rna", rna_i, rna_set, list(pvalue = 1e-2))
add_scan("rna", rna_i, rna_set, list(pvalue = 1e-3, RC = FALSE))

full_scan <- function(idx, seqs, threshold) {
  r <- suppressMessages(suppressWarnings(scan_sequences(
    scan_motifs[idx], seqs, threshold = threshold, threshold.type = "pvalue", RC = FALSE)))
  r <- as.data.frame(r)
  r <- r[order(r$motif.i, r$sequence.i, r$start), ]
  df_cols(r[, c("motif.i", "sequence.i", "start", "stop", "score", "match", "pvalue")])
}
other_scans <- list(
  list(set = "aa", motifs = I(aa_i - 1L), threshold = 1e-2,
       table = full_scan(aa_i, Biostrings::AAStringSet(aa_seqs), 1e-2)),
  list(set = "aa", motifs = I(aa_i - 1L), threshold = 1e-4,
       table = full_scan(aa_i, Biostrings::AAStringSet(aa_seqs), 1e-4)),
  list(set = "custom", motifs = I(cu_i - 1L), threshold = 1e-2,
       table = full_scan(cu_i, Biostrings::BStringSet(cu_seqs), 1e-2)))

# Enrichment against a fixed background, so nothing is shuffled.
en_tgt <- seq_set(DNA_BASES, motifs[dna_i[1:8]], 40, 80:200, 60)
en_bkg <- seq_set(DNA_BASES, motifs[dna_i], 50, 80:200, 0)
enrich_cases <- list()
for (args in list(list(qvalue = 1, pvalue = 1e-3),
                  list(qvalue = 1, pvalue = 1e-2, test = "sites"),
                  list(qvalue = 1, pvalue = 1e-3, RC = FALSE, pseudocount = 0),
                  list(qvalue = 1, pvalue = 1e-4, test = "sites", pseudocount = 2))) {
  r <- suppressMessages(do.call(enrich_motifs_lite, c(list(
    motifs = scan_motifs[dna_i], sequences = Biostrings::DNAStringSet(en_tgt),
    bkg.sequences = Biostrings::DNAStringSet(en_bkg)), args)))
  enrich_cases[[length(enrich_cases) + 1]] <- list(args = args, table = df_cols(r))
}

write_ref(list(dna_sequences = I(unname(dna_seqs)), rna_sequences = I(unname(rna_seqs)),
               aa_sequences = I(unname(aa_seqs)), custom_sequences = I(unname(cu_seqs)),
               scan = scan_cases, full_scan = other_scans,
               enrich_targets = I(unname(en_tgt)), enrich_background = I(unname(en_bkg)),
               enrich_motifs = I(dna_i - 1L), enrich = enrich_cases),
          "random_scan.json")

# ---------------------------------------------------------------------------
# Comparison and merging
#
# compare_motifs_lite() is the comparison the port follows. For amino acid
# and custom motifs, the full compare_motifs() with summed Pearson
# correlations and no IC filtering picks the same alignments.

cmp_cases <- list()
add_cmp <- function(set, idx, args) {
  r <- attempt(do.call(compare_motifs_lite, c(list(motifs = motifs[idx]), args)))
  case <- list(set = set, motifs = I(idx - 1L), args = args)
  if (!r$ok) {
    case$error <- r$error
  } else if (is.matrix(r$value)) {
    case$matrix <- rows(r$value)
  } else {
    d <- r$value
    case$long <- list(subject_i = I(d$subject.i), target_i = I(d$target.i),
                      offset = I(d$offset), strand = I(d$strand),
                      overlap = I(d$overlap), score = vec(d$score),
                      pvalue = vec(d$pvalue), qvalue = vec(d$qvalue))
  }
  cmp_cases[[length(cmp_cases) + 1]] <<- case
}
for (out in c("score", "pvalue", "qvalue"))
  for (rc in c(TRUE, FALSE))
    add_cmp("dna", dna_i, list(matrix.out = out, RC = rc,
                               min.overlap = sample(c(1L, 3L, 6L), 1)))
for (mo in c(1L, 4L, 12L))
  add_cmp("dna", dna_i, list(matrix.out = "pvalue", min.overlap = mo))
add_cmp("dna", dna_i, list(matrix.out = "pvalue", null = "parametric"))
add_cmp("dna", dna_i, list(matrix.out = "qvalue", null = "parametric", RC = FALSE,
                           bkg = setNames(rand_bkg(4), DNA_BASES)))
add_cmp("dna", dna_i, list(compare.to = sample.int(length(dna_i), 3), qvalue = 1))
add_cmp("dna", dna_i, list(compare.to = sample.int(length(dna_i), 2), qvalue = 1,
                           RC = FALSE, min.overlap = 3L, null = "parametric"))
add_cmp("rna", rna_i, list(matrix.out = "pvalue", min.overlap = 2L))
add_cmp("rna", rna_i, list(compare.to = 1L, qvalue = 1, RC = FALSE))

other_cmp <- function(set, idx, min.overlap) {
  s <- suppressMessages(compare_motifs(motifs[idx], method = "PCC", score.strat = "sum",
                                       min.overlap = min.overlap, min.mean.ic = 0,
                                       min.position.ic = 0, normalise.scores = FALSE,
                                       tryRC = FALSE))
  list(set = set, motifs = I(idx - 1L), min_overlap = min.overlap,
       sums = rows(as.matrix(s)))
}
other_compare <- list(other_cmp("aa", aa_i, 1L), other_cmp("aa", aa_i, 3L),
                      other_cmp("custom", cu_i, 2L), other_cmp("custom", cu_i, 5L))

merge_cases <- list()
for (k in 1:12) {
  pool <- if (k <= 9) dna_i else rna_i
  idx <- sample(pool, sample(2:4, 1))
  args <- list(min.overlap = sample(c(1L, 3L, 5L), 1), RC = runif(1) < 0.7,
               weighted = runif(1) < 0.4)
  r <- attempt(do.call(merge_motifs_lite, c(list(motifs = motifs[idx]), args)))
  merge_cases[[length(merge_cases) + 1]] <- c(list(motifs = I(idx - 1L), args = args),
                                              as_motif_result(r))
}

# Column metrics on random columns of 4, 5 and 20 letters.
metric_cases <- list()
for (k in 1:30) {
  n <- c(4, 5, 20)[k %% 3 + 1]
  x <- rand_col(n); y <- rand_col(n)
  b1 <- rand_bkg(n); b2 <- if (k %% 2) b1 else rand_bkg(n)
  ns <- sample(c(5, 13, 37, 100), 2, TRUE)
  for (method in COMPARE_METRICS) {
    r <- attempt(compare_columns(x, y, method, b1, b2, ns[1], ns[2]))
    metric_cases[[length(metric_cases) + 1]] <- c(
      list(x = vec(x), y = vec(y), method = method, bkg1 = vec(b1), bkg2 = vec(b2),
           nsites1 = ns[1], nsites2 = ns[2]),
      if (r$ok) list(value = r$value) else list(error = r$error))
  }
}

write_ref(list(compare = cmp_cases, full_compare = other_compare, merge = merge_cases,
               compare_columns = metric_cases), "random_compare.json")

# ---------------------------------------------------------------------------
# Backgrounds and complexity of random sequences, with low complexity
# stretches and runs mixed in.

cx_seqs <- vapply(1:6, function(k) {
  s <- rand_seq(DNA_BASES, sample(40:160, 1), prob = rand_bkg(4))
  for (r in 1:2) {
    run <- switch(sample(3, 1), strrep(sample(DNA_BASES, 1), sample(5:20, 1)),
                  strrep(paste(sample(DNA_BASES, 2), collapse = ""), sample(3:10, 1)),
                  strrep(paste(sample(DNA_BASES, 3), collapse = ""), sample(2:6, 1)))
    p <- sample.int(max(1, nchar(s) - nchar(run)), 1)
    substr(s, p, p + nchar(run) - 1) <- run
  }
  s
}, character(1))
names(cx_seqs) <- sprintf("c%d", 1:6)
cx_set <- Biostrings::DNAStringSet(cx_seqs)

bkg_cases <- list()
add_bkg <- function(args) {
  r <- attempt(do.call(get_bkg, c(list(sequences = cx_set), args)))
  bkg_cases[[length(bkg_cases) + 1]] <<- c(list(args = args),
    if (r$ok) list(table = df_cols(r$value)) else list(error = r$error))
}
add_bkg(list(k = 1:3))
add_bkg(list(k = 1:4, RC = TRUE))
add_bkg(list(k = 2L, merge.res = FALSE, pseudocount = 1L))
add_bkg(list(k = 1:2, window = TRUE, window.size = 17, window.overlap = 5))
add_bkg(list(k = 1:2, window = TRUE, merge.res = FALSE, window.size = 17,
             window.overlap = 5))
add_bkg(list(k = 3L, RC = TRUE, merge.res = FALSE))

seqcx_cases <- list()
for (args in list(list(),
                  list(window.size = 13, window.overlap = 4),
                  list(window.size = 0.3, window.overlap = 0.1, method = "WoottonFederhenFast"),
                  list(window.size = 25, window.overlap = 0, method = "Trifonov"),
                  list(window.size = 30, window.overlap = 12, method = "TrifonovFast",
                       trifonov.max.word.size = 5),
                  list(window.size = 24, window.overlap = 8, method = "DUST"))) {
  r <- do.call(sequence_complexity, c(list(seqs = cx_set), args))
  seqcx_cases[[length(seqcx_cases) + 1]] <- list(args = args, table = df_cols(r))
}

complexity_cases <- list()
for (s in c(cx_seqs, vapply(1:6, function(k) substr(cx_seqs[k], 1, sample(4:30, 1)),
                            character(1)))) {
  for (meth in c("WoottonFederhen", "WoottonFederhenFast", "Trifonov", "TrifonovFast",
                 "DUST")) {
    mw <- if (startsWith(meth, "Trifonov")) sample(2:7, 1) else 7L
    r <- attempt(calc_complexity(s, meth, "ACGT", mw))
    complexity_cases[[length(complexity_cases) + 1]] <- c(
      list(string = unname(s), method = meth, max_word = mw),
      if (r$ok) list(value = r$value) else list(error = r$error))
  }
}

write_ref(list(sequences = I(unname(cx_seqs)), names = I(names(cx_seqs)),
               bkg = bkg_cases, sequence_complexity = seqcx_cases,
               complexity = complexity_cases), "random_sequences.json")
