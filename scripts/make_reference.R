#!/usr/bin/env Rscript
#
# Export reference values from the universalmotif R package for this crate's
# test suite.
#
# Usage (from the root of this repository):
#   Rscript scripts/make_reference.R [path/to/universalmotif]
#
# universalmotif is loaded from source with pkgload so the values match the
# checked-out commit rather than whatever version happens to be installed.
# Output goes to tests/reference/*.json.

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

# Run an expression, returning list(ok = TRUE, value) or list(ok = FALSE,
# error) so failing cases are recorded rather than aborting the script.
attempt <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = suppressMessages(suppressWarnings(expr))),
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
}

as_motif_result <- function(r) {
  if (r$ok) list(motif = motif_record(r$value)) else list(error = r$error)
}

write_ref <- function(x, file) {
  x <- c(list(meta = list(universalmotif_commit = commit,
                          universalmotif_version = um_version,
                          generated = format(Sys.time(), "%Y-%m-%d"))),
         x)
  writeLines(toJSON(x, auto_unbox = TRUE, null = "null", na = "string",
                    digits = I(17), pretty = FALSE),
             file.path(out_dir, file))
  cat("wrote", file.path(out_dir, file), "\n")
}

# ---------------------------------------------------------------------------
# Column functions

ppms <- list(
  c(0.25, 0.25, 0.25, 0.25), c(1, 0, 0, 0), c(0.7, 0.1, 0.1, 0.1),
  c(0.5, 0.5, 0, 0), c(0.4, 0.3, 0.2, 0.1), c(0.33, 0.33, 0.33, 0.01),
  c(0.001, 0.997, 0.001, 0.001), c(0.6, 0, 0.4, 0), c(0.1, 0.45, 0.45, 0),
  c(0.3, 0.3, 0.2, 0.2), c(0.26, 0.26, 0.26, 0.22), c(0.05, 0.05, 0.1, 0.8),
  c(0.2, 0.3, 0.3, 0.2), c(0.3, 0, 0.3, 0.4)
)
pcms <- list(c(10, 0, 0, 0), c(3, 3, 2, 2), c(1, 2, 3, 4), c(50, 25, 25, 0),
             c(7, 1, 1, 1), c(0, 0, 5, 5))
pwms <- c(lapply(ppms, ppm_to_pwmC, bkg = numeric(), pseudocount = 1),
          list(c(1.2, -0.5, -3, 0.1), c(2, -Inf, -Inf, -Inf), c(-1, -1, -1, -1)))
icms <- c(lapply(ppms, ppm_to_icmC, bkg = numeric()), list(c(0.5, 0.2, 0.1, 0)))
cwms <- list(c(0.5, -0.2, 0.1, 0), c(0, 0, 0, 0), c(-0.3, -0.3, 0.6, 0.05))
bkgs <- list(numeric(), c(0.25, 0.25, 0.25, 0.25), c(0.3, 0.2, 0.2, 0.3),
             c(0.1, 0.4, 0.4, 0.1), c(0.5, 0.25, 0.15, 0.1))
pseudos <- c(0, 0.5, 1, 4)
nsites <- c(0, 1, 10, 37, 100, NaN)

cases <- list()
add <- function(fn, args, out) {
  cases[[length(cases) + 1]] <<- list(fn = fn, args = args, out = out)
}

for (p in pcms) for (ps in pseudos)
  add("pcm_to_ppm", list(pos = vec(p), pseudocount = ps), vec(pcm_to_ppmC(p, ps)))
for (p in ppms) for (n in nsites)
  add("ppm_to_pcm", list(pos = vec(p), nsites = n), vec(ppm_to_pcmC(p, n)))
for (p in ppms) for (b in bkgs) for (ps in pseudos) for (n in c(0, 37, 100))
  add("ppm_to_pwm", list(pos = vec(p), bkg = vec(b), pseudocount = ps, nsites = n),
      vec(ppm_to_pwmC(p, b, ps, n)))
for (p in pwms) for (b in bkgs)
  add("pwm_to_ppm", list(pos = vec(p), bkg = vec(b)), vec(pwm_to_ppmC(p, b)))
for (p in ppms) for (b in bkgs) for (re in c(FALSE, TRUE))
  add("ppm_to_icm", list(pos = vec(p), bkg = vec(b), relative_entropy = re),
      vec(ppm_to_icmC(p, b, re)))
for (p in icms)
  add("icm_to_ppm", list(pos = vec(p)), vec(icm_to_ppmC(p)))

by_type <- list(PCM = pcms, PPM = ppms, PWM = pwms, ICM = icms, CWM = cwms)
for (type in names(by_type)) for (p in by_type[[type]]) {
  for (b in bkgs[c(1, 3, 5)]) for (ps in c(0, 1)) for (n in c(0, 37)) for (re in c(FALSE, TRUE))
    add("position_icscore",
        list(pos = vec(p), bkg = vec(b), type = type, pseudocount = ps,
             nsites = n, relative_entropy = re),
        position_icscoreC(p, b, type, ps, n, re))
  for (alph in c("DNA", "RNA")) for (ps in c(0, 1))
    add("get_consensus",
        list(pos = vec(p), alphabet = alph, type = type, pseudocount = ps),
        get_consensusC(p, alph, type, ps))
}

aa_cols <- c(
  lapply(c(AA_STANDARD2, "B", "Z", "J", "X"), consensus_to_ppmAAC),
  list(rep(0.05, 20),
       c(0.3, 0.3, rep(0.4 / 18, 18)),
       c(0.15, rep(0.85 / 19, 19)),
       c(rep(0.01, 19), 0.81))
)
for (p in aa_cols) for (type in c("PPM", "PWM", "ICM")) {
  pos <- switch(type, PPM = p, PWM = ppm_to_pwmC(p, numeric()), ICM = ppm_to_icmC(p, numeric()))
  add("get_consensus_aa", list(pos = vec(pos), type = type, pseudocount = 0),
      get_consensusAAC(pos, type, 0))
}
for (p in list(c(5, rep(0, 19)), c(rep(1, 20)), c(4, 4, rep(1, 18))))
  add("get_consensus_aa", list(pos = vec(p), type = "PCM", pseudocount = 1),
      get_consensusAAC(p, "PCM", 1))

for (l in c(DNA_ALPHABET, "U"))
  add("consensus_to_ppm", list(letter = l), vec(consensus_to_ppmC(l)))
for (l in AA_ALPHABET)
  add("consensus_to_ppm_aa", list(letter = l), vec(consensus_to_ppmAAC(l)))

write_ref(list(cases = cases), "columns.json")

# ---------------------------------------------------------------------------
# Building motif objects (universalmotif_cpp)

mat <- function(cols, rows = c("A", "C", "G", "T")) {
  matrix(unlist(cols), nrow = length(rows), dimnames = list(rows, NULL))
}

build_inputs <- list(
  list(matrix = list(c(10, 0, 0, 0), c(3, 3, 2, 2))),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1), c(0.25, 0.25, 0.25, 0.25))),
  list(matrix = list(c(0.5, 0.3, 0.1, 0.05), c(0.2, 0.2, 0.2, 0.2)), type = "PPM"),
  list(matrix = list(c(1.5, 0.2, 0.1, 0), c(0.3, 0.3, 0.3, 0.3))),
  list(matrix = list(c(1.5, -0.2, -1, 0), c(-2, 0.3, 0.3, 0.3))),
  list(matrix = list(c(6, 2, 2, 0), c(1, 2, 3, 4)), nsites = 20),
  list(matrix = list(c(6, 2, 2, 0), c(1, 2, 3, 6))),
  list(matrix = list(c(6, 2, 2, 0), c(1, 2, 3, 4)), bkg = c(0.3, 0.2, 0.2, 0.3)),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), name = ""),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), alphabet = "RNA"),
  list(matrix = list(c(0.2, 0.5, 0.3), c(0.1, 0.1, 0.8)), alphabet = "XYZ"),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), strand = "-+"),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), strand = "-", pval = 1e-5,
       qval = 0.01, eval = 3, bkgsites = 500, altname = "alt", family = "fam",
       organism = "org", extrainfo = c(source = "test", "unnamed")),
  list(matrix = list(c(0.5, -0.2, 0.1, 0), c(0.1, 0.1, -0.6, 0)), type = "CWM"),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), nsites = 42),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), icscore = 3.5),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1), c(1, 0, 0, 0), c(0, 1, 0, 0),
                     c(0, 0, 0.5, 0.5)),
       isgapped = TRUE, gaploc = 2, mingap = 1, maxgap = 3),
  list(matrix = list(c(0.7, 0.1, 0.1, 0.1)), pseudocount = 0),
  list(matrix = list(c(2, 0, 0, 0), c(1, 1, 0, 0)), pseudocount = 0.5)
)
aa_mat <- lapply(c("A", "W", "B"), consensus_to_ppmAAC)
build_inputs[[length(build_inputs) + 1]] <- list(matrix = aa_mat, alphabet = "AA")

build_cases <- lapply(build_inputs, function(inp) {
  alph <- if (is.null(inp$alphabet)) "DNA" else inp$alphabet
  rows <- switch(alph, DNA = DNA_BASES, RNA = RNA_BASES, AA = AA_STANDARD2,
                 strsplit(alph, "")[[1]])
  m <- matrix(unlist(inp$matrix), nrow = length(rows), dimnames = list(rows, NULL))
  cargs <- list(motif = m, alphabet = alph)
  for (nm in c("name", "type", "nsites", "pseudocount", "bkg", "icscore",
               "strand", "pval", "qval", "eval", "bkgsites", "altname",
               "family", "organism", "extrainfo", "isgapped", "gaploc",
               "mingap", "maxgap"))
    if (!is.null(inp[[nm]])) cargs[[nm]] <- inp[[nm]]
  r <- attempt({
    x <- do.call(universalmotif_cpp, cargs)
    validObject_universalmotif(x)
    x
  })
  inp$matrix <- lapply(inp$matrix, vec)
  if (!is.null(inp[["bkg"]])) inp[["bkg"]] <- vec(inp[["bkg"]])
  if (!is.null(inp$extrainfo)) {
    ei <- inp$extrainfo
    inp$extrainfo <- lapply(seq_along(ei), function(i) I(c(names(ei)[i], unname(ei[i]))))
  }
  c(list(input = inp), as_motif_result(r))
})

write_ref(list(cases = build_cases), "build.json")

# ---------------------------------------------------------------------------
# create_motif()

create_cases <- list()
add_create <- function(input, args = list()) {
  call_args <- c(list(input = input$value), args)
  r <- attempt(do.call(create_motif, call_args))
  create_cases[[length(create_cases) + 1]] <<- c(
    list(input = input$record, args = args),
    as_motif_result(r)
  )
}
matrix_input <- function(cols, rows) {
  m <- matrix(unlist(cols), nrow = length(rows), dimnames = list(rows, NULL))
  list(value = m,
       record = list(kind = "matrix", matrix = lapply(cols, vec), rows = I(rows)))
}
strings_input <- function(x) {
  list(value = x, record = list(kind = "strings", strings = I(x)))
}

pcm_in <- matrix_input(list(c(5, 0, 0, 5), c(2, 2, 3, 3), c(0, 10, 0, 0)), DNA_BASES)
for (args in list(list(), list(type = "PCM"), list(type = "PWM"), list(type = "ICM"),
                  list(nsites = 20), list(pseudocount = 1),
                  list(type = "PWM", pseudocount = 1), list(name = "my motif")))
  add_create(pcm_in, args)
add_create(matrix_input(list(c(5, 0, 0, 5), c(2, 2, 3, 3)), c("T", "G", "C", "A")),
           list(type = "PCM"))
add_create(matrix_input(list(c(0.1, 0.2, 0.3, 0.4), c(0.25, 0.25, 0.25, 0.25)), DNA_BASES))
add_create(matrix_input(list(c(1.2, -0.5, -3, 0.1), c(-1, 1, -1, 1)), DNA_BASES))
add_create(matrix_input(list(c(1.2, -0.5, -3, 0.1), c(-1, 1, -1, 1)), DNA_BASES),
           list(type = "PWM"))
add_create(matrix_input(list(c(0.5, -0.2, 0.1, 0), c(0, 0, 0.3, -0.1)), DNA_BASES),
           list(type = "CWM"))
add_create(matrix_input(list(c(0.2, 0.5, 0.3), c(0.1, 0.1, 0.8)), c("x", "y", "z")))
add_create(matrix_input(list(c(3, 1, 0, 0)), RNA_BASES), list(type = "PCM"))
add_create(matrix_input(lapply(c("A", "K", "W"), consensus_to_ppmAAC), AA_STANDARD2))

for (s in c("TATAWAW", "ACGTN", "MRWSYKVHDBN", "ACGU", "acgt",
            "ACDEFGHIKLMNPQRSTVWY", "BZJX", "123"))
  for (args in list(list(), list(type = "PCM"), list(type = "PWM", pseudocount = 1),
                    list(type = "ICM"), list(nsites = 50)))
    add_create(strings_input(s), args)
add_create(strings_input("abba"), list(alphabet = "abc"))
add_create(strings_input("ACGT"), list(alphabet = "RNA"))
add_create(strings_input("ACGT"), list(alphabet = "AA"))

seq_sets <- list(c("AAGT", "ACGT", "ACCT", "TCGT"), c("ACGU", "ACGA", "UCGU"),
                 c("ACDE", "ACDF", "GCDE"), c("xyz", "xzz", "yyz"),
                 c("ACGN", "ACGT"))
for (s in seq_sets)
  for (args in list(list(), list(type = "PCM"), list(type = "PWM", pseudocount = 1),
                    list(type = "ICM", pseudocount = 1), list(nsites = 10),
                    list(bkg = NULL)))
    add_create(strings_input(s), args[!vapply(args, is.null, logical(1))])
add_create(strings_input(c("AAGT", "ACGT")), list(bkg = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)))

write_ref(list(cases = create_cases), "create.json")

# ---------------------------------------------------------------------------
# convert_type()

data(examplemotif, examplemotif2, ArabidopsisMotif, package = NULL,
     envir = environment(), lib.loc = NULL)
jaspar <- read_jaspar(file.path(um_path, "inst", "extdata", "jaspar.txt"))
base_motifs <- c(
  list(examplemotif, examplemotif2, ArabidopsisMotif),
  jaspar,
  list(
    create_motif("ACGUA", alphabet = "RNA", nsites = 30, pseudocount = 1),
    create_motif(c("ACDE", "ACDF", "GCDE", "ACKE"), type = "PCM", pseudocount = 1),
    create_motif(matrix(c(0.5, -0.2, 0.1, 0, 0, 0, 0.3, -0.1), nrow = 4,
                        dimnames = list(DNA_BASES, NULL)), type = "CWM"),
    create_motif("TTGACWS"),
    create_motif(c("AAGT", "ACGT", "ACCT", "TCGT"), type = "PCM",
                 bkg = c(A = 0.5, C = 0, G = 0, T = 0.5), pseudocount = 1),
    create_motif(c("xyz", "xzz", "yyz"), type = "PCM")
  )
)
base_motifs <- lapply(base_motifs, function(m) { m@multifreq <- list(); m })

types <- c("PCM", "PPM", "PWM", "ICM", "CWM")
convert_cases <- list()
for (b in base_motifs) {
  sources <- list(b)
  for (t in types[1:4]) {
    r <- attempt(convert_type(b, t))
    if (r$ok) sources[[length(sources) + 1]] <- r$value
  }
  for (src in sources) {
    for (t in types) {
      variants <- list(list())
      if (t %in% c("PPM", "PWM", "ICM")) variants <- c(variants, list(list(pseudocount = 2)))
      if (t == "ICM") variants <- c(variants, list(list(relative_entropy = TRUE)))
      for (v in variants) {
        r <- attempt(do.call(convert_type, c(list(motifs = src, type = t), v)))
        convert_cases[[length(convert_cases) + 1]] <- c(
          list(input = motif_record(src), args = c(list(type = t), v)),
          as_motif_result(r)
        )
      }
    }
  }
}

write_ref(list(cases = convert_cases), "convert.json")

# ---------------------------------------------------------------------------
# Motif operations

ops_cases <- list()
add_op <- function(op, m, args, expr) {
  r <- attempt(expr)
  ops_cases[[length(ops_cases) + 1]] <<- c(
    list(op = op, input = motif_record(m), args = args),
    as_motif_result(r)
  )
}

dna_motifs <- c(list(examplemotif, ArabidopsisMotif), jaspar,
                list(convert_type(examplemotif, "PWM"),
                     convert_type(jaspar[[1]], "ICM"),
                     base_motifs[[length(jaspar) + 7]],
                     base_motifs[[length(jaspar) + 8]]))
dna_motifs <- lapply(dna_motifs, function(m) { m@multifreq <- list(); m })
rna_motif <- base_motifs[[length(jaspar) + 4]]
cwm_motif <- create_motif(matrix(c(0.01, -0.02, 0, 0.01,
                                   0.5, -0.2, 0.1, 0,
                                   0, 0, 0.3, -0.1,
                                   0.02, 0, 0, 0.01), nrow = 4,
                                 dimnames = list(DNA_BASES, NULL)), type = "CWM")

for (m in c(dna_motifs, list(rna_motif)))
  add_op("motif_rc", m, list(), motif_rc(m))
for (m in dna_motifs) {
  for (ic in c(0.25, 0.5, 1))
    for (tf in c("both", "left", "right"))
      add_op("trim_motifs", m, list(min_ic = ic, trim_from = tf),
             trim_motifs(m, min.ic = ic, trim.from = tf))
  for (tol in c(0.05, 0.1))
    add_op("round_motif", m, list(pct_tolerance = tol), round_motif(m, tol))
  add_op("normalize", m, list(), normalize(m))
  m0 <- m; m0@pseudocount <- 0
  add_op("normalize", m0, list(), normalize(m0))
}
for (th in c(0.1, 0.3, 0.9))
  for (tf in c("both", "left", "right"))
    add_op("trim_cwm", cwm_motif, list(trim_threshold = th, trim_from = tf),
           trim_cwm(cwm_motif, trim.threshold = th, trim.from = tf))
add_op("trim_cwm", cwm_motif, list(trim_threshold = 0.3, abs_threshold = 0.25,
                                   trim_from = "both"),
       trim_cwm(cwm_motif, abs.threshold = 0.25))

write_ref(list(cases = ops_cases), "ops.json")

# ---------------------------------------------------------------------------
# Motif file I/O
#
# Input files go to tests/reference/files: the universalmotif example files,
# plus small files covering cases the examples do not.

files_dir <- file.path(out_dir, "files")
dir.create(files_dir, showWarnings = FALSE)
extdata <- file.path(um_path, "inst", "extdata")
for (f in c("cisbp.txt", "dreme.txt", "hocomoco.txt", "homer.txt", "jaspar.txt",
            "meme_custom_alph.txt", "meme_full.txt", "meme_minimal.txt",
            "transfac.txt", "uniprobe_full.txt", "uniprobe_minimal.txt"))
  file.copy(file.path(extdata, f), file.path(files_dir, f), overwrite = TRUE)

put <- function(name, lines) writeLines(lines, file.path(files_dir, name))

put("jaspar_unlabelled.txt", c(
  ">MA0001.1 AGL3",
  "0 3 79 40 66 48 65 11 65 0",
  "94 75 4 3 1 2 5 2 3 3",
  "1 0 3 4 1 0 5 3 28 88",
  "2 19 11 50 29 47 22 81 1 6",
  ">MA0002.1",
  "  10 0 0",
  "  0 10 0",
  "  0 0 10",
  "  0 0 0"))
put("jaspar_rna.txt", c(">rna1 alt", "A [ 4 0 1 ]", "C [ 0 4 1 ]",
                        "G [ 0 0 1 ]", "U [ 0 0 1 ]"))
aa_rows <- vapply(seq_along(AA_STANDARD2), function(i) {
  paste(AA_STANDARD2[i], "[", paste((i + 0:2) %% 7, collapse = " "), "]")
}, character(1))
put("jaspar_aa.txt", c(">aa1 protein", aa_rows))
put("homer_family.txt", c(
  paste(">GCTGTGGTTT", "RUNX(Runt)/HPC7-Runx1-ChIP-Seq(GSE22178)/Homer",
        "6.509", "-1.104e+03", "0", "T:2981.0(39.61%),B:1949.1(4.13%),P:1e-479",
        sep = "\t"),
  "0.331\t0.226\t0.240\t0.203",
  "0.001\t0.001\t0.997\t0.001",
  "0.001\t0.997\t0.001\t0.001",
  paste(">ACGT", "plain", "4.2", sep = "\t"),
  "0.97\t0.01\t0.01\t0.01",
  "0.01\t0.97\t0.01\t0.01",
  "0.01\t0.01\t0.97\t0.01",
  "0.01\t0.01\t0.01\t0.97"))
put("transfac_meta.txt", c(
  "VV  TRANSFAC MATRIX TABLE",
  "XX",
  "//",
  "AC  M00001",
  "XX",
  "ID  V$MYOD_01",
  "XX",
  "NA  MyoD",
  "HC  bHLH",
  "OS  mouse",
  "P0      A      C      G      T",
  "01      1      2      2      0      S",
  "02      2      1      2      0      R",
  "03      3      0      1      1      A",
  "04      0      5      0      0      C",
  "05      5      0      0      0      A",
  "XX",
  "//",
  "AC  M00002",
  "NA  E47",
  "PO      A      C      G      T",
  "01    1.5    2.5      0      1",
  "02      0      5      0      0",
  "03      4      0      1      0",
  "XX",
  "//",
  "ID  only_id",
  "P0  A  C  G  T",
  "1  4  0  0  0",
  "2  0  4  0  0",
  "//"))
put("cisbp_nometa.txt", c(
  "Pos\tA\tC\tG\tT",
  "1\t0.1\t0.2\t0.3\t0.4",
  "2\t0.7\t0.1\t0.1\t0.1",
  "",
  "Pos\tA\tC\tG\tT",
  "1\t0.25\t0.25\t0.25\t0.25",
  "2\t0.0\t1.0\t0.0\t0.0"))
put("cisbp_rna.txt", c(
  "TF\tT1",
  "TF Name\tRBFOX",
  "Family\tRRM",
  "Species\tHomo_sapiens",
  "Pos\tA\tC\tG\tU",
  "1\t0.1\t0.2\t0.3\t0.4",
  "2\t0.7\t0.1\t0.1\t0.1"))
put("matrix_cols.txt", c(
  "first", "0.1 0.2 0.7", "0.2 0.2 0.1", "0.3 0.3 0.1", "0.4 0.3 0.1", "",
  "second", "1 2", "3 4", "5 6", "7 8"))
put("matrix_rownames.txt", c(
  "A 0.1 0.7", "C 0.2 0.1", "G 0.3 0.1", "T 0.4 0.1"))

read_cases <- list()
add_read <- function(format, file, args, expr) {
  r <- attempt(expr)
  case <- list(format = format, file = file, args = args)
  if (!r$ok) {
    case$error <- r$error
  } else {
    v <- r$value
    if (is(v, "universalmotif")) v <- list(v)
    if (is.list(v) && !is.null(v$motifs)) {
      mots <- v$motifs
      if (is(mots, "universalmotif")) mots <- list(mots)
      case$motifs <- lapply(unname(mots), motif_record)
      if (!is.null(v$sites)) {
        sites <- v$sites
        if (!is.list(sites)) sites <- list(sites)
        case$sites <- lapply(unname(sites), function(s)
          list(names = I(names(s)), sites = I(as.character(s))))
      }
      if (!is.null(v$sites.meta)) {
        sm <- v$sites.meta
        if (is.data.frame(sm)) sm <- list(sm)
        case$sites_meta <- lapply(unname(sm), function(d)
          list(sequence = I(d$Sequence), position = vec(d$Position),
               pvalue = vec(d$Pvalue), site = I(d$Site)))
      }
      if (!is.null(v$sites.meta.combined)) {
        d <- v$sites.meta.combined
        case$combined <- list(sequence = I(d$Sequence),
                              pvalue = vec(d$Combined.Pvalue),
                              diagram = I(d$Diagram))
      }
    } else {
      case$motifs <- lapply(unname(v), motif_record)
    }
  }
  read_cases[[length(read_cases) + 1]] <<- case
}
fp <- function(f) file.path(files_dir, f)

for (f in c("jaspar.txt", "jaspar_unlabelled.txt", "jaspar_rna.txt", "jaspar_aa.txt"))
  add_read("jaspar", f, list(), read_jaspar(fp(f)))
for (f in c("homer.txt", "homer_family.txt"))
  add_read("homer", f, list(), read_homer(fp(f)))
add_read("transfac", "transfac.txt", list(), read_transfac(fp("transfac.txt")))
# universalmotif cannot read files that start with the usual "VV ... //"
# header block, so its values come from a copy without the header.
tf_lines <- readLines(fp("transfac_meta.txt"))
put("transfac_meta_noheader.txt", tf_lines[-(1:3)])
add_read("transfac", "transfac_meta.txt", list(),
         read_transfac(fp("transfac_meta_noheader.txt")))
for (f in c("cisbp.txt", "cisbp_nometa.txt", "cisbp_rna.txt"))
  add_read("cisbp", f, list(), read_cisbp(fp(f)))
for (f in c("uniprobe_full.txt", "uniprobe_minimal.txt"))
  add_read("uniprobe", f, list(), read_uniprobe(fp(f)))
for (f in c("meme_full.txt", "meme_minimal.txt", "dreme.txt", "meme_custom_alph.txt"))
  add_read("meme", f, list(), read_meme(fp(f)))
add_read("meme", "meme_minimal.txt", list(cwm = TRUE), read_meme(fp("meme_minimal.txt"), CWM = TRUE))
for (f in c("meme_full.txt", "meme_custom_alph.txt"))
  add_read("meme", f, list(readsites = TRUE, readsites_meta = TRUE),
           read_meme(fp(f), readsites = TRUE, readsites.meta = TRUE))
add_read("matrix", "hocomoco.txt", list(headers = ">", positions = "rows"),
         read_matrix(fp("hocomoco.txt"), headers = ">", positions = "rows"))
add_read("matrix", "matrix_cols.txt", list(),
         read_matrix(fp("matrix_cols.txt")))
add_read("matrix", "matrix_cols.txt", list(type = "PPM"),
         read_matrix(fp("matrix_cols.txt"), type = "PPM"))
add_read("matrix", "matrix_rownames.txt", list(headers = FALSE, rownames = TRUE),
         read_matrix(fp("matrix_rownames.txt"), headers = FALSE, rownames = TRUE))

# Motifs to write. Each writer call becomes a case with the exact file text.
w_dna <- c(list(examplemotif, ArabidopsisMotif), jaspar[1:2],
           list(create_motif("TTGACWS", name = "with altname", altname = "a=b c",
                             family = "WRKY", organism = "Arabidopsis thaliana",
                             nsites = 37, pseudocount = 1, eval = 1.5e-12)))
w_dna <- lapply(w_dna, function(m) { m@multifreq <- list(); m })
w_dna[[2]]@bkg <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
w_rna <- list(create_motif(c("ACGU", "ACGA", "UCGU"), name = "rna", nsites = 3))
w_aa <- list(create_motif(c("ACDE", "ACDF", "GCDE", "ACKE"), type = "PCM", name = "aa"))
w_custom <- list(create_motif(c("ACGm", "ACMm", "ACGT"), name = "meth",
                              alphabet = "ACGMTm"))
w_cwm <- list(create_motif(matrix(c(0.5, -0.2, 0.1, 0, 0, 0, 0.3, -0.1), nrow = 4,
                                  dimnames = list(DNA_BASES, NULL)), type = "CWM",
                           name = "cwm"))
w_meta <- list(create_motif("ACGTWS", name = "full", altname = "alt", family = "fam",
                            organism = "org", nsites = 20, bkgsites = 400,
                            pval = 1e-8, qval = 1e-6, eval = 3e-5, strand = "+",
                            extrainfo = c(source = "test", score = "4.25"),
                            bkg = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3),
                            pseudocount = 1))

write_cases <- list()
add_write <- function(format, motifs, args, fn) {
  f <- tempfile()
  r <- attempt({ suppressWarnings(fn(motifs, f)); readChar(f, file.size(f), useBytes = TRUE) })
  case <- list(format = format, args = args, motifs = lapply(motifs, motif_record))
  if (r$ok) case$text <- r$value else case$error <- r$error
  write_cases[[length(write_cases) + 1]] <<- case
  unlink(f)
}

for (ms in list(w_dna, w_rna, w_aa, w_custom))
  add_write("meme", ms, list(), function(m, f) write_meme(m, f))
add_write("meme", w_dna[1:2], list(version = 4, strand = "+"),
          function(m, f) write_meme(m, f, version = 4, strand = "+"))
add_write("meme", w_dna[3:4], list(bkg = I(c(0.3, 0.2, 0.2, 0.3))),
          function(m, f) write_meme(m, f, bkg = c(0.3, 0.2, 0.2, 0.3)))
add_write("meme", w_cwm, list(cwm = TRUE), function(m, f) write_meme(m, f, CWM = TRUE))
add_write("meme", w_cwm, list(), function(m, f) write_meme(m, f))
add_write("meme", w_meta, list(), function(m, f) write_meme(m, f))

for (ms in list(w_dna, w_rna, w_aa, w_custom))
  add_write("jaspar", ms, list(), function(m, f) write_jaspar(m, f))

for (ms in list(w_dna, w_meta))
  add_write("homer", ms, list(), function(m, f) write_homer(m, f))
add_write("homer", w_dna[1:3], list(threshold = I(0.5)),
          function(m, f) write_homer(m, f, threshold = 0.5))
add_write("homer", w_dna[1:3], list(threshold = I(1e-4), threshold_type = "pvalue"),
          function(m, f) write_homer(m, f, threshold = 1e-4, threshold.type = "pvalue"))
add_write("homer", w_dna[1:3], list(threshold = I(c(5, 6.5, 7)), threshold_type = "logodds_abs"),
          function(m, f) write_homer(m, f, threshold = c(5, 6.5, 7), threshold.type = "logodds.abs"))

for (ms in list(w_dna, w_rna, w_meta))
  add_write("transfac", ms, list(), function(m, f) write_transfac(m, f))
add_write("transfac", w_dna[c(1, 5)], list(name_tag = "AC", altname_tag = "ID"),
          function(m, f) write_transfac(m, f, name.tag = "AC", altname.tag = "ID"))

matrix_args <- list(
  list(),
  list(positions = "rows"),
  list(rownames = TRUE),
  list(positions = "rows", rownames = TRUE),
  list(headers = FALSE),
  list(headers = ">"),
  list(sep = "//"),
  list(sep = NULL),
  list(digits = 3),
  list(type = "PCM")
)
for (a in matrix_args) {
  rec_args <- a
  if ("sep" %in% names(a) && is.null(a$sep)) rec_args <- list(sep_none = TRUE)
  add_write("matrix", w_dna[c(1, 3)], rec_args,
            function(m, f) do.call(write_matrix, c(list(motifs = m, file = f), a)))
}
add_write("matrix", w_aa, list(), function(m, f) write_matrix(m, f))

write_ref(list(read = read_cases, write = write_cases), "io.json")

# ---------------------------------------------------------------------------
# P-values (method = "dynamic") and sequence scanning

pv_motifs <- c(
  list(examplemotif, ArabidopsisMotif),
  jaspar[1:3],
  list(
    convert_type(examplemotif, "PWM", pseudocount = 1),
    create_motif("TTGACWS"),
    create_motif(c("ACGU", "ACGA", "UCGU"), nsites = 3, pseudocount = 1),
    create_motif(c("ACDE", "ACDF", "GCDE", "ACKE"), type = "PCM", pseudocount = 1),
    create_motif(c("xyz", "xzz", "yyz"), type = "PCM", pseudocount = 1),
    create_motif(matrix(c(0.5, -0.2, 0.1, 0, 0, 0, 0.3, -0.1), nrow = 4,
                        dimnames = list(DNA_BASES, NULL)), type = "CWM")
  )
)
pv_motifs <- lapply(pv_motifs, function(m) { m@multifreq <- list(); m })
pvals_in <- c(1, 0.9, 0.5, 0.1, 0.01, 1e-3, 1e-4, 1e-6, 1e-8, 1e-12, 1e-300)

pv_cases <- list()
for (m in pv_motifs) {
  pwm <- tryCatch(suppressMessages(convert_type(normalize(m), "PWM"))@motif,
                  error = function(e) convert_type(m, "PWM")@motif)
  pwm[!is.finite(pwm)] <- -10
  lo <- sum(apply(pwm, 2, min)); hi <- sum(apply(pwm, 2, max))
  scores <- c(lo - 1, seq(lo, hi, length.out = 11), hi + 1, 0, 1.001, -1.001)
  bkgs <- list(NULL)
  if (m@alphabet == "DNA") bkgs <- c(bkgs, list(c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)))
  for (b in bkgs) {
    args <- list()
    if (!is.null(b)) args$bkg <- vec(b)
    r1 <- attempt(if (is.null(b)) motif_pvalue(m, score = scores)
                  else motif_pvalue(m, score = scores, bkg.probs = b))
    r2 <- attempt(if (is.null(b)) motif_pvalue(m, pvalue = pvals_in)
                  else motif_pvalue(m, pvalue = pvals_in, bkg.probs = b))
    pv_cases[[length(pv_cases) + 1]] <- list(
      input = motif_record(m), args = args,
      scores = vec(scores), pvalues = if (r1$ok) vec(r1$value) else NULL,
      error = if (!r1$ok) r1$error else NULL,
      pvalue_in = vec(pvals_in), score_out = if (r2$ok) vec(r2$value) else NULL
    )
  }
}
write_ref(list(cases = pv_cases), "pvalue.json")

data(ArabidopsisPromoters, package = NULL, envir = environment(), lib.loc = NULL)
proms <- ArabidopsisPromoters[1:12]
Biostrings::writeXStringSet(proms, fp("promoters.fa"))
extra <- Biostrings::DNAStringSet(c(
  with_n = "ACGTTTGACCAAAAAAAANNNNNAAAAAAAAAAAAGTCAAACGT",
  short = "ACGT",
  plain = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"))
Biostrings::writeXStringSet(extra, fp("extra.fa"))

scan_motifs <- c(list(ArabidopsisMotif, examplemotif), jaspar[1:3],
                 list(create_motif("TTGACWS")))
scan_motifs <- lapply(scan_motifs, function(m) { m@multifreq <- list(); m })
scan_cases <- list()
add_scan <- function(file, seqs, args) {
  r <- attempt(do.call(scan_sequences_lite,
                       c(list(motifs = scan_motifs, sequences = seqs,
                              return.granges = FALSE), args)))
  case <- list(file = file, args = args)
  if (r$ok) {
    d <- r$value
    case$hits <- list(
      motif = I(d$motif), motif_i = I(d$motif.i), sequence = I(d$sequence),
      sequence_i = I(d$sequence.i), start = I(d$start), end = I(d$end),
      strand = I(d$strand), score = vec(d$score), score_pct = vec(d$score.pct),
      match = I(d$match), pvalue = vec(d$pvalue))
  } else case$error <- r$error
  scan_cases[[length(scan_cases) + 1]] <<- case
}
add_scan("promoters.fa", proms, list())
add_scan("promoters.fa", proms, list(pvalue = 1e-3))
add_scan("promoters.fa", proms, list(pvalue = 1e-3, RC = FALSE))
add_scan("promoters.fa", proms, list(pvalue = 1e-3, no.overlaps = TRUE))
add_scan("promoters.fa", proms, list(pvalue = 1e-3, no.overlaps = TRUE,
                                     no.overlaps.by = "score"))
add_scan("promoters.fa", proms, list(pvalue = 1e-3, no.overlaps = TRUE,
                                     no.overlaps.by.strand = TRUE,
                                     no.overlaps.by.motif = TRUE))
add_scan("extra.fa", extra[c(1, 3)], list(pvalue = 0.05))
write_ref(list(motifs = lapply(scan_motifs, motif_record), cases = scan_cases), "scan.json")
