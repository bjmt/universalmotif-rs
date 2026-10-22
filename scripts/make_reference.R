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

# ---------------------------------------------------------------------------
# Comparison, merging and trees (lite functions)

homer_m <- read_homer(fp("homer.txt"))
meme_m <- read_meme(fp("meme_full.txt"))
transfac_m <- read_transfac(fp("transfac.txt"))
cmp_dna <- c(jaspar, list(examplemotif, ArabidopsisMotif), homer_m, meme_m,
             transfac_m[1:2])
cmp_dna <- lapply(cmp_dna, function(m) { m@multifreq <- list(); m })
cmp_rna <- list(
  create_motif(c("ACUU", "ACUA", "UCUC"), name = "r1", nsites = 3),
  create_motif("UUGACUAA", name = "r2", nsites = 20, pseudocount = 1),
  create_motif("CCGUAAGCU", name = "r3", nsites = 20, pseudocount = 1))
rows <- function(m) lapply(seq_len(nrow(m)), function(i) vec(m[i, ]))

cmp_cases <- list()
add_cmp <- function(set, motifs, args) {
  r <- attempt(do.call(compare_motifs_lite, c(list(motifs = motifs), args)))
  case <- list(set = set, args = args)
  if (!r$ok) {
    case$error <- r$error
  } else if (is.matrix(r$value)) {
    case$matrix <- rows(r$value)
  } else {
    d <- r$value
    case$long <- list(subject_i = I(d$subject.i), target_i = I(d$target.i),
                      offset = I(d$offset), strand = I(d$strand),
                      overlap = I(d$overlap), score = vec(d$score),
                      subject_consensus = I(d$subject.consensus),
                      target_consensus = I(d$target.consensus),
                      pvalue = vec(d$pvalue), qvalue = vec(d$qvalue))
  }
  cmp_cases[[length(cmp_cases) + 1]] <<- case
}
for (out in c("score", "pvalue", "qvalue"))
  for (rc in c(TRUE, FALSE))
    add_cmp("dna", cmp_dna, list(matrix.out = out, RC = rc))
for (mo in c(3L, 12L)) add_cmp("dna", cmp_dna, list(matrix.out = "pvalue", min.overlap = mo))
add_cmp("dna", cmp_dna, list(matrix.out = "pvalue", null = "parametric"))
add_cmp("dna", cmp_dna, list(matrix.out = "qvalue", null = "parametric",
                             bkg = c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)))
add_cmp("dna", cmp_dna, list(compare.to = c(1L, 7L), qvalue = 1))
add_cmp("dna", cmp_dna, list(compare.to = c(3L, 9L, 12L)))
add_cmp("dna", cmp_dna, list(compare.to = 2L, qvalue = 1, RC = FALSE, null = "parametric"))
add_cmp("rna", cmp_rna, list(matrix.out = "pvalue"))
add_cmp("rna", cmp_rna, list(compare.to = 1L, qvalue = 1))

merge_cases <- list()
add_merge <- function(set, idx, args) {
  ms <- (if (set == "dna") cmp_dna else cmp_rna)[idx]
  r <- attempt(do.call(merge_motifs_lite, c(list(motifs = ms), args)))
  merge_cases[[length(merge_cases) + 1]] <<- c(
    list(set = set, idx = I(idx), args = args), as_motif_result(r))
}
add_merge("dna", 1:3, list())
add_merge("dna", 1:3, list(weighted = TRUE))
add_merge("dna", 1:3, list(RC = FALSE, new.name = "joined"))
add_merge("dna", c(6, 8, 9, 10), list())
add_merge("dna", c(6, 8, 9, 10), list(min.overlap = 8L))
add_merge("dna", c(13, 14, 15), list(weighted = TRUE))
add_merge("dna", 4L, list(new.name = "single"))
add_merge("rna", 1:3, list(min.overlap = 3L))

similar_cases <- list()
for (args in list(list(), list(qvalue = 0.1), list(qvalue = 0.5),
                  list(qvalue = 0.1, null = "parametric", weighted = TRUE),
                  list(qvalue = 0.1, RC = FALSE, min.overlap = 6L))) {
  cl <- do.call(merge_similar_lite, c(list(motifs = cmp_dna, return.clusters = TRUE), args))
  ms <- do.call(merge_similar_lite, c(list(motifs = cmp_dna), args))
  similar_cases[[length(similar_cases) + 1]] <- list(
    args = args, clusters = I(cl$cluster), motifs = lapply(ms, motif_record))
}

tree_cases <- list()
methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median",
             "centroid", "ward.D2")
add_tree <- function(source, dmat, method, labels, extra = list()) {
  h <- hclust(as.dist(dmat), method = method)
  h$labels <- labels
  tree_cases[[length(tree_cases) + 1]] <<- c(list(
    source = source, method = method, dist = rows(dmat), labels = I(labels),
    merge = rows(h$merge), height = vec(h$height), order = I(h$order),
    newick = ape::write.tree(ape::as.phylo(h))), extra)
}
for (rc in c(TRUE, FALSE)) for (mo in c(6L, 4L)) {
  s <- compare_motifs_lite(cmp_dna, matrix.out = "score", RC = rc, min.overlap = mo)
  for (meth in c("complete", "average"))
    add_tree("motifs", (1 - s) / 2, meth, rownames(s),
             list(args = list(RC = rc, min.overlap = mo)))
}
set.seed(42)
for (k in 1:3) {
  pts <- matrix(round(runif(14 * 2) * 4), ncol = 2)  # many tied distances
  dmat <- as.matrix(dist(pts))
  for (meth in methods)
    add_tree("random", dmat, meth, paste0("leaf ", seq_len(nrow(dmat))))
}

write_ref(list(dna = lapply(cmp_dna, motif_record), rna = lapply(cmp_rna, motif_record),
               compare = cmp_cases, merge = merge_cases, similar = similar_cases,
               tree = tree_cases), "compare.json")

# ---------------------------------------------------------------------------
# k-lets, backgrounds, complexity and motif discovery
#
# Shuffling and sequence generation use R's random number generator, which
# the port does not reproduce, so those are covered by property tests in the
# crate instead.

df_cols <- function(d) {
  d <- as.data.frame(d)
  lapply(d, function(x) if (is.numeric(x)) vec(x) else I(as.character(x)))
}

klet_cases <- list(
  list(letters = I(c("A", "C", "G")), k = 2L, klets = I(get_klets(c("A", "C", "G"), 2))),
  list(letters = I(c("A", "C", "G", "T")), k = 3L,
       klets = I(get_klets(c("A", "C", "G", "T"), 3)))
)
count_cases <- list()
for (x in list(list("AACGTTAACG", 1L, "ACGT"), list("AACGTTAACG", 2L, "ACGT"),
               list("AACGTTAACG", 3L, "ACG"), list("MKVLAAGLLK", 2L, NULL))) {
  r <- if (is.null(x[[3]])) count_klets(x[[1]], x[[2]]) else count_klets(x[[1]], x[[2]], x[[3]])
  count_cases[[length(count_cases) + 1]] <- list(
    string = x[[1]], k = x[[2]], alph = x[[3]], klets = I(r$klets), counts = vec(r$counts))
}

window_cases <- list()
for (x in list(c(10, 3, 0), c(10, 3, 1), c(10, 4, 2), c(11, 4, 1), c(100, 7, 3),
               c(5, 5, 0), c(9, 3, 2))) {
  for (incomp in c(TRUE, FALSE)) {
    w <- calc_windows(x[1], x[2], x[3], incomp)
    window_cases[[length(window_cases) + 1]] <- list(
      n = x[1], window = x[2], overlap = x[3], incomplete = incomp,
      start = vec(w$start), stop = vec(w$stop))
  }
}

cx_strings <- c("ACGTACGTAC", "AAAAAAAAAT", "ACGTTGCAACGGTCAGTCGATCGATTTAGC", "AAAA",
                "MKVLAAGLLKWQRSTEE", "ATATATATATATCGCGCGCG")
cx_methods <- c("WoottonFederhen", "WoottonFederhenFast", "Trifonov", "TrifonovFast", "DUST")
complexity_cases <- list()
for (s in cx_strings) for (meth in cx_methods) {
  if (meth == "DUST" && grepl("[^ACGT]", s)) next
  for (alph in list(NULL, "ACGT")) {
    if (!is.null(alph) && grepl("[^ACGT]", s)) next
    for (mw in c(7L, 3L)) {
      if (mw != 7L && !startsWith(meth, "Trifonov")) next
      v <- calc_complexity(s, meth, alph, mw)
      complexity_cases[[length(complexity_cases) + 1]] <- list(
        string = s, method = meth, alph = alph, max_word = mw, value = v)
    }
  }
}

set.seed(7)
bkg_seqs <- Biostrings::DNAStringSet(vapply(1:4, function(i)
  paste(sample(c("A", "C", "G", "T"), 24, TRUE, prob = c(0.4, 0.1, 0.2, 0.3)), collapse = ""),
  character(1)))
names(bkg_seqs) <- c("s3", "s1", "s10", "s2")
cx_seqs <- Biostrings::DNAStringSet(c(
  one = "ACGTTGCAACGGTCAGTCGATCGATTTAGCAAAAAAAAAAAAAAACGT",
  two = "ATATATATATATCGCGCGCGTTTTTGGG",
  three = "ACGTACGTAC"))

bkg_cases <- list()
add_bkg <- function(args) {
  r <- attempt(do.call(get_bkg, c(list(sequences = bkg_seqs), args)))
  bkg_cases[[length(bkg_cases) + 1]] <<- c(list(args = args),
    if (r$ok) list(table = df_cols(r$value)) else list(error = r$error))
}
add_bkg(list())
add_bkg(list(k = 1:2, merge.res = FALSE))
add_bkg(list(k = 2L, RC = TRUE))
add_bkg(list(k = 1:2, pseudocount = 1L))
add_bkg(list(k = 1L, alphabet = "ACGTN"))
add_bkg(list(k = 1:2, window = TRUE, window.size = 6, window.overlap = 2))
add_bkg(list(k = 1L, window = TRUE, window.size = 0.25, window.overlap = 0))
add_bkg(list(k = 1:2, window = TRUE, merge.res = FALSE, window.size = 8,
             window.overlap = 0.1))
add_bkg(list(k = 1L, RC = TRUE, merge.res = FALSE, pseudocount = 2L))

meme_bkg <- paste(capture.output(get_bkg(bkg_seqs, k = 1:3, pseudocount = 1,
                                          to.meme = stdout())), collapse = "\n")

seqcx_cases <- list()
for (args in list(list(), list(window.size = 10, window.overlap = 0),
                  list(window.size = 12, method = "Trifonov"),
                  list(window.size = 15, window.overlap = 5, method = "TrifonovFast",
                       trifonov.max.word.size = 4),
                  list(window.size = 8, method = "WoottonFederhenFast"),
                  list(window.size = 20, method = "DUST"))) {
  r <- do.call(sequence_complexity, c(list(seqs = cx_seqs), args))
  seqcx_cases[[length(seqcx_cases) + 1]] <- list(args = args, table = df_cols(r))
}

# motif_finder() with explicit negatives. Two motifs are planted in 60 of 80
# random sequences, and the negatives are shuffled in R so the crate gets the
# exact same input.
set.seed(11)
rand_seq <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
plant <- function(s, m) {
  p <- sample.int(nchar(s) - nchar(m), 1)
  substr(s, p, p + nchar(m) - 1) <- m
  s
}
mf_pos <- vapply(1:80, function(i) {
  s <- rand_seq(150)
  if (i <= 60) s <- plant(s, sample(c("TTGACGTCAA", "TTGACGTCAT", "ATGACGTCAA"), 1))
  if (i %% 2 == 0) s <- plant(s, "GGGGCCCAT")
  s
}, character(1))
mf_neg <- as.character(shuffle_sequences(Biostrings::DNAStringSet(mf_pos), k = 2,
                                         rng.seed = 3))
finder_cases <- list()
for (args in list(list(min.width = 6L, max.width = 10L, qvalue = 1),
                  list(min.width = 8L, max.width = 8L, nmotifs = 2L, qvalue = 1, RC = FALSE),
                  list(min.width = 5L, max.width = 12L, hit.pvalue = 1e-3,
                       dedup.overlap = 0.2, pseudocount = 4L),
                  list(min.width = 6L, max.width = 7L, stop.pvalue = 1e-30))) {
  r <- do.call(motif_finder, c(list(sequences = Biostrings::DNAStringSet(mf_pos),
                                    bkg.sequences = Biostrings::DNAStringSet(mf_neg)), args))
  # to_list() copies the extra data frame columns into extrainfo, which the
  # crate's motifs leave empty because the statistics are returned alongside.
  ms <- lapply(if (nrow(r)) to_list(r) else list(), function(m) {
    m@extrainfo <- character()
    m
  })
  stats <- as.data.frame(r)[, c("rank", "width", "seqs_pos", "seqs_neg", "sites_pos",
                                "sites_neg", "n_pos", "n_neg", "pvalue", "qvalue")]
  finder_cases[[length(finder_cases) + 1]] <- list(
    args = args, motifs = lapply(ms, motif_record), stats = df_cols(stats))
}

write_ref(list(klets = klet_cases, count = count_cases, windows = window_cases,
               complexity = complexity_cases,
               bkg_sequences = I(as.character(bkg_seqs)), bkg_names = I(names(bkg_seqs)),
               bkg = bkg_cases, meme_bkg = meme_bkg,
               complexity_sequences = I(as.character(cx_seqs)),
               complexity_names = I(names(cx_seqs)), sequence_complexity = seqcx_cases,
               finder_positives = I(mf_pos), finder_negatives = I(mf_neg),
               finder = finder_cases), "sequences.json")

# ---------------------------------------------------------------------------
# Statistics, cutree, enrichment, co-occurrence and positional enrichment

stat_cases <- list(
  phyper = lapply(list(c(3, 10, 7, 8), c(0, 5, 5, 5), c(12, 30, 2000, 40),
                       c(40, 50, 60, 70), c(900, 1000, 99000, 1200),
                       c(-1, 4, 4, 4), c(5, 5, 5, 5), c(2, 400000, 3e6, 30),
                       c(5, 20, 30, 25)), function(a) {
    list(x = a[1], m = a[2], n = a[3], k = a[4],
         lower = phyper(a[1], a[2], a[3], a[4]),
         upper = phyper(a[1], a[2], a[3], a[4], lower.tail = FALSE),
         density = if (a[1] >= 0) dhyper(a[1], a[2], a[3], a[4]) else 0)
  }),
  pbinom = lapply(list(c(0, 10, 0.3), c(5, 10, 0.3), c(60, 100, 0.2), c(1, 1000, 1e-5),
                       c(400, 1000, 0.41), c(2, 2, 1), c(150, 150, 0.9),
                       c(3000, 10000, 0.29), c(80, 500, 0.02), c(11, 10, 0.5)),
                  function(a) list(k = a[1], n = a[2], p = a[3],
                                   upper = pbinom(a[1] - 1, a[2], a[3], lower.tail = FALSE),
                                   density = dbinom(a[1], a[2], a[3]))),
  fisher = lapply(list(c(10, 3, 2, 15), c(1, 1, 1, 1), c(0, 5, 5, 0), c(30, 70, 50, 50),
                       c(500, 120, 99500, 99880), c(81, 2, 1, 80), c(3, 0, 0, 3)),
                  function(a) list(table = I(a), pvalue = fisher.test(
                    matrix(a, 2, byrow = TRUE), alternative = "greater")$p.value)),
  bonferroni = list(p = I(c(0.01, 0.2, 0.003, 0.5)),
                    adjusted = I(p.adjust(c(0.01, 0.2, 0.003, 0.5), "bonferroni")))
)

set.seed(5)
cut_cases <- list()
for (meth in c("complete", "average", "single", "ward.D2")) {
  pts <- matrix(round(runif(12 * 2) * 5), ncol = 2)
  dmat <- as.matrix(dist(pts))
  h <- hclust(as.dist(dmat), method = meth)
  hs <- c(0, sort(unique(h$height))[c(1, 3, 5)], mean(range(h$height)), 100)
  cut_cases[[length(cut_cases) + 1]] <- list(
    method = meth, dist = rows(dmat),
    k = lapply(1:12, function(k) I(unname(cutree(h, k)))),
    h = I(hs), by_h = lapply(hs, function(x) I(unname(cutree(h, h = x)))))
}

# Motifs planted into random sequences, plus a background drawn in R so the
# crate gets the exact same input.
set.seed(21)
en_motifs <- list(
  create_motif("TTGACGTCAA", nsites = 20, pseudocount = 1, name = "creb"),
  create_motif("GGGGCGGGG", nsites = 20, pseudocount = 1, name = "sp1"),
  create_motif("CACGTG", nsites = 20, pseudocount = 1, name = "ebox"),
  create_motif("TATAAAA", nsites = 20, pseudocount = 1, name = "tata"),
  create_motif("TGASTCA", nsites = 30, pseudocount = 1, name = "ap1"),
  create_motif("GATAAG", nsites = 20, pseudocount = 1, name = "creb")
)
rand_dna <- function(n, len) Biostrings::DNAStringSet(
  vapply(seq_len(n), function(i) paste(sample(c("A", "C", "G", "T"), len, TRUE),
                                       collapse = ""), character(1)))
en_tgt <- rand_dna(60, 200)
en_tgt <- implant_motifs(en_motifs[c(1, 3)], en_tgt, n.per.seq = 1)
en_tgt[1:30] <- implant_motifs(en_motifs[2], en_tgt[1:30], n.per.seq = 2,
                               min.spacing = 3)
en_bkg <- rand_dna(70, 180)
enrich_cases <- list()
for (args in list(list(qvalue = 1), list(qvalue = 1, test = "sites"),
                  list(qvalue = 1, RC = FALSE, pvalue = 1e-3),
                  list(qvalue = 1, test = "sites", pseudocount = 0, pvalue = 1e-3),
                  list(qvalue = 0.05, pseudocount = 3))) {
  r <- do.call(enrich_motifs_lite, c(list(motifs = en_motifs, sequences = en_tgt,
                                          bkg.sequences = en_bkg), args))
  enrich_cases[[length(enrich_cases) + 1]] <- list(args = args, table = df_cols(r))
}

# Co-occurrence: through the internal scan, and from a hit table.
co_cases <- list()
add_coocc <- function(args) {
  r <- do.call(motif_coocc, c(list(motifs = en_motifs, sequences = en_tgt), args))
  co_cases[[length(co_cases) + 1]] <<- list(args = args, table = df_cols(r))
}
add_coocc(list(pvalue = 1e-3))
add_coocc(list(pvalue = 1e-3, max.distance = 50L))
add_coocc(list(pvalue = 1e-3, self.pairs = TRUE, max.distance = 100L, RC = FALSE))
add_coocc(list(pvalue = 1e-2, min.coocc = 5L, pseudocount = 1L))
add_coocc(list(pvalue = 1e-3, no.overlaps = TRUE, max.distance = 30L, self.pairs = TRUE))
co_hits <- data.frame(
  motif.i = c(1, 2, 1, 2, 3, 1, 3, 2, 2, 4, 1, 2, 3, 1, 1),
  sequence.i = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 8),
  start = c(10, 25, 5, 400, 7, 100, 150, 30, 90, 12, 60, 61, 3, 15, 15))
co_hit_cases <- lapply(list(list(), list(max.distance = 20L),
                            list(self.pairs = TRUE, max.distance = 1000L),
                            list(pseudocount = 2L, min.coocc = 0L)), function(args) {
  r <- do.call(motif_coocc, c(list(motifs = en_motifs[1:4], hits = co_hits,
                                   n.sequences = 10L), args))
  list(args = args, table = df_cols(r))
})

# Positional enrichment: a motif planted near the middle of centred
# sequences, another planted near one end, and the rest at random.
set.seed(33)
pk_seqs <- rand_dna(150, 300)
pk_seqs <- implant_motifs(en_motifs[1], pk_seqs, positions = lapply(1:150, function(i)
  if (i <= 100) 146L + sample.int(8, 1) else integer(0)))
pk_seqs <- implant_motifs(en_motifs[3], pk_seqs, positions = lapply(1:150, function(i)
  if (i %% 2 == 0) 40L + sample.int(20, 1) else integer(0)))
names(pk_seqs) <- sprintf("s%03d", seq_along(pk_seqs))
pk_hits <- scan_sequences_lite(en_motifs[c(1, 3, 4, 5)], pk_seqs, pvalue = 1e-3,
                               return.granges = FALSE)
pk_strand <- setNames(ifelse(seq_along(pk_seqs) %% 3 == 0, "-", "+"), names(pk_seqs))
peak_cases <- list()
for (args in list(list(qvalue = 1), list(qvalue = 1, mode = "local"),
                  list(qvalue = 1, mode = "local", seq.strand = pk_strand,
                       position.step = 3L, min.window = 6L, window.step = 4L),
                  list(qvalue = 0.05, max.window = 300L, window.step = 25L),
                  list(qvalue = 1, min.window = 1L, max.window = 40L, window.step = 1L))) {
  r <- do.call(motif_peaks, c(list(hits = pk_hits, seq.length = 300L), args))
  t <- df_cols(r[, setdiff(names(r), "centers")])
  t$centers <- lapply(r$centers, function(x) I(x))
  if (!is.null(args$seq.strand)) args$seq.strand <- as.list(args$seq.strand)
  peak_cases[[length(peak_cases) + 1]] <- list(args = args, table = t)
}

write_ref(list(stats = stat_cases, cutree = cut_cases,
               motifs = lapply(en_motifs, motif_record),
               targets = I(as.character(en_tgt)), background = I(as.character(en_bkg)),
               enrich = enrich_cases, coocc = co_cases,
               coocc_hits = df_cols(co_hits), coocc_hit_cases = co_hit_cases,
               peak_sequences = I(as.character(pk_seqs)),
               peak_hits = df_cols(pk_hits), peaks = peak_cases), "enrich.json")

# ---------------------------------------------------------------------------
# Plotting: logo polygons and the aligned matrices behind multi-motif logos

poly_cols <- function(d) {
  d <- as.data.frame(d)
  keep <- !is.na(d$x)
  list(x = vec(d$x[keep]), y = vec(d$y[keep]), group = I(as.character(d$group[keep])),
       letter_id = I(as.character(d$letter.id[keep])),
       motif = if ("motif.id" %in% names(d)) I(as.character(d$motif.id[keep])))
}
lmat <- function(vals, rows) matrix(vals, nrow = length(rows), dimnames = list(rows, NULL))
logo_mats <- list(
  plain = lmat(c(0.1, 0.6, 0.2, 0.1, 0.25, 0.25, 0.25, 0.25, 0.9, 0, 0.05, 0.05, 0, 0, 0, 0),
               c("A", "C", "G", "T")),
  signed = lmat(c(1.2, -0.8, 0.3, -1.5, -0.2, 0.9, -0.01, 0.4, 0.5, 0.5, -0.6, -0.6),
                c("A", "C", "G", "T")),
  multi = lmat(c(0.5, 1, 0.25, 0.3, 0, 2, 0.4, 0.7), c("AA", "C", "GT", "T")),
  tall = lmat(c(3, 1, 0.5, 0.2, 0.02, 2, 2, 2), c("A", "C", "G", "T"))
)
logo_args <- list(list(), list(sort.positions = TRUE), list(flip.neg = TRUE, sort.positions = TRUE),
                  list(min.height = 0.1, x.spacer = 0.1, y.spacer = 0.05),
                  list(fit.to.height = 1, sort.positions = TRUE, sort.positions.decreasing = FALSE))
logo_cases <- list()
for (nm in names(logo_mats)) for (a in logo_args) {
  m <- logo_mats[[nm]]
  d <- suppressWarnings(do.call(view_logo, c(list(m), a))$data)
  logo_cases[[length(logo_cases) + 1]] <- list(
    matrix = nm, args = a, polygons = poly_cols(d))
}

vm <- list(
  create_motif("TTGACGTCAA", nsites = 20, pseudocount = 1, name = "creb"),
  create_motif("TGACGTCA", nsites = 30, pseudocount = 1, name = "creb2"),
  create_motif(c("TTGACGTC", "TTGACCTC", "ATGACGTA"), name = "sites"),
  create_motif("ATTGTCAAA", nsites = 15, pseudocount = 1, name = "rc"),
  create_motif("GGGGCCCCAT", nsites = 25, pseudocount = 1, name = "other"),
  create_motif("TTGACG", nsites = 10, pseudocount = 1, name = "creb")
)
vm[[3]] <- convert_type(vm[[3]], "PCM")
vm[[5]]["bkg"] <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
vm_rna <- lapply(vm[1:3], function(m) switch_alph(m))
raw_cases <- list()
add_raw <- function(set, idx, args) {
  ms <- (if (set == "dna") vm else vm_rna)[idx]
  r <- suppressWarnings(do.call(view_motifs_lite, c(list(ms, return.raw = TRUE), args)))
  raw_cases[[length(raw_cases) + 1]] <<- list(
    set = set, idx = I(idx), args = args, names = I(names(r)),
    matrices = unname(lapply(r, function(m) lapply(seq_len(ncol(m)), function(j) vec(m[, j])))))
}
add_raw("dna", 1:6, list())
add_raw("dna", 1:6, list(use.type = "PPM", sort.by = "ic"))
add_raw("dna", 1:6, list(use.type = "PWM", sort.by = "similarity"))
add_raw("dna", 1:6, list(use.type = "PCM", tryRC = FALSE))
add_raw("dna", 1:6, list(relative_entropy = TRUE, min.overlap = 4))
add_raw("dna", c(2, 4), list(min.overlap = 12))
add_raw("dna", 3L, list(use.type = "PCM"))
add_raw("rna", 1:3, list(sort.by = "similarity", RC.text = "(rc)"))

single_cases <- list()
for (ut in c("ICM", "PPM", "PCM", "PWM")) for (i in c(1, 3, 5)) {
  r <- attempt(view_motifs_lite(vm[[i]], use.type = ut)$data)
  single_cases[[length(single_cases) + 1]] <- c(list(idx = i, use_type = ut),
    if (r$ok) list(polygons = poly_cols(r$value)) else list(error = r$error))
}
multi_data <- poly_cols(view_motifs_lite(vm[1:4], sort.by = "ic")$data)

write_ref(list(logo_matrices = lapply(logo_mats, function(m) list(
                 rows = I(rownames(m)), columns = lapply(seq_len(ncol(m)), function(j) vec(m[, j])))),
               logo = logo_cases, motifs = lapply(vm, motif_record),
               raw = raw_cases, single = single_cases, multi = multi_data), "plot.json")
