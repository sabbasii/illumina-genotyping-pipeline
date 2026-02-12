#!/usr/bin/env Rscript
# scripts/eqtl/utils/inspect_helpers.R
# Shared helpers for scripts/eqtl/05_inspect_results/*.R
#
# Path convention:
#   REPO_ROOT/output/eqtl/

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
get_eqtl_paths <- function() {
  REPO_ROOT <- Sys.getenv("REPO_ROOT")
  if (is.null(REPO_ROOT) || REPO_ROOT == "") {
    stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")
  }

  base_eqtl   <- file.path(REPO_ROOT, "output", "eqtl")
  results_dir <- file.path(base_eqtl, "results")
  inspect_dir <- file.path(results_dir, "inspect")

  qc_dir        <- file.path(inspect_dir, "qc")
  manhattan_dir <- file.path(inspect_dir, "manhattan")
  pde_dir       <- file.path(inspect_dir, "pde")

  # Ensure dirs exist
  dirs <- c(results_dir, inspect_dir, qc_dir, manhattan_dir, pde_dir)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  list(
    REPO_ROOT = REPO_ROOT,
    base_eqtl = base_eqtl,
    results_dir = results_dir,
    inspect_dir = inspect_dir,
    qc_dir = qc_dir,
    manhattan_dir = manhattan_dir,
    pde_dir = pde_dir,

    # Common inputs
    eqtl_all   = file.path(results_dir, "eqtl_all.tsv"),
    eqtl_cis   = file.path(results_dir, "eqtl_cis.tsv"),
    eqtl_trans = file.path(results_dir, "eqtl_trans.tsv"),
    snpsloc    = file.path(base_eqtl, "snpsloc.txt"),
    geneloc    = file.path(base_eqtl, "geneloc.txt"),
    SNP_overlap = file.path(base_eqtl, "SNP_overlap.txt"),
    GE_overlap  = file.path(base_eqtl, "GE_overlap.txt")
  )
}

assert_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  invisible(path)
}

# -----------------------------
# Column pickers (robust)
# -----------------------------
pick_p_col <- function(dt) {
  # MatrixEQTL output sometimes uses `p-value`, sometimes `p.value`
  cand <- c("p-value", "p.value", "pval", "P", "p")
  found <- cand[cand %in% names(dt)][1]
  if (is.na(found)) {
    stop("No p-value column found. Expected one of: ", paste(cand, collapse = ", "))
  }
  found
}

pick_fdr_col <- function(dt) {
  cand <- c("FDR", "fdr", "qvalue", "q.value", "q-value")
  found <- cand[cand %in% names(dt)][1]
  if (is.na(found)) {
    stop("No FDR/q-value column found. Expected one of: ", paste(cand, collapse = ", "))
  }
  found
}

pick_beta_col <- function(dt) {
  cand <- c("beta", "Beta", "BETA", "slope", "Slope", "effect", "Effect")
  found <- cand[cand %in% names(dt)][1]
  if (is.na(found)) {
    stop("No beta/effect column found. Expected one of: ", paste(cand, collapse = ", "))
  }
  found
}

# -----------------------------
# snpsloc standardization
# -----------------------------
read_snpsloc <- function(path) {
  assert_file(path)
  snploc <- fread(path)

  if (ncol(snploc) < 3) stop("snpsloc must have at least 3 columns (snpid, chr, pos).")

  # If already has SNP/CHR/POS, keep; otherwise rename first 3.
  if (!all(c("SNP", "CHR", "POS") %in% names(snploc))) {
    setnames(snploc, names(snploc)[1:3], c("SNP", "CHR", "POS"))
  }

  snploc[, CHR := as.character(CHR)]
  snploc[, CHR := gsub("^chr", "", CHR, ignore.case = TRUE)]

  # Parse chromosome as integer (quietly); keep autosomes 1-22 by default
  snploc[, CHR_NUM := as.integer(CHR)]
  snploc[, POS := as.numeric(POS)]  # numeric to avoid downstream integer overflow

  snploc <- snploc[is.finite(CHR_NUM) & is.finite(POS)]
  snploc <- snploc[CHR_NUM %between% c(1L, 22L)]  # autosomes by default
  snploc[, .(SNP, CHR_NUM, POS)]
}

# -----------------------------
# Manhattan cumulative position (SAFE: no integer overflow)
# -----------------------------
add_cumpos <- function(dt, chr_col = "CHR_NUM", pos_col = "POS") {
  if (!(chr_col %in% names(dt))) stop("Missing chr column: ", chr_col)
  if (!(pos_col %in% names(dt))) stop("Missing pos column: ", pos_col)

  x <- copy(dt)

  # Force safe types
  x[, CHR_NUM := as.integer(get(chr_col))]
  x[, POS := as.numeric(get(pos_col))]

  x <- x[is.finite(CHR_NUM) & is.finite(POS)]
  x <- x[CHR_NUM %between% c(1L, 22L)]

  setorder(x, CHR_NUM, POS)

  # Per-chromosome "length" in the data (max observed POS)
  chr_sizes <- x[, .(CHR_LEN = max(POS, na.rm = TRUE)), by = CHR_NUM]
  setorder(chr_sizes, CHR_NUM)

  # Numeric offsets (avoid 32-bit integer overflow)
  chr_sizes[, CHR_START := as.numeric(shift(cumsum(as.numeric(CHR_LEN)), fill = 0))]

  x <- merge(
    x,
    chr_sizes[, .(CHR_NUM, CHR_START, CHR_LEN)],
    by = "CHR_NUM",
    all.x = TRUE,
    sort = FALSE
  )

  # Cumulative position (numeric)
  x[, BP_CUM := as.numeric(POS) + as.numeric(CHR_START)]

  # Axis centers: start + (end-start)/2 (overflow-safe)
  axis_df <- chr_sizes[, {
    start <- as.numeric(CHR_START)
    end   <- as.numeric(CHR_START) + as.numeric(CHR_LEN)
    .(CENTER = start + (end - start) / 2)
  }, by = CHR_NUM]

  list(df = x, axis = axis_df)
}

# -----------------------------
# LogP helper
# -----------------------------
add_logp <- function(dt, p_col = NULL, out_col = "LOGP") {
  x <- copy(dt)
  if (is.null(p_col)) p_col <- pick_p_col(x)

  # Compute row-wise; invalid p become NA (quiet coercion)
  pv <- suppressWarnings(as.numeric(as.character(x[[p_col]])))
  ok <- is.finite(pv) & pv > 0 & pv <= 1

  x[[out_col]] <- NA_real_
  x[[out_col]][ok] <- -log10(pv[ok])

  x
}

# -----------------------------
# Safe filenames
# -----------------------------
safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
}