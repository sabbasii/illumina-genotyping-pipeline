#!/usr/bin/env Rscript
# scripts/eqtl/utils/inspect_helpers.R
# Shared helpers for scripts/eqtl/05_inspect_results/*.R

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------
get_eqtl_paths <- function(run = Sys.getenv("RUN", "genotype_run1")) {
  REPO_ROOT <- Sys.getenv("REPO_ROOT")
  if (is.null(REPO_ROOT) || REPO_ROOT == "") {
    stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")
  }

  base_eqtl   <- file.path(REPO_ROOT, "output", run, "eqtl")
  results_dir <- file.path(base_eqtl, "results")
  inspect_dir <- file.path(results_dir, "inspect")

  qc_dir        <- file.path(inspect_dir, "qc")
  manhattan_dir <- file.path(inspect_dir, "manhattan")
  effects_dir   <- file.path(inspect_dir, "effects")
  pde_dir       <- file.path(inspect_dir, "pde")

  # Ensure dirs exist
  dirs <- c(results_dir, inspect_dir, qc_dir, manhattan_dir, effects_dir, pde_dir)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  list(
    REPO_ROOT = REPO_ROOT,
    RUN = run,
    base_eqtl = base_eqtl,
    results_dir = results_dir,
    inspect_dir = inspect_dir,
    qc_dir = qc_dir,
    manhattan_dir = manhattan_dir,
    effects_dir = effects_dir,
    pde_dir = pde_dir,

    # Common inputs
    eqtl_all  = file.path(results_dir, "eqtl_all.tsv"),
    eqtl_cis  = file.path(results_dir, "eqtl_cis.tsv"),
    eqtl_trans= file.path(results_dir, "eqtl_trans.tsv"),
    snpsloc   = file.path(base_eqtl, "snpsloc.txt"),
    geneloc   = file.path(base_eqtl, "geneloc.txt"),
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

  # Expect 3 cols: snpid, chr, pos (but be tolerant)
  if (ncol(snploc) < 3) stop("snpsloc must have at least 3 columns (snpid, chr, pos).")

  # If already has SNP/CHR/POS, keep; otherwise rename first 3.
  if (!all(c("SNP", "CHR", "POS") %in% names(snploc))) {
    setnames(snploc, names(snploc)[1:3], c("SNP", "CHR", "POS"))
  }

  snploc[, CHR := as.character(CHR)]
  snploc[, CHR := gsub("^chr", "", CHR, ignore.case = TRUE)]
  snploc[, CHR_NUM := suppressWarnings(as.integer(CHR))]
  snploc[, POS := suppressWarnings(as.integer(POS))]

  snploc <- snploc[is.finite(CHR_NUM) & is.finite(POS)]
  snploc <- snploc[CHR_NUM %between% c(1L, 22L)]  # autosomes by default
  snploc[, .(SNP, CHR_NUM, POS)]
}

# -----------------------------
# Manhattan cumulative position
# -----------------------------
add_cumpos <- function(dt, chr_col = "CHR_NUM", pos_col = "POS") {
  if (!(chr_col %in% names(dt))) stop("Missing chr column: ", chr_col)
  if (!(pos_col %in% names(dt))) stop("Missing pos column: ", pos_col)

  x <- copy(dt)
  x[, CHR_NUM := as.integer(get(chr_col))]
  x[, POS := as.integer(get(pos_col))]
  x <- x[is.finite(CHR_NUM) & is.finite(POS)]
  x <- x[CHR_NUM %between% c(1L, 22L)]

  setorder(x, CHR_NUM, POS)

  chr_sizes <- x[, .(CHR_LEN = max(POS, na.rm = TRUE)), by = CHR_NUM]
  setorder(chr_sizes, CHR_NUM)
  chr_sizes[, CHR_START := shift(cumsum(CHR_LEN), fill = 0L)]

  x <- merge(x, chr_sizes[, .(CHR_NUM, CHR_START)], by = "CHR_NUM", all.x = TRUE, sort = FALSE)
  x[, BP_CUM := POS + CHR_START]

  axis_df <- x[, .(CENTER = (min(BP_CUM) + max(BP_CUM)) / 2), by = CHR_NUM]

  list(df = x, axis = axis_df)
}

# -----------------------------
# LogP helper
# -----------------------------
add_logp <- function(dt, p_col = NULL, out_col = "LOGP") {
  x <- copy(dt)
  if (is.null(p_col)) p_col <- pick_p_col(x)

  p <- suppressWarnings(as.numeric(x[[p_col]]))
  p <- p[is.finite(p) & p > 0 & p <= 1]

  if (length(p) == 0) {
    x[[out_col]] <- NA_real_
    return(x)
  }

  # Write back row-wise (keep NA for invalid)
  x[[out_col]] <- suppressWarnings(-log10(as.numeric(x[[p_col]])))
  x
}

# -----------------------------
# Safe filenames
# -----------------------------
safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
}
