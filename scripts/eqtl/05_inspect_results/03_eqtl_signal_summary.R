#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# 03_eqtl_signal_summary.R
#
# Summarize eQTL signal before per-hit genotype→expression plots.
# Produces a compact overview of:
#   - number of tested pairs (or SNPs if collapsed)
#   - counts passing p-value / FDR thresholds
#   - top hits table
#   - top genes / top SNPs by hit burden (pair-level only)
#   - hits per chromosome (requires snpsloc.txt)
#
# Modes:
#   --mode cis : uses results/eqtl_cis.tsv (default; no collapsing)
#   --mode all : uses results/eqtl_all.tsv (default: collapse to min p per SNP)
#
# Output naming is mode-first so cis outputs sort together and all outputs sort together:
#   eqtl_cis_*   and   eqtl_all_*
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

source(file.path(REPO_ROOT, "scripts/eqtl/utils/inspect_helpers.R"))
paths <- get_eqtl_paths()

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) return(default)
  args[w[1] + 1]
}

mode <- tolower(get_arg("--mode", "cis"))
if (!(mode %in% c("cis", "all"))) stop("Invalid --mode. Use: cis or all")

collapse_arg <- tolower(get_arg("--collapse", "auto"))
collapse <- FALSE
if (collapse_arg == "auto") {
  collapse <- (mode == "all")
} else if (collapse_arg %in% c("true", "t", "1", "yes", "y")) {
  collapse <- TRUE
} else if (collapse_arg %in% c("false", "f", "0", "no", "n")) {
  collapse <- FALSE
} else {
  stop("Invalid --collapse. Use: auto | TRUE | FALSE")
}

TOP_N <- suppressWarnings(as.integer(get_arg("--top", "100")))
if (!is.finite(TOP_N) || TOP_N <= 0) TOP_N <- 100L

P_CUTOFF <- suppressWarnings(as.numeric(get_arg("--p", "1e-5")))
FDR_CUTOFF <- suppressWarnings(as.numeric(get_arg("--fdr", "0.05")))

# ---- IO ----
res_file <- if (mode == "cis") paths$eqtl_cis else paths$eqtl_all
assert_file(res_file)

snploc <- NULL
if (file.exists(paths$snpsloc)) snploc <- read_snpsloc(paths$snpsloc) # SNP, CHR_NUM, POS

summary_dir <- file.path(paths$inspect_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- paste0("eqtl_", mode)
suffix <- if (collapse) "_minp_per_snp" else ""

out_txt  <- file.path(summary_dir, paste0(prefix, "_signal_summary", suffix, ".txt"))
out_hits <- file.path(summary_dir, paste0(prefix, "_top_hits", suffix, ".tsv"))
out_chr  <- file.path(summary_dir, paste0(prefix, "_hits_by_chr", suffix, ".tsv"))

# Pair-level only (do not add suffix)
out_genes <- file.path(summary_dir, paste0(prefix, "_top_genes.tsv"))
out_snps  <- file.path(summary_dir, paste0(prefix, "_top_snps.tsv"))

# ---- load + standardize ----
x <- fread(res_file)

if (!("SNP" %in% names(x)))  stop("Association table missing column: SNP")
if (!("gene" %in% names(x))) stop("Association table missing column: gene")

p_col <- pick_p_col(x)
fdr_col <- NULL
if (any(c("FDR", "fdr", "qvalue", "q.value", "q-value") %in% names(x))) {
  fdr_col <- pick_fdr_col(x)
}

x[, P := suppressWarnings(as.numeric(get(p_col)))]
x <- x[is.finite(P) & P > 0 & P <= 1]

if (!is.null(fdr_col)) {
  x[, FDRv := suppressWarnings(as.numeric(get(fdr_col)))]
} else {
  x[, FDRv := NA_real_]
}

# ---- collapse (optional) ----
# Collapse = one row per SNP using min p-value across genes (and best/lowest FDR if available).
if (collapse) {
  if (!is.null(fdr_col)) {
    y <- x[, .(
      P = min(P, na.rm = TRUE),
      FDRv = min(FDRv, na.rm = TRUE)
    ), by = SNP]
  } else {
    y <- x[, .(P = min(P, na.rm = TRUE)), by = SNP]
    y[, FDRv := NA_real_]
  }
  y[, gene := NA_character_]
  x_use <- y
} else {
  x_use <- x[, .(SNP, gene, P, FDRv)]
}

x_use[, LOGP := -log10(P)]

# ---- counts ----
n_rows <- nrow(x_use)
n_snps <- length(unique(x_use$SNP))
n_genes <- if (!collapse) length(unique(x_use$gene)) else NA_integer_

n_p <- sum(x_use$P <= P_CUTOFF, na.rm = TRUE)
n_fdr <- if (!all(is.na(x_use$FDRv))) sum(x_use$FDRv <= FDR_CUTOFF, na.rm = TRUE) else NA_integer_

# ---- top hits ----
setorder(x_use, P)
top_hits <- x_use[1:min(TOP_N, .N)]
fwrite(top_hits, out_hits, sep = "\t")

# ---- hits by chromosome (if snpsloc available) ----
if (!is.null(snploc)) {
  chr_dt <- merge(x_use[, .(SNP, P)], snploc[, .(SNP, CHR_NUM)], by = "SNP", all.x = TRUE)
  chr_dt <- chr_dt[!is.na(CHR_NUM) & CHR_NUM %between% c(1L, 22L)]
  hits_by_chr <- chr_dt[, .(
    n = .N,
    n_p = sum(P <= P_CUTOFF, na.rm = TRUE),
    best_p = min(P, na.rm = TRUE)
  ), by = CHR_NUM][order(CHR_NUM)]
  fwrite(hits_by_chr, out_chr, sep = "\t")
} else {
  hits_by_chr <- NULL
}

# ---- top genes / top SNPs by burden (pair-level only) ----
if (!collapse) {
  top_genes <- x_use[, .(
    n_pairs = .N,
    unique_snps = uniqueN(SNP),
    best_p = min(P, na.rm = TRUE),
    best_fdr = if (!all(is.na(FDRv))) min(FDRv, na.rm = TRUE) else NA_real_
  ), by = gene][order(-n_pairs, best_p)]
  fwrite(top_genes, out_genes, sep = "\t")

  top_snps <- x_use[, .(
    n_pairs = .N,
    unique_genes = uniqueN(gene),
    best_p = min(P, na.rm = TRUE),
    best_fdr = if (!all(is.na(FDRv))) min(FDRv, na.rm = TRUE) else NA_real_
  ), by = SNP][order(-n_pairs, best_p)]
  fwrite(top_snps, out_snps, sep = "\t")
}

# ---- write summary text ----
lines <- c(
  "eQTL signal summary",
  "-------------------",
  paste0("Input: ", res_file),
  paste0("Mode:  ", mode),
  paste0("Collapse to min p per SNP: ", ifelse(collapse, "TRUE", "FALSE")),
  "",
  paste0("Rows analyzed: ", n_rows),
  paste0("Unique SNPs:   ", n_snps),
  if (!collapse) paste0("Unique genes:  ", n_genes) else NULL,
  "",
  paste0("Thresholds: p <= ", format(P_CUTOFF, scientific = TRUE),
         " ; FDR <= ", format(FDR_CUTOFF, scientific = TRUE)),
  paste0("Count p-pass:   ", n_p),
  if (!is.na(n_fdr)) paste0("Count FDR-pass: ", n_fdr) else "Count FDR-pass: NA (no FDR column)",
  "",
  paste0("Top hits:       ", out_hits),
  if (!is.null(snploc)) paste0("Hits by chr:    ", out_chr) else "Hits by chr:    (snpsloc.txt missing; skipped)",
  if (!collapse) paste0("Top genes:      ", out_genes) else "Top genes:      (collapsed; skipped)",
  if (!collapse) paste0("Top SNPs:       ", out_snps) else "Top SNPs:       (collapsed; skipped)"
)

writeLines(lines, out_txt)

cat("Wrote:\n  ", out_txt, "\n  ", out_hits, "\n", sep = "")
if (!is.null(snploc)) cat("  ", out_chr, "\n", sep = "")
if (!collapse) cat("  ", out_genes, "\n  ", out_snps, "\n", sep = "")

# ---- Run ----
# source scripts/00_config.sh
# Rscript scripts/eqtl/05_inspect_results/03_eqtl_signal_summary.R --mode cis
# Rscript scripts/eqtl/05_inspect_results/03_eqtl_signal_summary.R --mode all
#
# Optional:
#   --collapse TRUE|FALSE|auto   (default: auto; collapses only for --mode all)
#   --top 100                   (default: 100)
#   --p 1e-5                    (default: 1e-5)
#   --fdr 0.05                  (default: 0.05)
