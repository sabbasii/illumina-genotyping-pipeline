#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 05_target_gene_eqtl_summary.R
#
# Targeted eQTL inspection for a gene set.
#
# Filters eQTL association results to a target gene set defined by:
#   - a regular expression (--gene-regex), and/or
#   - a plain-text gene list file (--gene-list)
#
# Works with either cis or all results:
#   --mode cis  -> results/eqtl_cis.tsv
#   --mode all  -> results/eqtl_all.tsv
#
# Outputs are written under:
#   REPO_ROOT/output/eqtl/results/inspect/targets/<target_name>/<mode>/
#
# Files:
#   - eqtl_<mode>_target_hits.tsv         (filtered hits, sorted)
#   - eqtl_<mode>_target_gene_summary.tsv (per-gene counts/best p)
#   - eqtl_<mode>_target_snp_summary.tsv  (per-SNP counts/best p)
#   - eqtl_<mode>_target_summary.txt      (quick readable summary)
#   - eqtl_<mode>_target_gene_counts.png  (barplot: top genes by hit count)
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

target_name <- get_arg("--target-name", "target")
if (!nzchar(target_name)) target_name <- "target"

gene_regex <- get_arg("--gene-regex", "")
gene_list_file <- get_arg("--gene-list", "")

TOP_GENE_PLOT <- suppressWarnings(as.integer(get_arg("--top-genes", "25")))
if (!is.finite(TOP_GENE_PLOT) || TOP_GENE_PLOT <= 0) TOP_GENE_PLOT <- 25L

RANK_BY_FDR <- tolower(get_arg("--rank-by-fdr", "false")) %in% c("true", "t", "1", "yes", "y")

FDR_FILTER_ON <- tolower(get_arg("--fdr-filter", "false")) %in% c("true", "t", "1", "yes", "y")
FDR_CUTOFF <- suppressWarnings(as.numeric(get_arg("--fdr", "0.05")))

P_FILTER_ON <- tolower(get_arg("--p-filter", "false")) %in% c("true", "t", "1", "yes", "y")
P_CUTOFF <- suppressWarnings(as.numeric(get_arg("--p", "1e-5")))

# ---- IO ----
res_file <- if (mode == "cis") paths$eqtl_cis else paths$eqtl_all
assert_file(res_file)

out_dir <- file.path(paths$inspect_dir, "targets", safe_name(target_name), mode)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_hits   <- file.path(out_dir, paste0("eqtl_", mode, "_target_hits.tsv"))
out_genes  <- file.path(out_dir, paste0("eqtl_", mode, "_target_gene_summary.tsv"))
out_snps   <- file.path(out_dir, paste0("eqtl_", mode, "_target_snp_summary.tsv"))
out_txt    <- file.path(out_dir, paste0("eqtl_", mode, "_target_summary.txt"))
out_png    <- file.path(out_dir, paste0("eqtl_", mode, "_target_gene_counts.png"))

# ---- load association results ----
x <- fread(res_file)

need_min <- c("SNP", "gene")
miss_min <- setdiff(need_min, names(x))
if (length(miss_min) > 0) stop("Association table missing columns: ", paste(miss_min, collapse = ", "))

p_col <- pick_p_col(x)

fdr_col <- NULL
if (any(c("FDR", "fdr", "qvalue", "q.value", "q-value") %in% names(x))) {
  fdr_col <- pick_fdr_col(x)
}

beta_col <- NULL
beta_candidates <- c("beta", "Beta", "BETA", "slope", "Slope", "effect", "Effect")
if (any(beta_candidates %in% names(x))) beta_col <- pick_beta_col(x)

# numeric p/FDR/beta
x[, P := suppressWarnings(as.numeric(get(p_col)))]
x <- x[is.finite(P) & P > 0 & P <= 1]

if (!is.null(fdr_col)) {
  x[, FDRv := suppressWarnings(as.numeric(get(fdr_col)))]
} else {
  x[, FDRv := NA_real_]
}

if (!is.null(beta_col)) {
  x[, BETAv := suppressWarnings(as.numeric(get(beta_col)))]
} else {
  x[, BETAv := NA_real_]
}

# ---- gene filtering (regex and/or list file) ----
genes_keep <- NULL

if (nzchar(gene_list_file)) {
  if (!file.exists(gene_list_file)) stop("Missing gene list: ", gene_list_file)
  gl <- fread(gene_list_file, header = FALSE, sep = "\n", data.table = FALSE)[, 1]
  gl <- trimws(gl)
  gl <- gl[nzchar(gl)]
  genes_keep <- unique(gl)
}

if (nzchar(gene_regex)) {
  if (is.null(genes_keep)) {
    x <- x[grepl(gene_regex, gene)]
  } else {
    x <- x[gene %in% genes_keep & grepl(gene_regex, gene)]
  }
} else if (!is.null(genes_keep)) {
  x <- x[gene %in% genes_keep]
} else {
  stop("No gene filter provided. Use --gene-regex and/or --gene-list.")
}

if (nrow(x) == 0) stop("No hits left after gene filtering. Check your regex/list and gene naming.")

# ---- optional significance filters ----
if (FDR_FILTER_ON && !all(is.na(x$FDRv))) x <- x[FDRv <= FDR_CUTOFF]
if (P_FILTER_ON) x <- x[P <= P_CUTOFF]

if (nrow(x) == 0) stop("No hits left after applying significance filters. Relax thresholds or disable filters.")

# ---- rank + write hits ----
if (RANK_BY_FDR && !all(is.na(x$FDRv))) setorder(x, FDRv, P) else setorder(x, P, FDRv)

hits <- x[, .(SNP, gene, P, FDRv, BETAv)]
fwrite(hits, out_hits, sep = "\t")

# ---- summaries ----
gene_sum <- hits[, .(
  n_pairs = .N,
  unique_snps = uniqueN(SNP),
  best_p = min(P, na.rm = TRUE),
  best_fdr = if (!all(is.na(FDRv))) min(FDRv, na.rm = TRUE) else NA_real_
), by = gene][order(-n_pairs, best_p)]

snp_sum <- hits[, .(
  n_pairs = .N,
  unique_genes = uniqueN(gene),
  best_p = min(P, na.rm = TRUE),
  best_fdr = if (!all(is.na(FDRv))) min(FDRv, na.rm = TRUE) else NA_real_
), by = SNP][order(-n_pairs, best_p)]

fwrite(gene_sum, out_genes, sep = "\t")
fwrite(snp_sum, out_snps, sep = "\t")

# ---- plot (top genes by hit count) ----
plot_dt <- gene_sum[1:min(TOP_GENE_PLOT, .N)]
if (nrow(plot_dt) > 0) {
  p <- ggplot(plot_dt, aes(x = reorder(gene, n_pairs), y = n_pairs)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Target eQTL burden (", mode, ")"),
      subtitle = paste0("Target: ", target_name, " | genes shown: ", nrow(plot_dt)),
      x = "",
      y = "Number of SNP–gene pairs"
    )
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
} else {
  out_png <- "(no plot written)"
}

# ---- summary text ----
lines <- c(
  "Target eQTL summary",
  "-------------------",
  paste0("Input:  ", res_file),
  paste0("Mode:   ", mode),
  paste0("Target: ", target_name),
  paste0("Gene regex: ", if (nzchar(gene_regex)) gene_regex else "(none)"),
  paste0("Gene list:  ", if (nzchar(gene_list_file)) gene_list_file else "(none)"),
  "",
  paste0("Rows kept:        ", nrow(hits)),
  paste0("Unique genes:     ", uniqueN(hits$gene)),
  paste0("Unique SNPs:      ", uniqueN(hits$SNP)),
  "",
  paste0("Filters: p-filter=", ifelse(P_FILTER_ON, "ON", "OFF"),
         " (p <= ", format(P_CUTOFF, scientific = TRUE), "), ",
         "FDR-filter=", ifelse(FDR_FILTER_ON, "ON", "OFF"),
         " (FDR <= ", format(FDR_CUTOFF, scientific = TRUE), ")"),
  "",
  "Outputs:",
  paste0("  Hits:        ", out_hits),
  paste0("  Gene summary:", out_genes),
  paste0("  SNP summary: ", out_snps),
  paste0("  Plot:        ", out_png)
)

writeLines(lines, out_txt)

cat("Wrote:\n  ", out_txt, "\n  ", out_hits, "\n  ", out_genes, "\n  ", out_snps, "\n", sep = "")
if (out_png != "(no plot written)") cat("  ", out_png, "\n", sep = "")

# ---- Run ----
# source scripts/00_config.sh
#
# PDE via regex:
# Rscript scripts/eqtl/05_inspect_results/05_target_gene_eqtl_summary.R \
#   --mode all --target-name pde --gene-regex "^PDE[0-9]"
#
# PDE via gene list:
# Rscript scripts/eqtl/05_inspect_results/05_target_gene_eqtl_summary.R \
#   --mode cis --target-name pde --gene-list input_data/target_lists/PDE_genes.txt
#
# Both (intersection):
# Rscript scripts/eqtl/05_inspect_results/05_target_gene_eqtl_summary.R \
#   --mode cis --target-name pde \
#   --gene-regex "^PDE[0-9]" --gene-list input_data/target_lists/PDE_genes.txt
#
# Optional significance filters:
#   --p-filter true --p 1e-6
#   --fdr-filter true --fdr 0.05
