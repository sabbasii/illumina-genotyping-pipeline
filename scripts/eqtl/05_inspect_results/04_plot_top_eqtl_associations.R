#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 04_plot_top_eqtl_associations.R
#
# Genotype → expression association plots for top eQTL hits.
#
# For each selected SNP–gene pair, this script pulls:
#   - genotype values (0/1/2 ALT allele counts) from SNP_overlap.txt
#   - expression values from GE_overlap.txt
# then plots expression by genotype group (boxplot + jitter).
#
# Modes:
#   --mode cis : selects top cis pairs from results/eqtl_cis.tsv
#   --mode all : selects top pairs from results/eqtl_all.tsv
#
# Inputs (under REPO_ROOT/output/eqtl/):
#   - results/eqtl_cis.tsv or results/eqtl_all.tsv
#   - SNP_overlap.txt  (first col: snpid; remaining cols: sample IDs)
#   - GE_overlap.txt   (first col: geneid; remaining cols: sample IDs)
#
# Outputs:
#   - REPO_ROOT/output/eqtl/results/inspect/associations/cis/
#   - REPO_ROOT/output/eqtl/results/inspect/associations/all/
#   - eqtl_<mode>_top_hits.tsv (in the corresponding mode folder)
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

TOP_N <- suppressWarnings(as.integer(get_arg("--top", "20")))
if (!is.finite(TOP_N) || TOP_N <= 0) TOP_N <- 20L

UNIQUE_BY_GENE <- tolower(get_arg("--unique-by-gene", if (mode == "cis") "true" else "false"))
UNIQUE_BY_GENE <- UNIQUE_BY_GENE %in% c("true", "t", "1", "yes", "y")

RANK_BY_FDR <- tolower(get_arg("--rank-by-fdr", "false")) %in% c("true", "t", "1", "yes", "y")

FDR_FILTER_ON <- tolower(get_arg("--fdr-filter", "false")) %in% c("true", "t", "1", "yes", "y")
FDR_CUTOFF <- suppressWarnings(as.numeric(get_arg("--fdr", "0.05")))

P_FILTER_ON <- tolower(get_arg("--p-filter", "false")) %in% c("true", "t", "1", "yes", "y")
P_CUTOFF <- suppressWarnings(as.numeric(get_arg("--p", "1e-5")))

# ---- IO ----
res_file <- if (mode == "cis") paths$eqtl_cis else paths$eqtl_all
assert_file(res_file)

snp_file <- paths$SNP_overlap
ge_file  <- paths$GE_overlap
assert_file(snp_file)
assert_file(ge_file)

assoc_dir <- file.path(paths$inspect_dir, "associations", mode)
dir.create(assoc_dir, recursive = TRUE, showWarnings = FALSE)

hits_out <- file.path(assoc_dir, paste0("eqtl_", mode, "_top_hits.tsv"))

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

# beta column optional (for subtitle)
beta_col <- NULL
beta_candidates <- c("beta", "Beta", "BETA", "slope", "Slope", "effect", "Effect")
if (any(beta_candidates %in% names(x))) beta_col <- pick_beta_col(x)

# standard numeric cols
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

# optional filters
if (FDR_FILTER_ON && !all(is.na(x$FDRv))) x <- x[FDRv <= FDR_CUTOFF]
if (P_FILTER_ON) x <- x[P <= P_CUTOFF]

if (nrow(x) == 0) stop("No hits after filtering. Try relaxing thresholds or disabling filters.")

# ranking
if (RANK_BY_FDR && !all(is.na(x$FDRv))) {
  setorder(x, FDRv, P)
} else {
  setorder(x, P, FDRv)
}

# optionally keep only best hit per gene
if (UNIQUE_BY_GENE) {
  x <- x[, .SD[1], by = gene]
  if (RANK_BY_FDR && !all(is.na(x$FDRv))) setorder(x, FDRv, P) else setorder(x, P, FDRv)
}

hits <- x[1:min(TOP_N, .N), .(SNP, gene, P, FDRv, BETAv)]
fwrite(hits, hits_out, sep = "\t")

# ---- load matrices ----
SNP <- fread(snp_file)
GE  <- fread(ge_file)

if (!("snpid" %in% names(SNP))) stop("SNP_overlap.txt must have first column named 'snpid'")
if (!("geneid" %in% names(GE)))  stop("GE_overlap.txt must have first column named 'geneid'")

setnames(SNP, "snpid", "SNP")

snp_samples <- names(SNP)[-1]
ge_samples  <- names(GE)[-1]
common_samples <- intersect(snp_samples, ge_samples)
if (length(common_samples) < 5) stop("Too few overlapping sample columns: ", length(common_samples))

# ---- plotting ----
for (i in seq_len(nrow(hits))) {
  snp_id  <- hits$SNP[i]
  gene_id <- hits$gene[i]
  pval    <- hits$P[i]
  fdr     <- hits$FDRv[i]
  beta    <- hits$BETAv[i]

  snp_row <- SNP[SNP == snp_id]
  ge_row  <- GE[geneid == gene_id]

  if (nrow(snp_row) == 0) { message("Skipping (SNP not found): ", snp_id); next }
  if (nrow(ge_row)  == 0) { message("Skipping (gene not found): ", gene_id); next }

  g <- suppressWarnings(as.numeric(unlist(snp_row[, ..common_samples])))
  e <- suppressWarnings(as.numeric(unlist(ge_row[, ..common_samples])))

  df <- data.frame(sample = common_samples, genotype = g, expression = e)
  df <- df[is.finite(df$genotype) & is.finite(df$expression), ]

  df$genotype_f <- factor(df$genotype, levels = sort(unique(df$genotype)))

  if (length(unique(df$genotype_f)) < 2) {
    message("Skipping (only one genotype group): ", snp_id, " -> ", gene_id)
    next
  }

  subtitle_parts <- c(
    paste0("p=", format(pval, scientific = TRUE, digits = 3)),
    if (!is.na(fdr)) paste0("FDR=", format(fdr, scientific = TRUE, digits = 3)) else NULL,
    if (!is.na(beta)) paste0("beta=", signif(beta, 3)) else NULL,
    paste0("n=", nrow(df))
  )

  p <- ggplot(df, aes(x = genotype_f, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, size = 1.4, alpha = 0.75) +
    labs(
      title = paste0("eQTL association: ", snp_id, " \u2192 ", gene_id),
      subtitle = paste(subtitle_parts, collapse = " | "),
      x = "Genotype (0/1/2)",
      y = "Expression"
    ) +
    theme_classic(base_size = 13)

  out_png <- file.path(
    assoc_dir,
    paste0("eqtl_", mode, "_assoc_", sprintf("%02d", i), "_",
           safe_name(snp_id), "__", safe_name(gene_id), ".png")
  )

  ggsave(out_png, p, width = 7.5, height = 5.0, dpi = 300)
  message("Wrote: ", out_png)
}

message("Done.")
message("Hits table: ", hits_out)

# ---- Run ----
# source scripts/00_config.sh
#
# Cis:
# Rscript scripts/eqtl/05_inspect_results/04_plot_top_eqtl_associations.R --mode cis --top 10
#
# All:
# Rscript scripts/eqtl/05_inspect_results/04_plot_top_eqtl_associations.R --mode all --top 10
#
# Optional:
#   --unique-by-gene true|false
#   --rank-by-fdr true
#   --p-filter true  --p 1e-6
#   --fdr-filter true --fdr 0.05
