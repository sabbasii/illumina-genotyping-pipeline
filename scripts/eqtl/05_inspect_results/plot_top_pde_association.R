#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# plot_top_pde_association.R
#
# For each PDE gene in eqtl_all_pde.tsv, select the top SNP hit
# (smallest p-value) and plot expression vs genotype (0/1/2):
#   boxplot + jitter
#
# Inputs:
#   output/eqtl/results/inspect/pde/eqtl_all_pde.tsv
#   output/eqtl/SNP_overlap.txt
#   output/eqtl/GE_overlap.txt
#
# Outputs:
#   output/eqtl/results/inspect/pde/associations/
#     pde_assoc__<GENE>__<SNP>.png
#
# Run:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/05_inspect_results/plot_top_pde_association.R
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT not set. Did you source scripts/00_config.sh?")

in_hits <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/eqtl_all_pde.tsv")
snp_file <- file.path(REPO_ROOT, "output/eqtl/SNP_overlap.txt")
ge_file  <- file.path(REPO_ROOT, "output/eqtl/GE_overlap.txt")
out_dir  <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/associations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(in_hits, snp_file, ge_file)) if (!file.exists(f)) stop("Missing file: ", f)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# ---- load PDE hits ----
hits <- fread(in_hits)
if (!all(c("SNP", "gene") %in% names(hits))) stop("eqtl_all_pde.tsv must have columns: SNP, gene")

# detect p-value column
p_col <- NULL
for (cand in c("pvalue","p-value","p.value","PValue","P","p")) {
  if (cand %in% names(hits)) { p_col <- cand; break }
}
if (is.null(p_col)) stop("Could not find a p-value column in eqtl_all_pde.tsv")

hits[, P := suppressWarnings(as.numeric(get(p_col)))]
hits <- hits[is.finite(P) & P > 0 & P <= 1]
if (nrow(hits) == 0) stop("No valid p-values in eqtl_all_pde.tsv")

# pick top hit per gene
setorder(hits, gene, P)
top <- hits[, .SD[1], by = gene]
setorder(top, P)

cat("Genes with at least one hit:", nrow(top), "\n")

# ---- load matrices ----
SNP <- fread(snp_file)
GE  <- fread(ge_file)

if (!("snpid" %in% names(SNP))) stop("SNP_overlap.txt must have first column 'snpid'")
if (!("geneid" %in% names(GE))) stop("GE_overlap.txt must have first column 'geneid'")

setnames(SNP, "snpid", "SNP")

snp_samples <- names(SNP)[-1]
ge_samples  <- names(GE)[-1]
common_samples <- intersect(snp_samples, ge_samples)
if (length(common_samples) < 5) stop("Too few overlapping sample columns: ", length(common_samples))

# ---- plot loop ----
written <- 0L
skipped <- 0L

for (i in seq_len(nrow(top))) {
  snp_id  <- top$SNP[i]
  gene_id <- top$gene[i]
  pval    <- top$P[i]

  snp_row <- SNP[SNP == snp_id]
  ge_row  <- GE[geneid == gene_id]

  if (nrow(snp_row) == 0 || nrow(ge_row) == 0) {
    skipped <- skipped + 1L
    next
  }

  g <- suppressWarnings(as.numeric(unlist(snp_row[, ..common_samples])))
  e <- suppressWarnings(as.numeric(unlist(ge_row[, ..common_samples])))

  df <- data.frame(sample = common_samples, genotype = g, expression = e)
  df <- df[is.finite(df$genotype) & is.finite(df$expression), ]
  if (nrow(df) < 8) { skipped <- skipped + 1L; next }

  df$genotype_f <- factor(df$genotype, levels = c(0,1,2))
  df$genotype_f <- droplevels(df$genotype_f)
  if (nlevels(df$genotype_f) < 2) { skipped <- skipped + 1L; next }

  p <- ggplot(df, aes(x = genotype_f, y = expression)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, size = 1.3, alpha = 0.75) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0(gene_id, " ~ ", snp_id),
      subtitle = paste0("Top eQTL per gene (from eqtl_all_pde.tsv), p=",
                        format(pval, scientific = TRUE, digits = 3),
                        " | n=", nrow(df)),
      x = "Genotype (ALT count: 0/1/2)",
      y = "Expression"
    )

  out_png <- file.path(out_dir, paste0("pde_assoc__", safe_name(gene_id), "__", safe_name(snp_id), ".png"))
  ggsave(out_png, p, width = 7.2, height = 5.0, dpi = 300)
  written <- written + 1L
}

cat("Done.\n")
cat("Plots written: ", written, "\n", sep = "")
cat("Skipped:       ", skipped, "\n", sep = "")
cat("Output dir:    ", out_dir, "\n", sep = "")