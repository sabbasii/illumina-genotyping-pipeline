#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Rtsne)
})

# ============================================================
# 03_plot_overall_expression_tsne_by_top_pde_snp.R
#
# Purpose
#   For a given PDE gene (e.g. PDE3A), identify its top eQTL SNP from:
#     output/eqtl/results/inspect/pde/eqtl_all_pde.tsv
#   then run t-SNE on curated gene-set expression (samples x genes),
#   and plot samples colored by genotype (0/1/2) for that top SNP.
#
# Grouping logic
#   Genotype groups are kept separate: 0, 1, 2
#
# How to run
#   source scripts/00_config.sh
# Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/03_plot_overall_expression_tsne_by_top_pde_snp.R \
# --gene PDE7A
#
# Optional
#   --gene-set input_data/target_lists/stroke_inflammation_signature_human.txt
#   --perplexity 30
#   --seed 1
#
# Output
#   output/eqtl/results/inspect/pde/overall_expression_by_top_pde_snp/<PDE_GENE>/
#     <PDE_GENE>_<TOP_SNP>_overall_expression_tsne.png
# ============================================================

DEFAULT_GENE_SET <- file.path("input_data", "target_lists", "stroke_inflammation_signature_human.txt")
DEFAULT_PERP <- 30L
DEFAULT_SEED <- 1L

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) stop(paste0("Missing value after ", flag))
  args[w[1] + 1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (!nzchar(REPO_ROOT)) stop("REPO_ROOT not set. Source scripts/00_config.sh first.")

gene_in <- get_arg("--gene", "")
if (!nzchar(gene_in)) stop("Provide --gene (e.g. PDE3A)")
PDE_GENE <- toupper(gene_in)

gene_set_path <- get_arg("--gene-set", DEFAULT_GENE_SET)
if (!grepl("^/", gene_set_path)) gene_set_path <- file.path(REPO_ROOT, gene_set_path)

PERP <- suppressWarnings(as.integer(get_arg("--perplexity", DEFAULT_PERP)))
if (!is.finite(PERP) || PERP <= 1) PERP <- DEFAULT_PERP

SEED <- suppressWarnings(as.integer(get_arg("--seed", DEFAULT_SEED)))
if (!is.finite(SEED)) SEED <- DEFAULT_SEED

# inputs
eqtl_pde_file <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/eqtl_all_pde.tsv")
SNP_file <- file.path(REPO_ROOT, "output/eqtl/SNP_overlap.txt")
GE_file  <- file.path(REPO_ROOT, "output/eqtl/GE_overlap.txt")

for (f in c(eqtl_pde_file, SNP_file, GE_file, gene_set_path)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

# ------------------------------------------------------------
# Identify top SNP for PDE_GENE
# ------------------------------------------------------------
hits <- fread(eqtl_pde_file)
if (!all(c("SNP", "gene") %in% names(hits))) stop("eqtl_all_pde.tsv must have columns: SNP, gene")

p_col <- NULL
for (cand in c("pvalue", "p-value", "p.value", "PValue", "P", "p")) {
  if (cand %in% names(hits)) { p_col <- cand; break }
}
if (is.null(p_col)) stop("No p-value column found in: ", eqtl_pde_file)

hits[, geneU := toupper(gene)]
hits[, P := suppressWarnings(as.numeric(get(p_col)))]
hits <- hits[is.finite(P) & P > 0 & P <= 1]

sub <- hits[geneU == PDE_GENE]
if (nrow(sub) == 0) stop("No PDE hits found for gene: ", PDE_GENE)

setorder(sub, P)
top_snp <- sub$SNP[1]
top_p   <- sub$P[1]

# output dir + filename
out_dir <- file.path(
  REPO_ROOT,
  "output/eqtl/results/inspect/pde",
  "overall_expression_by_top_pde_snp",
  PDE_GENE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, paste0(PDE_GENE, "_", top_snp, "_overall_expression_tsne.png"))

# ------------------------------------------------------------
# Load matrices + curated set
# ------------------------------------------------------------
SNP <- fread(SNP_file)
GE  <- fread(GE_file)

if (!("snpid" %in% names(SNP))) stop("SNP_overlap.txt must have first column named 'snpid'")
if (!("geneid" %in% names(GE)))  stop("GE_overlap.txt must have first column named 'geneid'")

setnames(SNP, "snpid", "SNP")

snp_row <- SNP[SNP == top_snp]
if (nrow(snp_row) == 0) stop("Top SNP not found in SNP_overlap.txt: ", top_snp)

gs <- fread(gene_set_path, header = FALSE, sep = "\n", data.table = FALSE)[, 1]
gs <- unique(toupper(trimws(gs)))
gs <- gs[nzchar(gs)]

GE[, geneU := toupper(geneid)]
GEc <- GE[geneU %in% gs]
if (nrow(GEc) == 0) stop("None of curated genes found in GE_overlap.txt (check symbol matching).")

# sample overlap
snp_samples <- names(SNP)[-1]
ge_samples  <- names(GE)[!(names(GE) %in% c("geneid", "geneU"))]
common_samples <- intersect(snp_samples, ge_samples)
if (length(common_samples) < 10) stop("Too few common sample columns between SNP and GE: ", length(common_samples))

# genotype (keep separate 0/1/2)
g <- suppressWarnings(as.numeric(unlist(snp_row[, ..common_samples])))
df_g <- data.table(sample = common_samples, genotype = g)
df_g <- df_g[is.finite(genotype)]
df_g <- df_g[genotype %in% c(0, 1, 2)]
if (nrow(df_g) < 10) stop("Too few samples with valid genotype after filtering.")

geno_counts <- df_g[, .N, by = genotype][order(genotype)]
if (uniqueN(df_g$genotype) < 2) stop("Need at least 2 genotype groups (0/1/2) to plot.")

keep_samples <- intersect(common_samples, df_g$sample)
if (length(keep_samples) < 10) stop("Too few samples after filtering: ", length(keep_samples))

# curated expression matrix: genes x samples
GEc2 <- GEc[, c("geneid", keep_samples), with = FALSE]
expr_mat <- as.matrix(GEc2[, -"geneid"])
rownames(expr_mat) <- GEc2$geneid
mode(expr_mat) <- "numeric"

# drop genes with all NA or zero variance
ok_gene <- apply(expr_mat, 1, function(x) {
  x <- x[is.finite(x)]
  length(x) >= 3 && stats::sd(x) > 0
})
expr_mat <- expr_mat[ok_gene, , drop = FALSE]
if (nrow(expr_mat) < 3) stop("Too few usable curated genes after filtering (need >=3).")

# ------------------------------------------------------------
# t-SNE on samples x genes (z-scored per gene)
# ------------------------------------------------------------
X <- t(expr_mat)          # samples x genes
Xs <- scale(X)            # z-score per gene across samples
Xs <- Xs[complete.cases(Xs), , drop = FALSE]
if (nrow(Xs) < 10) stop("Too few complete samples for t-SNE after scaling/filtering.")

# align genotype to Xs rows
df_g2 <- df_g[match(rownames(Xs), sample)]
if (any(is.na(df_g2$genotype))) stop("Internal error: sample alignment failed.")

# perplexity must be < (n-1)/3
max_perp <- floor((nrow(Xs) - 1) / 3)
if (max_perp < 2) stop("Not enough samples for t-SNE (need more samples).")
if (PERP > max_perp) {
  message("Perplexity too large for n=", nrow(Xs), "; using perplexity=", max_perp)
  PERP <- max_perp
}

set.seed(SEED)
ts <- Rtsne(
  Xs,
  dims = 2,
  perplexity = PERP,
  check_duplicates = FALSE,
  pca = TRUE,
  verbose = FALSE
)

pdat <- data.table(
  sample = rownames(Xs),
  tSNE1 = ts$Y[, 1],
  tSNE2 = ts$Y[, 2],
  genotype = factor(df_g2$genotype, levels = c(0, 1, 2))
)

# ------------------------------------------------------------
# Plot (colored by genotype 0/1/2)
# ------------------------------------------------------------
p <- ggplot(pdat, aes(x = tSNE1, y = tSNE2, color = genotype)) +
  geom_point(size = 2.4, alpha = 0.9) +
  theme_classic(base_size = 13) +
  labs(
    title = paste0(PDE_GENE, " ", top_snp, " overall expression t-SNE"),
    subtitle = paste0(
      "Curated genes used: ", nrow(expr_mat),
      " | Genotype counts: ", paste0(geno_counts$genotype, "=", geno_counts$N, collapse = ", "),
      " | Top SNP p=", format(top_p, scientific = TRUE, digits = 3),
      " | perplexity=", PERP, ", seed=", SEED
    ),
    color = "Genotype (0/1/2)",
    x = "t-SNE1",
    y = "t-SNE2"
  )

ggsave(out_png, p, width = 7.2, height = 5.2, dpi = 300)

cat("PDE gene: ", PDE_GENE, "\n", sep = "")
cat("Top SNP:  ", top_snp, " (p=", format(top_p, scientific = TRUE, digits = 3), ")\n", sep = "")
cat("Curated genes used (usable after QC): ", nrow(expr_mat), "\n", sep = "")
cat("Genotype counts: ", paste0(geno_counts$genotype, "=", geno_counts$N, collapse = ", "), "\n", sep = "")
cat("Wrote:\n  ", out_png, "\n", sep = "")