#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ============================================================
# 04_plot_volcano_overall_expression_by_top_pde_snp.R
#
# Purpose
#   Volcano plot from:
#     <PDE_GENE>_<TOP_SNP>_overall_expression_limma.tsv
#
# Comparison
#   Group2 (ALT/ALT) vs Group1 (ALT/REF + REF/REF)
#   i.e., genotype 2 vs (genotype 1 + genotype 0)
#
# Coloring
#   - Significant Up   (logFC > 0): red
#   - Significant Down (logFC < 0): blue
#   - Not significant  : grey
#
# Thresholds (defaults)
#   --p  0.05   (uses P.Value by default; set --use-fdr to use adj.P.Val)
#   --fc 0.50   (absolute logFC cutoff)
#
# Labels
#   Labels top N most significant genes (default 10) with leader lines.
#
# How to run
#   source scripts/00_config.sh
# Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/04_plot_volcano_overall_expression_by_top_pde_snp.R \
#   --file output/eqtl/results/inspect/pde/overall_expression_by_top_pde_snp/PDE3A/PDE3A_rs6433690_overall_expression_limma.tsv

# Output
#   Saves next to input:
#     <PDE_GENE>_<TOP_SNP>_overall_expression_volcano.png
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) stop(paste0("Missing value after ", flag))
  args[w[1] + 1]
}

infile <- get_arg("--file", "")
if (!nzchar(infile)) stop("Provide --file <..._overall_expression_limma.tsv>")
if (!file.exists(infile)) stop("File not found: ", infile)

p_thr <- suppressWarnings(as.numeric(get_arg("--p", "0.05")))
if (!is.finite(p_thr) || p_thr <= 0 || p_thr > 1) p_thr <- 0.05

fc_thr <- suppressWarnings(as.numeric(get_arg("--fc", "0.2")))
if (!is.finite(fc_thr) || fc_thr < 0) fc_thr <- 0.5

label_n <- suppressWarnings(as.integer(get_arg("--label-n", "20")))
if (!is.finite(label_n) || label_n < 0) label_n <- 20L

use_fdr <- tolower(get_arg("--use-fdr", "false")) %in% c("true", "t", "1", "yes", "y")

# ggrepel (auto-install if missing)
if (label_n > 0) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(ggrepel))
}

df <- fread(infile)

need <- c("geneid", "logFC", "P.Value")
if (!all(need %in% names(df))) stop("Input must contain columns: geneid, logFC, P.Value")

stat_col <- if (use_fdr) {
  if (!("adj.P.Val" %in% names(df))) stop("Requested --use-fdr but adj.P.Val column not found.")
  "adj.P.Val"
} else {
  "P.Value"
}

df[, stat := suppressWarnings(as.numeric(get(stat_col)))]
df <- df[is.finite(stat) & stat > 0 & stat <= 1]
df[, neglog := -log10(stat)]

# categorize points
df[, Direction := fifelse(logFC > 0, "Up", "Down")]
df[, Sig := stat < p_thr & abs(logFC) >= fc_thr]

df[, Class := fifelse(!Sig, "Not significant",
               fifelse(Direction == "Up", "Up (sig)", "Down (sig)"))]

# pick labels: most significant among significant genes
top_lab <- data.table()
if (label_n > 0) {
  top_lab <- df[Sig == TRUE][order(stat)][1:min(label_n, .N)]
}

# Extract PDE_GENE and TOP_SNP from filename if possible
base <- basename(infile)
base <- sub("_overall_expression_limma\\.tsv$", "", base)
# base should now look like PDE3A_rs6433690
title_main <- paste0(base, " overall expression volcano")

# group counts: try to read from columns if present, else NA
# Your limma results table had optional metadata columns earlier; handle both cases.
group_line <- "Group2 (ALT/ALT) vs Group1 (ALT/REF + REF/REF)"
n_line <- ""
if ("group_counts" %in% names(df)) {
  # group_counts is repeated per row; use first non-NA
  gc <- df[!is.na(group_counts)][1]$group_counts
  if (!is.na(gc) && nzchar(gc)) n_line <- paste0("n: ", gc)
}

subtitle_txt <- group_line
if (nzchar(n_line)) subtitle_txt <- paste0(group_line, " | ", n_line)

# Volcano plot
p <- ggplot(df, aes(x = logFC, y = neglog)) +
  geom_point(aes(color = Class), alpha = 0.75, size = 1.5) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(p_thr), linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = c(
    "Up (sig)" = "red",
    "Down (sig)" = "blue",
    "Not significant" = "grey70"
  )) +
  theme_classic(base_size = 14) +
  labs(
    title = title_main,
    subtitle = subtitle_txt,
    x = "logFC (Group2 − Group1)",
    y = paste0("-log10(", stat_col, ")"),
    color = ""
  )

# Add labels with leader lines
if (label_n > 0 && nrow(top_lab) > 0) {
  p <- p +
    ggrepel::geom_text_repel(
      data = top_lab,
      aes(label = geneid),
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf,
      segment.color = "black",
      min.segment.length = 0
    )
}

# Output filename
outfile <- sub(
  "_overall_expression_limma\\.tsv$",
  "_overall_expression_volcano.png",
  infile
)

ggsave(outfile, p, width = 7.4, height = 5.4, dpi = 300)

cat("Input:   ", infile, "\n", sep = "")
cat("Stat:    ", stat_col, " (thr=", p_thr, ")\n", sep = "")
cat("Effect:  |logFC| >= ", fc_thr, "\n", sep = "")
cat("Labels:  ", if (label_n > 0) label_n else 0, "\n", sep = "")
cat("Output:  ", outfile, "\n", sep = "")

# Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/04_plot_volcano_overall_expression_by_top_pde_snp.R \
#   --file output/eqtl/results/inspect/pde/overall_expression_by_top_pde_snp/PDE7A/PDE7A_rs1353749_overall_expression_limma.tsv
