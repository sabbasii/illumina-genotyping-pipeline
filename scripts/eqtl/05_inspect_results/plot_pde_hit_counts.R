#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 06_plot_pde_hit_counts.R
#
# Bar plot: number of eQTL hits per PDE gene
#
# Input:
#   output/eqtl/results/inspect/pde/eqtl_all_pde.tsv
#
# Output:
#   output/eqtl/results/inspect/pde/pde_eqtl_hit_counts.png
#
# Run:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/05_inspect_results/plot_pde_hit_counts.R
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT not set")

in_file <- file.path(REPO_ROOT,
                     "output/eqtl/results/inspect/pde/eqtl_all_pde.tsv")

out_png <- file.path(REPO_ROOT,
                     "output/eqtl/results/inspect/pde/pde_eqtl_hit_counts.png")

x <- fread(in_file)

if (!("gene" %in% names(x))) stop("Missing column: gene")
if (nrow(x) == 0) stop("Input file is empty")

gene_counts <- x[, .(n_hits = .N), by = gene][order(-n_hits)]

p <- ggplot(gene_counts,
            aes(x = reorder(gene, n_hits), y = n_hits)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(
    title = "PDE genes ranked by number of eQTL associations",
    x = "",
    y = "Number of SNP–gene associations"
  )

ggsave(out_png, p, width = 7, height = 5, dpi = 300)

cat("Wrote:", out_png, "\n")