#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# 01_qq_eqtl.R
#
# QQ plots for MatrixEQTL output tables.
# Produces one QQ plot per results file (all / cis / trans) using the
# p-value column, with invalid values filtered out.
#
# Inputs (under REPO_ROOT/output/eqtl/results/):
#   - eqtl_all.tsv
#   - eqtl_cis.tsv
#   - eqtl_trans.tsv
#
# Outputs (under REPO_ROOT/output/eqtl/results/inspect/qc/):
#   - qq__all.png
#   - qq__cis.png
#   - qq__trans.png
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

source(file.path(REPO_ROOT, "scripts/eqtl/utils/inspect_helpers.R"))
paths <- get_eqtl_paths()

read_pvals <- function(path) {
  assert_file(path)
  x <- fread(path)

  p_col <- pick_p_col(x)
  p <- suppressWarnings(as.numeric(x[[p_col]]))
  p <- p[is.finite(p) & p > 0 & p <= 1]

  p
}

qq_points <- function(p) {
  p <- sort(p)
  n <- length(p)
  exp_p <- -log10(ppoints(n))
  obs_p <- -log10(p)
  list(exp = exp_p, obs = obs_p, n = n)
}

plot_one <- function(label, p, out_png) {
  q <- qq_points(p)

  png(out_png, width = 900, height = 900, res = 150)
  plot(
    q$exp, q$obs,
    xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)",
    main = paste0("QQ plot: ", label, " (n=", q$n, ")"),
    pch = 16, cex = 0.5
  )
  abline(0, 1)
  dev.off()
}

files <- list(
  all   = paths$eqtl_all,
  cis   = paths$eqtl_cis,
  trans = paths$eqtl_trans
)

for (nm in names(files)) {
  if (!file.exists(files[[nm]])) {
    message("Skipping ", nm, " (missing: ", files[[nm]], ")")
    next
  }

  p <- read_pvals(files[[nm]])
  if (length(p) == 0) {
    message("Skipping ", nm, " (no valid p-values).")
    next
  }

  out_png <- file.path(paths$qc_dir, paste0("qq__", nm, ".png"))
  plot_one(nm, p, out_png)
  message("Wrote: ", out_png)
}

message("Done.")

# ---- Run ----
# source scripts/00_config.sh
# Rscript scripts/eqtl/05_inspect_results/01_qq_eqtl.R