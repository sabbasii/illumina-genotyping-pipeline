#!/usr/bin/env Rscript

# Plots: p-value histogram + volcano for rs966221 (ADD_G, DOM_G)
# Inputs (from env): EXPR_OUT_DIR/assoc_rs966221/limma_{ADD_G,DOM_G}_results.csv
# Outputs: (see "# Outputs" at end)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
stopifnot(nzchar(EXPR_OUT_DIR))
assoc_dir <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
plot_dir  <- file.path(assoc_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: safe read and column checks ----
read_results <- function(path){
  if (!file.exists(path)) stop("[ERR] Missing results file: ", path)
  df <- read_csv(path, show_col_types = FALSE)
  req <- c("P.Value","logFC")
  missing <- setdiff(req, names(df))
  if (length(missing)) stop("[ERR] Missing columns in ", basename(path), ": ", paste(missing, collapse=", "))
  # Choose label col
  label_col <- if ("GeneSymbol" %in% names(df)) "GeneSymbol" else if ("ProbeID" %in% names(df)) "ProbeID" else NULL
  list(df=df, label_col=label_col)
}

# ---- plotting helpers (base R only) ----
save_p_hist <- function(df, outfile, title){
  p <- df$P.Value
  png(outfile, width=1200, height=900)
  par(mar=c(5,5,4,2)+0.1)
  hist(p, breaks=50, main=paste(title, "- P-value histogram"), xlab="P-value", ylab="Count")
  mtext(sprintf("n=%d", sum(!is.na(p))), side=3, line=0.5, adj=1, cex=0.9)
  dev.off()
}

save_volcano <- function(
  df, outfile, title,
  label_col = NULL,
  p1e2 = 0.01, p1e3 = 0.001,
  label_angle = 0   # degrees: 0, 30, 45, etc.
){
  y <- -log10(df$P.Value)

  # color/size by nominal P
  col <- ifelse(df$P.Value < p1e3, "red",
                ifelse(df$P.Value < p1e2, "orange", "gray40"))
  cex <- ifelse(df$P.Value < p1e3, 1.2, ifelse(df$P.Value < p1e2, 0.9, 0.6))

  # add a bit of room so labels remain inside bounds
  x_rng <- range(df$logFC, na.rm=TRUE)
  y_rng <- range(y, na.rm=TRUE)
  x_pad <- 0.08 * diff(x_rng); if (!is.finite(x_pad)) x_pad <- 1
  y_pad <- 0.08 * diff(y_rng); if (!is.finite(y_pad)) y_pad <- 1
  xlim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  ylim <- c(y_rng[1], y_rng[2] + y_pad)

  png(outfile, width=1500, height=1050)  # a touch bigger for readability
  par(mar=c(5,6,4,2)+0.1)
  plot(df$logFC, y,
       pch=20, cex=cex, col=col,
       xlab="logFC per allele (effect size)",
       ylab="-log10(P)",
       main=paste(title, "- Volcano"),
       xlim=xlim, ylim=ylim)

  # Threshold lines (nominal P)
  abline(h=-log10(p1e2), lty=2)  # p=0.01
  abline(h=-log10(p1e3), lty=2)  # p=0.001

  # FDR=0.05 line (if adj.P.Val exists), label kept inside figure
  if ("adj.P.Val" %in% names(df)) {
    fdr_cut <- 0.05
    yh <- -log10(fdr_cut)
    abline(h = yh, lty=3)
    tx <- xlim[2] - 0.02*diff(xlim)
    ty <- yh + 0.02*diff(ylim)
    text(tx, ty, "FDR=0.05", cex=1.1, font=2, adj=c(1,0))
  }

  # Legend (color meaning)
  legend("topright",
         legend = c("P < 0.001", "0.001 ≤ P < 0.01", "P ≥ 0.01"),
         col    = c("red", "orange", "gray40"),
         pch    = 20,
         pt.cex = c(1.2, 0.9, 0.6),
         bty    = "n",
         cex    = 1.05)

  # Label ALL red points, labels NEXT TO dots (no leader lines)
  if (!is.null(label_col) && label_col %in% names(df)) {
    red_idx <- which(df$P.Value < p1e3)
    if (length(red_idx)) {
      # emphasize red points with thin outline
      points(df$logFC[red_idx], y[red_idx], pch=21, bg="red", col="black", cex=1.3, lwd=0.7)
      # put labels on the side: right for +logFC, left for -logFC; allow rotation
      pos_vec <- ifelse(df$logFC[red_idx] >= 0, 4, 2)  # 4=right, 2=left
      text(df$logFC[red_idx], y[red_idx],
           labels = df[[label_col]][red_idx],
           pos = pos_vec, offset = 0.55, cex = 1.2, srt = label_angle)
    }
  }

  # Summary
  msg <- sprintf("n=%d | p<0.01: %d | p<0.001: %d",
                 nrow(df),
                 sum(df$P.Value<0.01, na.rm=TRUE),
                 sum(df$P.Value<0.001, na.rm=TRUE))
  if ("adj.P.Val" %in% names(df)) {
    msg <- paste0(msg, sprintf(" | FDR<0.05: %d", sum(df$adj.P.Val < 0.05, na.rm=TRUE)))
  }
  mtext(msg, side=3, line=0.5, adj=1, cex=1.0)

  dev.off()
}

# ---- run for both models ----
files <- list(
  ADD  = file.path(assoc_dir, "limma_ADD_G_adjcov_results.csv"),  # limma_ADD_G_results.csv / limma_ADD_G_adjcov_results.csv
  DOM  = file.path(assoc_dir, "limma_DOM_G_adjcov_results.csv")   # limma_DOM_G_results.csv / limma_DOM_G_adjcov_results.csv
)

# configure label angle here (degrees)
LABEL_ANGLE <- 30

for (name in names(files)) {
  obj <- read_results(files[[name]])
  df  <- obj$df
  lab <- obj$label_col

  # Drop rows with missing P.Value/logFC
  df <- df[complete.cases(df$P.Value, df$logFC), , drop=FALSE]

  # === Save red + orange hits (P < 0.01) to CSV ===
  p1e2 <- 0.01
  p1e3 <- 0.001
  tier <- ifelse(df$P.Value < p1e3, "red",
                 ifelse(df$P.Value < p1e2, "orange", NA))
  keep_idx <- which(!is.na(tier))

  if (length(keep_idx)) {
    # Choose label column (fallbacks)
    label_vec <-
      if (!is.null(lab) && lab %in% names(df)) df[[lab]] else
      if ("ProbeID" %in% names(df)) df[["ProbeID"]] else
      as.character(seq_len(nrow(df)))

    out_df <- df[keep_idx, , drop=FALSE] |>
      mutate(
        Label = label_vec[keep_idx],
        tier  = tier[keep_idx],
        neglog10P = -log10(P.Value)
      ) |>
      # put useful columns up front; keep others if present
      select(Label, tier, P.Value,
             adj.P.Val = any_of("adj.P.Val"),
             logFC, neglog10P, everything()) |>
      arrange(P.Value, desc(abs(logFC)))

    out_csv <- file.path(assoc_dir, paste0(name, "_sig_hits_p001_p01.csv"))
    write_csv(out_df, out_csv)
  }

  # Volcano
  save_volcano(df,
               file.path(plot_dir, paste0(name, "_volcano_adj.png")),
               paste(name, "model"),
               label_col = lab,
               label_angle = LABEL_ANGLE)
}

cat("[OK] Plots written to: ", plot_dir, "\n")
cat("[OK] Sig-hit tables (if any) written next to: ", assoc_dir, "\n")

# ----------------------------
# Outputs
# ----------------------------
# Plots:
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/ADD_volcano_adj.png
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/DOM_volcano_adj.png
#
# Significant hits (red + orange; P < 0.01):
#   <EXPR_OUT_DIR>/assoc_rs966221/ADD_sig_hits_p001_p01.csv
#   <EXPR_OUT_DIR>/assoc_rs966221/DOM_sig_hits_p001_p01.csv
#
# (Optional, if you enable p-hist):
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/ADD_p_hist.png
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/DOM_p_hist.png
#
# To run:
#   Rscript scripts/rs966221_plot_volcano.R
