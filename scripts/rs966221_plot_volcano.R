#!/usr/bin/env Rscript

# Plots: p-value histogram + volcano for rs966221 (ADD_G, DOM_G)
# Inputs (from env): EXPR_OUT_DIR/assoc_rs966221/limma_{ADD_G,DOM_G}_results.csv
# Outputs: EXPR_OUT_DIR/assoc_rs966221/plots/*.png

suppressPackageStartupMessages({
  library(readr)
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

save_volcano <- function(df, outfile, title, label_col=NULL, p1e2=0.01, p1e3=0.001){
  # Color by nominal p thresholds
  y <- -log10(df$P.Value)
  col <- ifelse(df$P.Value < p1e3, "red",
                ifelse(df$P.Value < p1e2, "orange", "gray40"))
  cex <- ifelse(df$P.Value < p1e3, 1.2, ifelse(df$P.Value < p1e2, 0.9, 0.6))

  png(outfile, width=1200, height=900)
  par(mar=c(5,5,4,2)+0.1)
  plot(df$logFC, y,
       pch=20, cex=cex, col=col,
       xlab="logFC per allele (effect size)",
       ylab="-log10(P)",
       main=paste(title, "- Volcano"))
  # Threshold lines
  abline(h=-log10(p1e2), lty=2)  # p=0.01
  abline(h=-log10(p1e3), lty=2)  # p=0.001

  # Annotate top 10 by P (if we have labels)
  if (!is.null(label_col) && label_col %in% names(df)) {
    ord <- order(df$P.Value)
    top_idx <- ord[seq_len(min(10, length(ord)))]
    with(df[top_idx, ], {
      text(logFC, -log10(P.Value), labels=get(label_col), pos=4, cex=0.8, offset=0.5)
    })
  }
  mtext(sprintf("n=%d | p<0.01: %d | p<0.001: %d",
                nrow(df), sum(df$P.Value<0.01, na.rm=TRUE), sum(df$P.Value<0.001, na.rm=TRUE)),
        side=3, line=0.5, adj=1, cex=0.9)
  dev.off()
}

# ---- run for both models ----
files <- list(
  ADD  = file.path(assoc_dir, "limma_ADD_G_results.csv"),
  DOM  = file.path(assoc_dir, "limma_DOM_G_results.csv")
)

for (name in names(files)) {
  obj <- read_results(files[[name]])
  df  <- obj$df
  lab <- obj$label_col

  # Drop rows with missing P.Value/logFC
  df <- df[complete.cases(df$P.Value, df$logFC), , drop=FALSE]

  # P-hist
  save_p_hist(df, file.path(plot_dir, paste0(name, "_p_hist.png")), paste(name, "ADD_G"[name=="ADD"], "DOM_G"[name=="DOM"]))
  # Volcano
  save_volcano(df, file.path(plot_dir, paste0(name, "_volcano.png")), paste(name, "model"), label_col=lab)
}

cat("[OK] Plots written to: ", plot_dir, "\n")
