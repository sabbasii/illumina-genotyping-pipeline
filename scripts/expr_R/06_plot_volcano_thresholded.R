#!/usr/bin/env Rscript
# Thresholded + Annotated Volcano Plot
# Reads: $OUT_DIR/expr/dge_unadjusted/dge_unadjusted_all.csv
# Writes: volcano_thresholded_[raw|fdr]_fc{X}_p{Y}.png/.pdf

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(ggrepel)
  library(stringr); library(tibble)
})

# ---------------- paths ----------------
OUT_DIR <- Sys.getenv("OUT_DIR")
stopifnot(nzchar(OUT_DIR))
INDIR   <- file.path(OUT_DIR, "expr", "dge_unadjusted")
INCSV   <- file.path(INDIR, "dge_unadjusted_all.csv")

OUTDIR  <- INDIR  # write next to results
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("[I] Input  :", INCSV, "\n")
cat("[I] Out dir:", OUTDIR, "\n")
stopifnot(file.exists(INCSV))

# ---------------- args (optional) ----------------
# Args: [1]=use_fdr ("fdr" or "raw"), [2]=fc_abs_threshold (default 1),
#       [3]=p_threshold (default 0.05), [4]=label_top_n (default 20)
args <- commandArgs(trailingOnly = TRUE)
use_fdr        <- if (length(args) >= 1) tolower(args[1]) else "raw"  # "raw" or "fdr"
fc_thresh      <- if (length(args) >= 2) as.numeric(args[2]) else 1.0
p_thresh_input <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
label_top_n    <- if (length(args) >= 4) as.integer(args[4]) else 20L

if (!use_fdr %in% c("raw","fdr")) use_fdr <- "raw"
if (is.na(fc_thresh) || fc_thresh <= 0) fc_thresh <- 1.0
if (is.na(label_top_n) || label_top_n < 0) label_top_n <- 20L

# ---------------- load data ----------------
df <- readr::read_csv(INCSV, show_col_types = FALSE)

# Requi
 columns sanity
req_cols <- c("FeatureID","GeneSymbol","logFC","P.Value","adj.P.Val")
stopifnot(all(req_cols %in% names(df)))

# Label column fallback
label_col <- if ("GeneSymbol" %in% names(df)) "GeneSymbol" else "FeatureID"

# Pick p metric
if (use_fdr == "fdr") {
  p_col <- "adj.P.Val"
  p_thresh <- p_thresh_input
  ylab_txt <- "-log10(FDR)"
  df <- df %>% mutate(neglog10p = -log10(pmax(!!as.name(p_col), .Machine$double.xmin)))
} else {
  p_col <- "P.Value"
  p_thresh <- p_thresh_input
  ylab_txt <- "-log10(p-value)"
  df <- df %>% mutate(neglog10p = -log10(pmax(!!as.name(p_col), .Machine$double.xmin)))
}

# Significance calls
df <- df %>%
  mutate(
    sig_fc   = abs(logFC) >= fc_thresh,
    sig_p    = !!as.name(p_col) < p_thresh,
    status   = case_when(
      sig_fc & sig_p & logFC > 0  ~ "Up",
      sig_fc & sig_p & logFC < 0  ~ "Down",
      TRUE                        ~ "NS"
    )
  )

# Choose top labels: prioritize significant; else top by -log10p
df <- df %>%
  mutate(label_use = dplyr::case_when(
    status != "NS" ~ !!as.name(label_col),
    TRUE           ~ NA_character_
  ))

# If fewer than label_top_n significant, fill with best by p-value
sig_n <- sum(df$status != "NS", na.rm = TRUE)
if (label_top_n > 0) {
  if (sig_n >= label_top_n) {
    lab_idx <- which(df$status != "NS")
    # keep the top label_top_n by neglog10p among significant
    keep <- order(df$neglog10p[lab_idx], decreasing = TRUE)[seq_len(min(label_top_n, length(lab_idx)))]
    keep_idx <- lab_idx[keep]
  } else {
    # take all significant + (label_top_n - sig_n) best remaining by p
    non_idx <- which(df$status == "NS")
    fill_n  <- max(0, label_top_n - sig_n)
    fill_idx <- if (fill_n > 0) non_idx[order(df$neglog10p[non_idx], decreasing = TRUE)][seq_len(min(fill_n, length(non_idx)))] else integer(0)
    keep_idx <- c(which(df$status != "NS"), fill_idx)
  }
  lab_vec <- rep(NA_character_, nrow(df))
  lab_vec[keep_idx] <- df[[label_col]][keep_idx]
  df$label_use <- lab_vec
}

# Axis and thresholds
xlab_txt <- "log2 Fold Change"
vline_pos <- c(-fc_thresh, fc_thresh)
hline_pos <- -log10(p_thresh)

# Colors
cols <- c("Down" = "#377eb8", "Up" = "#e41a1c", "NS" = "grey70")

# portable, no custom font
base_theme <- theme_minimal(base_size = 12)

title_txt <- sprintf("Volcano (unadjusted) — %s threshold",
                     ifelse(use_fdr=="fdr","FDR","p-value"))
subtitle_txt <- sprintf("|log2FC|≥%.2f, %s<%.3g",
                        fc_thresh, ifelse(use_fdr=="fdr","FDR","p"), p_thresh)

p <- ggplot(df, aes(x = logFC, y = neglog10p, color = status)) +
  geom_point(alpha = 0.8, size = 1.6) +
  scale_color_manual(values = cols, breaks = c("Up","Down","NS")) +
  geom_vline(xintercept = vline_pos, linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = hline_pos, linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = subset(df, !is.na(label_use)),
    aes(label = label_use),
    size = 3, max.overlaps = 100, min.segment.length = 0,
    box.padding = 0.25, point.padding = 0.15, seed = 42
  ) +
  labs(
    title = title_txt,
    subtitle = subtitle_txt,
    x = xlab_txt, y = ylab_txt, color = "Status"
  ) +
  base_theme

# filenames
tag <- sprintf("%s_fc%.2f_p%.3g", ifelse(use_fdr=="fdr","fdr","raw"), fc_thresh, p_thresh)
pngf <- file.path(OUTDIR, paste0("volcano_thresholded_", tag, ".png"))
# pdff <- file.path(OUTDIR, paste0("volcano_thresholded_", tag, ".pdf"))

ggsave(pngf, p, width = 7.5, height = 5.5, dpi = 200)
ggsave(pdff, p, width = 7.5, height = 5.5)

cat("[I] Wrote:\n - ", pngf, "\n - ", pdff, "\n", sep = "")

# ------ Installation ------
## mamba install -c conda-forge r-ggrepel

# ------ Run it ------
## Rscript scripts/expr/06_plot_volcano_thresholded.R


# ------------------------------------------------------------
# NOTES (how to use / where results go)
#
# • Coloring logic:
#     - "Up"  = (abs(logFC) >= fc_thresh) AND (p metric < p_thresh) AND logFC > 0
#     - "Down"= (abs(logFC) >= fc_thresh) AND (p metric < p_thresh) AND logFC < 0
#     - "NS"  = otherwise. If you see only grey dots, no gene passed BOTH thresholds.
#
# • Thresholds (CLI args):
#     arg1 = "raw" (use P.Value) or "fdr" (use adj.P.Val)
#     arg2 = fc_thresh (abs log2FC cutoff), default 1.0  # ≈2x
#     arg3 = p_thresh  (p or FDR cutoff),  default 0.05
#     arg4 = label_top_n (how many points to label), default 20
#
# • Examples:
#     # raw p < 0.05 and |log2FC| ≥ 1.0
#     # Rscript scripts/expr/plot_volcano_thresholded.R
#     # raw p < 0.05 and |log2FC| ≥ 0.58 (≈1.5x)
#     # Rscript scripts/expr/plot_volcano_thresholded.R raw 0.58 0.05 25
#     # FDR < 0.10 and |log2FC| ≥ 0.58
#     # Rscript scripts/expr/plot_volcano_thresholded.R fdr 0.58 0.10 25
#
# • Output files (next to DGE results):
#     $OUT_DIR/expr/dge_unadjusted/volcano_thresholded_<raw|fdr>_fc<FC>_p<P>.png
#     $OUT_DIR/expr/dge_unadjusted/volcano_thresholded_<raw|fdr>_fc<FC>_p<P>.pdf
#
# • Troubleshooting:
#     - Only grey points? Relax fc_thresh (e.g., 0.58) or p/FDR threshold.
#     - Font error ("invalid font type")? Remove any custom font theme lines.
# ------------------------------------------------------------
