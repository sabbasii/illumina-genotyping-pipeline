#!/usr/bin/env Rscript

# rs966221 volcano + p<0.05 hit tables (ADD and DOM models)
# ---------------------------------------------------------
# - Reads limma results for ADD_G and DOM_G models
# - Makes volcano plots (all points gray; darker if P < 0.05)
# - Saves one CSV per model with all genes where P < 0.05
#
# Inputs (from env):
#   $EXPR_OUT_DIR/assoc_rs966221/limma_ADD_G_adjcov_results.csv
#   $EXPR_OUT_DIR/assoc_rs966221/limma_DOM_G_adjcov_results.csv
#
# Outputs:
#   Plots:
#     $EXPR_OUT_DIR/assoc_rs966221/plots/ADD_volcano_adj.png
#     $EXPR_OUT_DIR/assoc_rs966221/plots/DOM_volcano_adj.png
#   Hit tables (P < 0.05):
#     $EXPR_OUT_DIR/assoc_rs966221/ADD_p_lt_0_05.csv
#     $EXPR_OUT_DIR/assoc_rs966221/DOM_p_lt_0_05.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- paths and directories -----------------------------------------------

EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
if (!nzchar(EXPR_OUT_DIR)) {
  stop("EXPR_OUT_DIR is not set. Did you source scripts/00_config.sh?")
}

assoc_dir <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
plot_dir  <- file.path(assoc_dir, "plots")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: read limma results safely -----------------------------------
# - Checks file exists
# - Ensures required columns are present
# - Chooses a label column for plotting (GeneSymbol or ProbeID)

read_results <- function(path) {
  if (!file.exists(path)) {
    stop("[ERR] Missing results file: ", path)
  }

  df <- read_csv(path, show_col_types = FALSE)

  req <- c("P.Value", "logFC")
  missing <- setdiff(req, names(df))
  if (length(missing)) {
    stop(
      "[ERR] Missing columns in ", basename(path), ": ",
      paste(missing, collapse = ", ")
    )
  }

  label_col <- if ("GeneSymbol" %in% names(df)) {
    "GeneSymbol"
  } else if ("ProbeID" %in% names(df)) {
    "ProbeID"
  } else {
    NULL
  }

  list(df = df, label_col = label_col)
}

# ---- helper: P-value histogram (optional QC) -----------------------------
# Uses base R hist() to inspect the P-value distribution.

save_p_hist <- function(df, outfile, title) {
  p <- df$P.Value

  png(outfile, width = 1200, height = 900)
  par(mar = c(5, 5, 4, 2) + 0.1)

  hist(
    p,
    breaks = 50,
    main   = paste(title, "- P-value histogram"),
    xlab   = "P-value",
    ylab   = "Count"
  )

  mtext(
    sprintf("n=%d", sum(!is.na(p))),
    side = 3, line = 0.5, adj = 1, cex = 0.9
  )

  dev.off()
}

# ---- helper: volcano plot (base R only) ----------------------------------
# - x-axis: logFC
# - y-axis: -log10(P)
# - Colors by P:
#     P < 0.001           -> red
#     0.001 <= P < 0.01   -> orange
#     0.01  <= P < 0.05   -> green
#     P >= 0.05           -> gray
# - Optional labels for very small P-values (P < 0.001)

save_volcano <- function(
  df,
  outfile,
  title,
  label_col   = NULL,
  label_angle = 0   # text rotation in degrees
) {
  # y-axis: -log10(P)
  y <- -log10(df$P.Value)

  # Color by P-value tier
  col <- ifelse(
    df$P.Value < 0.001, "red",
    ifelse(
      df$P.Value < 0.01, "orange",
      ifelse(
        df$P.Value < 0.05, "green",
        "gray75"
      )
    )
  )

  # Point size
  cex <- 0.7

  # Padding so points and labels are not clipped
  x_rng <- range(df$logFC, na.rm = TRUE)
  y_rng <- range(y, na.rm = TRUE)
  x_pad <- 0.08 * diff(x_rng)
  if (!is.finite(x_pad)) x_pad <- 1
  y_pad <- 0.08 * diff(y_rng)
  if (!is.finite(y_pad)) y_pad <- 1

  xlim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  ylim <- c(y_rng[1], y_rng[2] + y_pad)

  png(outfile, width = 1500, height = 1050)
  par(mar = c(5, 6, 4, 2) + 0.1)

  plot(
    df$logFC, y,
    pch = 20, cex = cex, col = col,
    xlab = "logFC per allele (effect size)",
    ylab = "-log10(P)",
    main = paste(title, "- Volcano"),
    xlim = xlim, ylim = ylim
  )

  # P-value threshold lines and labels
  p_cuts <- c(0.05, 0.01, 0.001)
  cut_labs <- c("p = 0.05", "p = 0.01", "p = 0.001")

  for (i in seq_along(p_cuts)) {
    yc <- -log10(p_cuts[i])
    abline(h = yc, lty = 2)
    tx <- xlim[2] - 0.02 * diff(xlim)
    ty <- yc + 0.02 * diff(ylim)
    text(
      tx, ty,
      cut_labs[i],
      cex = 0.95,
      adj = c(1, 0)
    )
  }

  # Legend for color tiers
  legend(
    "topright",
    legend = c(
      "P < 0.001",
      "0.001 <= P < 0.01",
      "0.01 <= P < 0.05",
      "P >= 0.05"
    ),
    col    = c("red", "orange", "green", "gray75"),
    pch    = 20,
    pt.cex = rep(0.7, 4),
    bty    = "n",
    cex    = 1.05
  )

  # Label points with very small P-values (P < 0.001), if we have labels
  if (!is.null(label_col) && label_col %in% names(df)) {
    lab_idx <- which(df$P.Value < 0.001)
    if (length(lab_idx)) {
      # Emphasize these points with a slightly larger outlined marker
      points(
        df$logFC[lab_idx], y[lab_idx],
        pch = 21,
        bg  = "red",
        col = "black",
        cex = 1.1,
        lwd = 0.7
      )

      # Right for positive logFC, left for negative logFC
      pos_vec <- ifelse(df$logFC[lab_idx] >= 0, 4, 2)

      text(
        df$logFC[lab_idx], y[lab_idx],
        labels = df[[label_col]][lab_idx],
        pos    = pos_vec,
        offset = 0.55,
        cex    = 1.2,
        srt    = label_angle
      )
    }
  }

  # Summary line at the top
  msg <- sprintf(
    "n=%d | p<0.05: %d | p<0.01: %d | p<0.001: %d",
    nrow(df),
    sum(df$P.Value < 0.05,  na.rm = TRUE),
    sum(df$P.Value < 0.01,  na.rm = TRUE),
    sum(df$P.Value < 0.001, na.rm = TRUE)
  )

  if ("adj.P.Val" %in% names(df)) {
    msg <- paste0(
      msg,
      sprintf(
        " | FDR<0.05: %d",
        sum(df$adj.P.Val < 0.05, na.rm = TRUE)
      )
    )
  }

  mtext(msg, side = 3, line = 0.5, adj = 1, cex = 1.0)

  dev.off()
}

# ---- input file list for both models -------------------------------------
# Use simple names "ADD" and "DOM" for the two genotype models.

files <- list(
  ADD = file.path(
    assoc_dir,
    "limma_ADD_G_adjcov_results.csv"
  ),
  DOM = file.path(
    assoc_dir,
    "limma_DOM_G_adjcov_results.csv"
  )
)

# Label rotation angle for volcano labels (degrees).
LABEL_ANGLE <- 10

# ---- main loop: process ADD and DOM results ------------------------------
# For each model:
#  - read limma results
#  - drop rows with missing P.Value/logFC
#  - save a CSV of all genes with P < 0.05
#  - make a volcano plot

for (name in names(files)) {
  obj <- read_results(files[[name]])
  df  <- obj$df
  lab <- obj$label_col

  # Keep only rows with non-missing P.Value and logFC.
  df <- df[
    complete.cases(df$P.Value, df$logFC),
    ,
    drop = FALSE
  ]

  # ---- save all genes with P < 0.05 to CSV -------------------------------

  p_cut <- 0.05
  hit_idx <- which(df$P.Value < p_cut & !is.na(df$P.Value))

  if (length(hit_idx)) {
    hits_df <- df[hit_idx, , drop = FALSE] |>
      mutate(
        neglog10P = -log10(P.Value)
      ) |>
      select(
        any_of(c("GeneSymbol", "ProbeID")),
        P.Value,
        adj.P.Val = any_of("adj.P.Val"),
        logFC,
        neglog10P,
        everything()
      ) |>
      arrange(P.Value, desc(abs(logFC)))

    out_csv <- file.path(
      assoc_dir,
      paste0(name, "_p_lt_0_05.csv")
    )

    write_csv(hits_df, out_csv)
  }

  # ---- volcano plot ------------------------------------------------------

  out_png <- file.path(
    plot_dir,
    paste0(name, "_volcano_adj.png")
  )

  save_volcano(
    df          = df,
    outfile     = out_png,
    title       = paste(name, "model"),
    label_col   = lab,
    label_angle = LABEL_ANGLE
  )
}

cat("[OK] Plots written to: ", plot_dir, "\n")
cat("[OK] P<0.05 hit tables written in: ", assoc_dir, "\n")

# To run from the project root:
#   Rscript scripts/rs966221_plot_volcano.R

##########################################################
##########################################################
##########################################################
##########################################################
##########################################################

# ---- extra: volcano with FC (fold change) on x-axis ----------------------
# Uses FC = 2^logFC, keeps the same P-value color tiers as before.
# Output files:
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/ADD_volcano_fc_adj.png
#   <EXPR_OUT_DIR>/assoc_rs966221/plots/DOM_volcano_fc_adj.png

for (name in names(files)) {
  obj <- read_results(files[[name]])
  df  <- obj$df

  # Keep rows with valid P.Value and logFC
  df <- df[
    complete.cases(df$P.Value, df$logFC),
    ,
    drop = FALSE
  ]

  # Fold change: FC = 2^logFC
  df$FC <- 2 ^ df$logFC

  y <- -log10(df$P.Value)

  # Color by P tiers (same as logFC volcano)
  col <- ifelse(
    df$P.Value < 0.001, "red",
    ifelse(
      df$P.Value < 0.01, "orange",
      ifelse(
        df$P.Value < 0.05, "green",
        "gray75"
      )
    )
  )

  cex <- 0.7

  # Use a trimmed range for FC to avoid extreme outliers dominating the axis
  fc_q <- quantile(df$FC, probs = c(0.01, 0.99), na.rm = TRUE)
  x_rng <- as.numeric(fc_q)
  y_rng <- range(y, na.rm = TRUE)

  x_pad <- 0.08 * diff(x_rng)
  if (!is.finite(x_pad)) x_pad <- 0.1
  y_pad <- 0.08 * diff(y_rng)
  if (!is.finite(y_pad)) y_pad <- 1

  xlim <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  ylim <- c(y_rng[1], y_rng[2] + y_pad)

  out_png_fc <- file.path(
    plot_dir,
    paste0(name, "_volcano_fc_adj.png")
  )

  png(out_png_fc, width = 1500, height = 1050)
  par(mar = c(5, 6, 4, 2) + 0.1)

  plot(
    df$FC, y,
    pch = 20, cex = cex, col = col,
    xlab = "Fold change (2^logFC)",
    ylab = "-log10(P)",
    main = paste(name, "model - Volcano (FC axis)"),
    xlim = xlim, ylim = ylim
  )

  # P-value threshold lines and labels
  p_fc_cuts  <- c(0.05, 0.01, 0.001)
  p_fc_labs  <- c("p = 0.05", "p = 0.01", "p = 0.001")

  for (i in seq_along(p_fc_cuts)) {
    yc <- -log10(p_fc_cuts[i])
    abline(h = yc, lty = 2)
    tx <- xlim[2] - 0.02 * diff(xlim)
    ty <- yc + 0.02 * diff(ylim)
    text(
      tx, ty,
      p_fc_labs[i],
      cex = 0.95,
      adj = c(1, 0)
    )
  }

  # Legend
  legend(
    "topright",
    legend = c(
      "P < 0.001",
      "0.001 <= P < 0.01",
      "0.01 <= P < 0.05",
      "P >= 0.05"
    ),
    col    = c("red", "orange", "green", "gray75"),
    pch    = 20,
    pt.cex = rep(0.7, 4),
    bty    = "n",
    cex    = 1.05
  )

  # Summary line
  msg_fc <- sprintf(
    "n=%d | p<0.05: %d | p<0.01: %d | p<0.001: %d",
    nrow(df),
    sum(df$P.Value < 0.05,  na.rm = TRUE),
    sum(df$P.Value < 0.01,  na.rm = TRUE),
    sum(df$P.Value < 0.001, na.rm = TRUE)
  )

  if ("adj.P.Val" %in% names(df)) {
    msg_fc <- paste0(
      msg_fc,
      sprintf(
        " | FDR<0.05: %d",
        sum(df$adj.P.Val < 0.05, na.rm = TRUE)
      )
    )
  }

  mtext(msg_fc, side = 3, line = 0.5, adj = 1, cex = 1.0)

  dev.off()
}
