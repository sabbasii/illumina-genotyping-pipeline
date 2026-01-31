#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

# ---- paths / inputs ----
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
if (!nzchar(EXPR_OUT_DIR)) stop("EXPR_OUT_DIR is not set. Did you 'source scripts/00_config.sh'?")

OUTDIR   <- file.path(dirname(EXPR_OUT_DIR), "explore_batch")
expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("[I] EXPR_OUT_DIR:", EXPR_OUT_DIR, "\n")
cat("[I] OUTDIR      :", OUTDIR, "\n")
cat("[I] expr_csv    :", expr_csv, "\n")
cat("[I] meta_csv    :", meta_csv, "\n\n")

stopifnot(file.exists(expr_csv), file.exists(meta_csv))

# ---- load data ----
expr <- readr::read_csv(expr_csv, guess_max = 200000, show_col_types = FALSE)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

# Find first UASG sample column by name (robust to annotation width)
sample_idx <- which(startsWith(colnames(expr), "UASG-"))
if (length(sample_idx) == 0) stop("No columns starting with 'UASG-' found in expression file.")

# Expression matrix: features x samples
X <- expr[, sample_idx, drop = FALSE] %>% as.data.frame()
X[] <- lapply(X, function(col) suppressWarnings(as.numeric(col)))

# ---- fix duplicate header names in metadata ----
# pick the FIRST column whose name starts with "UASG" and call it "UASG"
u_idx <- grep("^UASG(\\b|\\.{3}\\d+)", names(meta))
if (length(u_idx) < 1) stop("Could not find any 'UASG' column in metadata (even with duplicates).")
if (length(u_idx) > 1) message("[I] Multiple UASG-like columns detected: using the FIRST (", names(meta)[u_idx[1]], ").")
meta$UASG <- meta[[u_idx[1]]]

# similarly, standardize Scan Date → ScanDate (use first if duplicated)
sd_idx <- which(names(meta) == "Scan Date" | grepl("^Scan Date(\\.|\\.{3}\\d+)", names(meta)))
if (length(sd_idx) < 1) stop("Metadata is missing a 'Scan Date' column (or recognizable variant).")
if (length(sd_idx) > 1) message("[I] Multiple 'Scan Date' columns detected: using the FIRST (", names(meta)[sd_idx[1]], ").")
meta$ScanDate <- meta[[sd_idx[1]]]

# ---- align samples with metadata using the standardized UASG ----
common_ids <- intersect(colnames(X), meta$UASG)
if (length(common_ids) < 2) stop("Fewer than 2 overlapping UASG IDs between expression and metadata.")

# Subset and reorder to metadata order
meta_sub <- meta %>% dplyr::filter(UASG %in% common_ids)
X <- X[, meta_sub$UASG, drop = FALSE]  # columns ordered like meta_sub$UASG

# ---- quick sanity: log scale and variance filtering ----
x_max <- suppressWarnings(max(X, na.rm = TRUE))
did_log <- FALSE
if (is.finite(x_max) && x_max > 100) {
  X <- log2(X + 1)
  did_log <- TRUE
  cat("[I] Applied log2(x+1) transform (max value was >", x_max, ").\n")
}

# Remove rows (features) with zero or NA variance
var_ok <- apply(X, 1, function(r) {
  v <- stats::var(r, na.rm = TRUE)
  is.finite(v) && v > 0
})
Xf <- X[var_ok, , drop = FALSE]
cat("[I] Filtered features: kept", sum(var_ok), "of", length(var_ok), "with non-zero variance.\n")

# Replace remaining NAs with column medians (rare)
fill_na_with_col_med <- function(df) {
  for (j in seq_along(df)) {
    col <- df[[j]]
    if (anyNA(col)) {
      med <- stats::median(col, na.rm = TRUE)
      col[is.na(col)] <- med
      df[[j]] <- col
    }
  }
  df
}
Xf <- fill_na_with_col_med(Xf)

# ---- PCA (samples as rows) ----
Xt <- t(as.matrix(Xf))
pca <- prcomp(Xt, center = TRUE, scale. = TRUE)

# Variance explained
var_exp <- pca$sdev^2
var_exp <- var_exp / sum(var_exp)
var_df <- tibble::tibble(PC = paste0("PC", seq_along(var_exp)),
                         VarianceExplained = var_exp,
                         Percent = round(100 * var_exp, 2))
readr::write_csv(var_df, file.path(OUTDIR, "pca_variance.csv"))

# Scores (samples)
scores <- as.data.frame(pca$x)
scores$UASG <- rownames(scores)

# Merge Scan Date (batch) for coloring
if (!"Scan Date" %in% colnames(meta_sub)) {
  stop("Metadata is missing 'Scan Date' column; cannot assess batch by Scan Date.")
}
plot_df <- scores %>%
  left_join(meta_sub %>% select(UASG, `Scan Date`), by = "UASG") %>%
  rename(ScanDate = `Scan Date`)

# Save scores
readr::write_csv(plot_df, file.path(OUTDIR, "pca_scores_with_ScanDate.csv"))

# ---- Plots ----
# PC1 vs PC2 colored by ScanDate
p12 <- ggplot(plot_df, aes(PC1, PC2, color = ScanDate)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(title = "PCA: PC1 vs PC2 by Scan Date",
       subtitle = if (did_log) "log2(x+1), centered & scaled" else "centered & scaled",
       x = paste0("PC1 (", round(100*var_df$VarianceExplained[1], 2), "%)"),
       y = paste0("PC2 (", round(100*var_df$VarianceExplained[2], 2), "%)")) +
  theme_minimal(base_size = 12)
ggsave(file.path(OUTDIR, "PCA_PC1_PC2_by_ScanDate.png"), p12, width = 7, height = 5, dpi = 150)

# Scree plot
scree <- ggplot(var_df[1:20, ], aes(x = PC, y = Percent, group = 1)) +
  geom_col() + geom_point() + geom_line() +
  labs(title = "Scree: percent variance explained", y = "Percent", x = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTDIR, "PCA_scree.png"), scree, width = 8, height = 4.5, dpi = 150)

# ---- PERMANOVA (irrespective of other factors) ----
perm_out_txt <- file.path(OUTDIR, "permanova_scan_date.txt")
perm_done <- FALSE
suppressWarnings({
  if (requireNamespace("vegan", quietly = TRUE)) {
    library(vegan)
    # Euclidean distance on the PCA input (Xt is already centered/scaled by prcomp, but we used prcomp on Xt)
    # Use scaled Xt again to avoid leakage; scale here explicitly:
    Xt_scaled <- scale(Xt, center = TRUE, scale = TRUE)
    d <- dist(Xt_scaled, method = "euclidean")
    md <- plot_df %>% select(UASG, ScanDate)
    perm <- adonis2(d ~ ScanDate, data = md, permutations = 999)
    capture.output(perm, file = perm_out_txt)
    perm_done <- TRUE
    cat("[I] PERMANOVA completed. See:", perm_out_txt, "\n")
  } else {
    cat("[W] Package 'vegan' not installed; skipping PERMANOVA. Install with: install.packages('vegan')\n")
  }
})

# ---- PC-wise tests (ANOVA + Kruskal) on first 10 PCs ----
pc_n <- min(10, ncol(scores) - 1L)  # -1 for UASG column
anova_rows   <- list()
kruskal_rows <- list()

for (k in seq_len(pc_n)) {
  pc_name <- paste0("PC", k)
  dfk <- plot_df %>% select(ScanDate, !!pc_name)
  dfk$ScanDate <- as.factor(dfk$ScanDate)

  # ANOVA
  aov_res <- tryCatch({
    summary(aov(reformulate("ScanDate", response = pc_name), data = dfk))[[1]]
  }, error = function(e) NA)

  if (is.data.frame(aov_res) && "Pr(>F)" %in% colnames(aov_res)) {
    pval <- aov_res[1, "Pr(>F)"]
  } else {
    pval <- NA_real_
  }
  anova_rows[[k]] <- data.frame(PC = pc_name, Test = "ANOVA", p_value = pval, stringsAsFactors = FALSE)

  # Kruskal-Wallis
  kw <- tryCatch({
    kruskal.test(reformulate("ScanDate", response = pc_name), data = dfk)
  }, error = function(e) NA)
  kw_p <- if (inherits(kw, "htest")) kw$p.value else NA_real_
  kruskal_rows[[k]] <- data.frame(PC = pc_name, Test = "Kruskal", p_value = kw_p, stringsAsFactors = FALSE)
}

pc_tests <- bind_rows(do.call(rbind, anova_rows), do.call(rbind, kruskal_rows)) %>%
  arrange(PC, Test)
readr::write_csv(pc_tests, file.path(OUTDIR, "pc_batch_tests_scan_date.csv"))

# ---- final console notes ----
cat("\n=== Batch check summary ===\n")
cat("* Samples (n):", nrow(scores), "\n")
cat("* Features used (non-zero var):", nrow(Xf), "\n")
cat("* Log transform applied:", did_log, "\n")
if (perm_done) cat("* PERMANOVA: see permanova_scan_date.txt\n")
cat("* PC tests:  see pc_batch_tests_scan_date.csv\n")
cat("* Plots:     PCA_PC1_PC2_by_ScanDate.png, PCA_scree.png\n")

# --------- MAKE IT EXECUTABLE ---------
# chmod +x scripts/expr/check_batch_effect_scan_date.R
# --------- RUN IT ---------
# Rscript scripts/expr/check_batch_effect_scan_date.R
