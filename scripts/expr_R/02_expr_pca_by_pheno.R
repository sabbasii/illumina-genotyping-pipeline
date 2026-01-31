#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# ---------------- paths ----------------
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
OUT_DIR      <- Sys.getenv("OUT_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(OUT_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")
OUTDIR   <- file.path(OUT_DIR, "expr", "pca_pheno")
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("[I] expr:", expr_csv, "\n")
cat("[I] meta:", meta_csv, "\n")
cat("[I] out :", OUTDIR, "\n\n")
stopifnot(file.exists(expr_csv), file.exists(meta_csv))

# ---------------- load data ----------------
expr <- readr::read_csv(expr_csv, guess_max = 200000, show_col_types = FALSE)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

# ---- canonicalize key metadata columns (handle duplicates) ----
pick_first <- function(df, pattern_exact, pattern_regex=NULL, fallback=NA_character_) {
  idx <- integer(0)
  if (!is.null(pattern_exact)) idx <- c(idx, which(names(df) == pattern_exact))
  if (!is.null(pattern_regex)) idx <- c(idx, grep(pattern_regex, names(df)))
  idx <- unique(idx)
  if (length(idx) >= 1) df[[idx[1]]] else rep(fallback, nrow(df))
}

# UASG
meta$UASG <- pick_first(meta, "UASG", "^UASG(\\.|\\.{3}\\d+)")
stopifnot(any(!is.na(meta$UASG)))

# StrokeStatus (derive from Final Diagnosis if absent)
if (!"StrokeStatus" %in% names(meta)) {
  FD <- pick_first(meta, "Final Diagnosis", "^Final Diagnosis(\\.|\\.{3}\\d+)")
  if (all(is.na(FD))) {
    meta$StrokeStatus <- NA_character_
    message("[W] No StrokeStatus/Final Diagnosis; set to NA.")
  } else {
    meta$StrokeStatus <- ifelse(!is.na(FD) & tolower(trimws(FD)) == "control", "Control", "Non-control")
    message("[I] Derived StrokeStatus from Final Diagnosis.")
  }
}

# Sex
meta$Sex <- pick_first(meta, "Sex", "^Sex(\\.|\\.{3}\\d+)")

# Age (numeric, prefer 'Age' then 'Age At Onset')
Age1 <- suppressWarnings(as.numeric(pick_first(meta, "Age", "^Age(\\.|\\.{3}\\d+)$")))
Age2 <- suppressWarnings(as.numeric(pick_first(meta, "Age At Onset", "^Age At Onset(\\.|\\.{3}\\d+)$")))
meta$Age <- ifelse(!is.na(Age1), Age1, Age2)

# NEW: Hypertension & Diabetes (robust to dup headers/variants)
meta$Hypertension <- pick_first(meta, "Hypertension", "^Hypertension(\\.|\\.{3}\\d+)")
meta$Diabetes     <- pick_first(meta, "Diabetes",     "^Diabetes(\\.|\\.{3}\\d+)")
# normalize typical encodings (Yes/No, 1/0, Y/N, True/False)
normalize_bin <- function(x){
  y   <- trimws(as.character(x))
  yl  <- tolower(y)

  out <- ifelse(yl %in% c("yes","y","true","1"), "Yes",
         ifelse(yl %in% c("no","n","false","0"),  "No", NA_character_))

  # keep any other non-empty labels as-is (e.g., "Unknown")
  keep <- is.na(out) & !is.na(y) & y != ""
  out[keep] <- y[keep]

  factor(out, levels = c("No","Yes"))
}
if (!all(is.na(meta$Hypertension))) meta$Hypertension <- normalize_bin(meta$Hypertension)
if (!all(is.na(meta$Diabetes)))     meta$Diabetes     <- normalize_bin(meta$Diabetes)

# ---------------- build expression matrix ----------------
sample_idx <- which(startsWith(colnames(expr), "UASG-"))
if (length(sample_idx) == 0) stop("No columns starting with 'UASG-' in expression.")
X <- expr[, sample_idx, drop = FALSE] %>% as.data.frame()
X[] <- lapply(X, function(col) suppressWarnings(as.numeric(col)))

# Align to metadata sample order
common_ids <- intersect(colnames(X), meta$UASG)
if (length(common_ids) < 2) stop("Fewer than 2 overlapping UASG IDs.")
meta_sub <- meta %>% filter(UASG %in% common_ids)
X <- X[, meta_sub$UASG, drop = FALSE]   # reorder columns

# Log if needed
xmax <- suppressWarnings(max(X, na.rm = TRUE))
did_log <- FALSE
if (is.finite(xmax) && xmax > 100) {
  X <- log2(X + 1)
  did_log <- TRUE
  cat("[I] Applied log2(x+1) (max was >", xmax, ").\n")
}

# Remove zero/NA variance rows
var_ok <- apply(X, 1, function(r) { v <- stats::var(r, na.rm = TRUE); is.finite(v) && v > 0 })
Xf <- X[var_ok, , drop = FALSE]
cat("[I] Kept non-zero-variance features:", sum(var_ok), "of", length(var_ok), "\n")

# Fill remaining NA by column median
for (j in seq_len(ncol(Xf))) {
  col <- Xf[[j]]
  if (anyNA(col)) {
    med <- stats::median(col, na.rm = TRUE)
    col[is.na(col)] <- med
    Xf[[j]] <- col
  }
}

# ---------------- PCA ----------------
Xt  <- t(as.matrix(Xf))                 # samples x features
pca <- prcomp(Xt, center = TRUE, scale. = TRUE)

# Variance explained
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
var_df  <- tibble(PC = paste0("PC", seq_along(var_exp)),
                  VarianceExplained = var_exp,
                  Percent = round(100 * var_exp, 2))
write_csv(var_df, file.path(OUTDIR, "pca_variance.csv"))

# Scores (samples) + phenotypes
scores <- as.data.frame(pca$x)
scores$UASG <- rownames(scores)
plot_df <- scores %>%
  left_join(meta_sub %>% select(UASG, StrokeStatus, Sex, Age, Hypertension, Diabetes),
            by = "UASG")

write_csv(plot_df, file.path(OUTDIR, "pca_scores_with_pheno.csv"))

# ---------------- plots ----------------
subtitle_txt <- if (did_log) "log2(x+1), centered & scaled" else "centered & scaled"
xlab <- paste0("PC1 (", round(100 * var_exp[1], 2), "%)")
ylab <- paste0("PC2 (", round(100 * var_exp[2], 2), "%)")

save_scatter <- function(df, aes_col, title, file){
  p <- ggplot(df, aes(PC1, PC2, color = {{aes_col}})) +
    geom_point(size = 2, alpha = 0.9) +
    labs(title = title, subtitle = subtitle_txt, x = xlab, y = ylab) +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTDIR, file), p, width = 7, height = 5, dpi = 150)
}

save_scatter(plot_df, StrokeStatus, "PCA: PC1 vs PC2 by StrokeStatus", "PCA_PC1_PC2_by_StrokeStatus.png")
save_scatter(plot_df, Sex,          "PCA: PC1 vs PC2 by Sex",          "PCA_PC1_PC2_by_Sex.png")

# Age gradient if present
if (any(!is.na(plot_df$Age))) {
  p_age <- ggplot(plot_df, aes(PC1, PC2, color = Age)) +
    geom_point(size = 2, alpha = 0.9) +
    labs(title = "PCA: PC1 vs PC2 by Age", subtitle = subtitle_txt, x = xlab, y = ylab) +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTDIR, "PCA_PC1_PC2_by_Age.png"), p_age, width = 7, height = 5, dpi = 150)
}

# NEW: Hypertension / Diabetes plots (skip gracefully if all NA)
if (any(!is.na(plot_df$Hypertension))) {
  save_scatter(plot_df, Hypertension, "PCA: PC1 vs PC2 by Hypertension", "PCA_PC1_PC2_by_Hypertension.png")
} else {
  message("[W] Hypertension all NA — skipping plot.")
}
if (any(!is.na(plot_df$Diabetes))) {
  save_scatter(plot_df, Diabetes, "PCA: PC1 vs PC2 by Diabetes", "PCA_PC1_PC2_by_Diabetes.png")
} else {
  message("[W] Diabetes all NA — skipping plot.")
}

# Scree
scree <- ggplot(var_df[1:min(20, nrow(var_df)), ], aes(x = PC, y = Percent, group = 1)) +
  geom_col() + geom_point() + geom_line() +
  labs(title = "Scree: percent variance explained", y = "Percent", x = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTDIR, "PCA_scree.png"), scree, width = 8, height = 4.5, dpi = 150)

cat("\nDone. Wrote:\n",
    "- pca_variance.csv\n",
    "- pca_scores_with_pheno.csv\n",
    "- PCA_PC1_PC2_by_StrokeStatus.png\n",
    "- PCA_PC1_PC2_by_Sex.png\n",
    "- PCA_PC1_PC2_by_Age.png (if Age present)\n",
    "- PCA_PC1_PC2_by_Hypertension.png (if available)\n",
    "- PCA_PC1_PC2_by_Diabetes.png (if available)\n",
    "- PCA_scree.png\n", sep = "")


#--------- Make an output folder ---------
# mkdir -p "$OUT_DIR/expr/pca_pheno"

#--------- Make it executable ---------
# chmod +x scripts/expr/expr_pca_by_pheno.R

#--------- Run it ---------
# Rscript scripts/expr/expr_pca_by_pheno.R

#--------- Output files (all in $OUT_DIR/expr/pca_pheno):
# pca_variance.csv
# pca_scores_with_pheno.csv
# PCA_PC1_PC2_by_StrokeStatus.png
# PCA_PC1_PC2_by_Sex.png
# PCA_PC1_PC2_by_Age.png
# PCA_scree.png