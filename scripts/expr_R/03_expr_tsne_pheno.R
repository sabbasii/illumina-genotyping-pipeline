#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr)
  library(stringr); library(tibble)
})

# ---- deps ----
if (!requireNamespace("Rtsne", quietly = TRUE)) {
  stop("Package 'Rtsne' not installed. In R: install.packages('Rtsne')")
}

# ---- paths ----
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
OUT_DIR      <- Sys.getenv("OUT_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(OUT_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")
OUTDIR   <- file.path(OUT_DIR, "expr", "tsne_pheno")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("[I] expr:", expr_csv, "\n[I] meta:", meta_csv, "\n[I] out :", OUTDIR, "\n\n")
stopifnot(file.exists(expr_csv), file.exists(meta_csv))

# ---- helpers ----
pick_first <- function(df, exact, regex=NULL, fallback=NA_character_) {
  idx <- integer(0)
  if (!is.null(exact)) idx <- c(idx, which(names(df) == exact))
  if (!is.null(regex)) idx <- c(idx, grep(regex, names(df)))
  idx <- unique(idx)
  if (length(idx) >= 1) df[[idx[1]]] else rep(fallback, nrow(df))
}
normalize_bin <- function(x){
  y <- trimws(as.character(x))
  y <- dplyr::case_when(
    tolower(y) %in% c("yes","y","true","1")  ~ "Yes",
    tolower(y) %in% c("no","n","false","0") ~ "No",
    TRUE ~ ifelse(is.na(y) | y=="", NA_character_, y)
  )
  factor(y, levels = c("No","Yes"))
}

# ---- load ----
expr <- readr::read_csv(expr_csv, guess_max = 200000, show_col_types = FALSE)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

# Canonical metadata fields
meta$UASG <- pick_first(meta, "UASG", "^UASG(\\.|\\.{3}\\d+)")
stopifnot(any(!is.na(meta$UASG)))

FD <- pick_first(meta, "Final Diagnosis", "^Final Diagnosis(\\.|\\.{3}\\d+)")
if (!"StrokeStatus" %in% names(meta)) {
  meta$StrokeStatus <- ifelse(!is.na(FD) & tolower(trimws(FD))=="control", "Control", "Non-control")
}
meta$Sex          <- pick_first(meta, "Sex", "^Sex(\\.|\\.{3}\\d+)")
Age1 <- suppressWarnings(as.numeric(pick_first(meta, "Age", "^Age(\\.|\\.{3}\\d+)$")))
Age2 <- suppressWarnings(as.numeric(pick_first(meta, "Age At Onset", "^Age At Onset(\\.|\\.{3}\\d+)$")))
meta$Age          <- ifelse(!is.na(Age1), Age1, Age2)
meta$Hypertension <- normalize_bin(pick_first(meta, "Hypertension", "^Hypertension(\\.|\\.{3}\\d+)"))
meta$Diabetes     <- normalize_bin(pick_first(meta, "Diabetes", "^Diabetes(\\.|\\.{3}\\d+)"))
meta$ScanDate     <- pick_first(meta, "Scan Date", "^Scan Date(\\.|\\.{3}\\d+)")

# Expression matrix (samples = columns starting with UASG-)
sample_idx <- which(startsWith(colnames(expr), "UASG-"))
stopifnot(length(sample_idx) > 0)
X <- as.data.frame(expr[, sample_idx, drop = FALSE])
X[] <- lapply(X, function(z) suppressWarnings(as.numeric(z)))

# Align columns to metadata order
common_ids <- intersect(colnames(X), meta$UASG)
stopifnot(length(common_ids) >= 2)
meta_sub <- meta %>% filter(UASG %in% common_ids)
X <- X[, meta_sub$UASG, drop = FALSE]

# Log if needed, remove zero-var features, median-impute NA
xmax <- suppressWarnings(max(X, na.rm = TRUE)); did_log <- FALSE
if (is.finite(xmax) && xmax > 100) { X <- log2(X + 1); did_log <- TRUE }
var_ok <- apply(X, 1, function(r){ v <- stats::var(r, na.rm = TRUE); is.finite(v) && v > 0 })
Xf <- X[var_ok, , drop = FALSE]
for (j in seq_len(ncol(Xf))) {
  col <- Xf[[j]]; if (anyNA(col)) { m <- stats::median(col, na.rm = TRUE); col[is.na(col)] <- m; Xf[[j]] <- col }
}

# ---- t-SNE ----
set.seed(1234)
Xt <- t(as.matrix(Xf))                 # samples x features
n  <- nrow(Xt)

# heuristic perplexity in a safe range
px <- max(5, min(50, floor((n - 1) / 3)))
cat("[I] n samples:", n, "  perplexity:", px, "\n")
emb <- Rtsne::Rtsne(
  Xt, perplexity = px, theta = 0.5, dims = 2, pca = TRUE, verbose = TRUE,
  max_iter = 1000, check_duplicates = FALSE
)

tsne_df <- tibble::tibble(
  UASG = rownames(Xt),
  tSNE1 = emb$Y[,1],
  tSNE2 = emb$Y[,2]
) %>% left_join(meta_sub %>% select(UASG, StrokeStatus, Sex, Age, Hypertension, Diabetes, ScanDate),
                by = "UASG")
readr::write_csv(tsne_df, file.path(OUTDIR, "tsne_embeddings_with_pheno.csv"))

# ---- plots ----
subtitle_txt <- paste0("perplexity=", px, if (did_log) "; log2(x+1), centered & scaled" else "; centered & scaled")

plot_cat <- function(df, col, title, file){
  p <- ggplot(df, aes(tSNE1, tSNE2, color = {{col}})) +
    geom_point(size = 2, alpha = 0.9) +
    labs(title = title, subtitle = subtitle_txt, x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTDIR, file), p, width = 7, height = 5, dpi = 150)
}

plot_cat(tsne_df, StrokeStatus, "t-SNE by StrokeStatus", "tSNE_by_StrokeStatus.png")
plot_cat(tsne_df, Sex,          "t-SNE by Sex",          "tSNE_by_Sex.png")
if (any(!is.na(tsne_df$Hypertension))) plot_cat(tsne_df, Hypertension, "t-SNE by Hypertension", "tSNE_by_Hypertension.png")
if (any(!is.na(tsne_df$Diabetes)))     plot_cat(tsne_df, Diabetes,     "t-SNE by Diabetes",     "tSNE_by_Diabetes.png")
if (any(!is.na(tsne_df$ScanDate)))     plot_cat(tsne_df, ScanDate,     "t-SNE by Scan Date",    "tSNE_by_ScanDate.png")

# Age gradient if available
if (any(!is.na(tsne_df$Age))) {
  p_age <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Age)) +
    geom_point(size = 2, alpha = 0.9) +
    labs(title = "t-SNE by Age", subtitle = subtitle_txt, x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUTDIR, "tSNE_by_Age.png"), p_age, width = 7, height = 5, dpi = 150)
}

# save params
params <- tibble(n_samples = n, perplexity = px, theta = 0.5, dims = 2, max_iter = 1000, did_log = did_log)
write_csv(params, file.path(OUTDIR, "tsne_params.csv"))

cat("\nDone. Wrote:\n",
    "- tsne_embeddings_with_pheno.csv\n",
    "- tSNE_by_StrokeStatus.png\n",
    "- tSNE_by_Sex.png\n",
    "- tSNE_by_Age.png (if Age present)\n",
    "- tSNE_by_Hypertension.png / tSNE_by_Diabetes.png (if present)\n",
    "- tSNE_by_ScanDate.png (if present)\n",
    "- tsne_params.csv\n", sep = "")

# Required libraries
## in a regular terminal (or VS Code terminal): mamba install -c conda-forge r-rtsne

#--------- Make an output folder ---------
# mkdir -p "$OUT_DIR/expr/tsne_pheno"

#--------- Make it executable ---------
# chmod +x scripts/expr/expr_tsne_pheno.R

#--------- Run it ---------
# Rscript scripts/expr/expr_tsne_pheno.R

#--------- Output files (all in $OUT_DIR/expr/tsne_pheno):
# tsne_embeddings_with_pheno.csv
# tSNE_by_StrokeStatus.png
# tSNE_by_Sex.png
# tSNE_by_Age.png
# tSNE_by_Hypertension.png
# tSNE_by_Diabetes.png
# tSNE_by_ScanDate.png
# tsne_params.csv