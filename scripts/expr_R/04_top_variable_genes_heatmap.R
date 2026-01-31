#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
})

# ---- deps ----
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package 'pheatmap' not installed. In R: install.packages('pheatmap')")
}

# ---- paths / config ----
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
OUT_DIR      <- Sys.getenv("OUT_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(OUT_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")
OUTDIR   <- file.path(OUT_DIR, "expr", "topvar_heatmap")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# top-N (arg1 or env TOPN; default 500)
args <- commandArgs(trailingOnly = TRUE)
TOPN <- if (length(args) >= 1) as.integer(args[1]) else
        if (nzchar(Sys.getenv("TOPN"))) as.integer(Sys.getenv("TOPN")) else 500L
if (is.na(TOPN) || TOPN < 50L) TOPN <- 500L

cat("[I] expr:", expr_csv, "\n[I] meta:", meta_csv, "\n[I] out :", OUTDIR, "\n[I] TOPN:", TOPN, "\n\n")
stopifnot(file.exists(expr_csv), file.exists(meta_csv))

# ---- helpers ----
pick_first <- function(df, exact, regex=NULL, fallback=NA_character_) {
  idx <- integer(0)
  if (!is.null(exact)) idx <- c(idx, which(names(df) == exact))
  if (!is.null(regex)) idx <- c(idx, grep(regex, names(df)))
  idx <- unique(idx)
  if (length(idx) >= 1) df[[idx[1]]] else rep_len(fallback, nrow(df))
}
normalize_bin <- function(x){
  y   <- trimws(as.character(x)); yl <- tolower(y)
  out <- ifelse(yl %in% c("yes","y","true","1"), "Yes",
         ifelse(yl %in% c("no","n","false","0"),  "No", NA_character_))
  keep <- is.na(out) & !is.na(y) & y != ""
  out[keep] <- y[keep]
  factor(out, levels = c("No","Yes"))
}
coalesce_str <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  out <- ifelse(!is.na(a) & nzchar(a), a, b)
  out
}

# ---- load ----
expr <- readr::read_csv(expr_csv, guess_max = 200000, show_col_types = FALSE)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

# ---- locate ID & Gene Symbol columns in expression ----
feature_id_col <- names(expr)[1]                 # first column = ProbeID/GeneID
gene_symbol_candidates <- c("Gene Symbol","GeneSymbol","Gene_Symbol","Symbol","HGNC symbol","HGNC Symbol")
gs_idx <- which(names(expr) %in% gene_symbol_candidates)
if (length(gs_idx) == 0) {
  # fallback: try regex for 'symbol'
  gs_idx <- grep("symbol", names(expr), ignore.case = TRUE)
}
gene_symbol_col <- if (length(gs_idx) >= 1) names(expr)[gs_idx[1]] else NA_character_

# ---- canonical metadata columns ----
meta$UASG         <- pick_first(meta, "UASG", "^UASG(\\.|\\.{3}\\d+)")
FD                <- pick_first(meta, "Final Diagnosis", "^Final Diagnosis(\\.|\\.{3}\\d+)")
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

# ---- expression matrix (features x samples) ----
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

# ---- compute variance & select TOPN ----
feature_ids <- expr[[feature_id_col]][var_ok]
gene_symbols_all <- if (!is.na(gene_symbol_col)) expr[[gene_symbol_col]][var_ok] else rep(NA_character_, sum(var_ok))

vars <- apply(Xf, 1, function(r) stats::var(r, na.rm = TRUE))
ord <- order(vars, decreasing = TRUE)
top_idx <- ord[seq_len(min(TOPN, length(ord)))]
X_top <- Xf[top_idx, , drop = FALSE]

top_table <- tibble(
  FeatureID  = feature_ids[top_idx],
  GeneSymbol = gene_symbols_all[top_idx],
  Variance   = vars[top_idx]
)
# write ranked list
readr::write_csv(top_table, file.path(OUTDIR, sprintf("top_variable_genes_top%d.csv", nrow(X_top))))

# also write the numeric matrix with FeatureID + GeneSymbol in first two cols
mat_df <- as.data.frame(X_top)
mat_df <- cbind(FeatureID = feature_ids[top_idx],
                GeneSymbol = gene_symbols_all[top_idx],
                mat_df)
readr::write_csv(mat_df, file.path(OUTDIR, sprintf("top_variable_matrix_top%d.csv", nrow(X_top))))

# ---- z-score rows for heatmap ----
zscore_rows <- function(m){
  m <- as.matrix(m)
  m <- t(scale(t(m), center = TRUE, scale = TRUE))
  m[is.na(m)] <- 0
  m
}
Z <- zscore_rows(X_top)

# row names for plotting: prefer Symbol, else ID; ensure uniqueness
row_names <- coalesce_str(gene_symbols_all[top_idx], feature_ids[top_idx])
rownames(Z) <- make.unique(row_names)
colnames(Z) <- meta_sub$UASG  # already aligned

# ---- column annotations ----
ann <- meta_sub %>% select(UASG, StrokeStatus, Sex, Hypertension, Diabetes, ScanDate)
ann <- as.data.frame(ann); rownames(ann) <- ann$UASG; ann$UASG <- NULL

# ---- colors (heatmap) ----
# negatives ~ yellow/light; positives ~ dark blue (your preference)
cols <- colorRampPalette(c("#F6E58D", "#FFFFFF", "#4472C4"))(101)

# ---- annotation colors (robust) ----
ann_colors <- NULL
if (requireNamespace("RColorBrewer", quietly = TRUE)) {
  get_pal <- function(k) RColorBrewer::brewer.pal(max(3, k), "Set1")[1:k]
  ann_colors <- list()
  map_levels <- function(x) { lv <- levels(factor(x)); setNames(get_pal(length(lv)), lv) }
  if (!all(is.na(ann$StrokeStatus))) ann_colors$StrokeStatus <- map_levels(ann$StrokeStatus)
  if (!all(is.na(ann$Sex)))          ann_colors$Sex          <- map_levels(ann$Sex)
  if (!all(is.na(ann$ScanDate)))     ann_colors$ScanDate     <- map_levels(ann$ScanDate)
  ann_colors$Hypertension <- c(No = "#A1D99B", Yes = "#31A354")
  ann_colors$Diabetes     <- c(No = "#9ECAE1", Yes = "#3182BD")
} else {
  ann_colors <- list(
    Hypertension = c(No = "#A1D99B", Yes = "#31A354"),
    Diabetes     = c(No = "#9ECAE1", Yes = "#3182BD")
  )
}

# ---- draw heatmap (PNG + PDF) ----
pheatmap::pheatmap(
  Z,
  color = cols,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_col = ann,
  annotation_colors = ann_colors,
  main = sprintf("Top variable genes (n=%d)%s", nrow(Z), if (did_log) " — log2(x+1)" else ""),
  filename = file.path(OUTDIR, sprintf("heatmap_topvar_top%d.png", nrow(Z))),
  width = 12, height = 9, dpi = 150
)

pheatmap::pheatmap(
  Z,
  color = cols,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_col = ann,
  annotation_colors = ann_colors,
  main = sprintf("Top variable genes (n=%d)%s", nrow(Z), if (did_log) " — log2(x+1)" else ""),
  filename = file.path(OUTDIR, sprintf("heatmap_topvar_top%d.pdf", nrow(Z))),
  width = 12, height = 9
)

cat("\nWrote:\n",
    "- ", file.path(OUTDIR, sprintf("top_variable_genes_top%d.csv", nrow(Z))), "\n",
    "- ", file.path(OUTDIR, sprintf("top_variable_matrix_top%d.csv", nrow(Z))), "\n",
    "- ", file.path(OUTDIR, sprintf("heatmap_topvar_top%d.png", nrow(Z))), "\n",
    "- ", file.path(OUTDIR, sprintf("heatmap_topvar_top%d.pdf", nrow(Z))), "\n", sep = "")


# ------ RUN ------
# Rscript scripts/expr/top_variable_genes_heatmap.R 500
# or change 500 to another N (e.g., 1000)


# ------ Outputs ($OUT_DIR/expr/topvar_heatmap/) ------
## top_variable_genes_topN.csv (ranked list: feature ID + variance)
## top_variable_matrix_topN.csv (numeric matrix used)
## heatmap_topvar_topN.png and .pdf (z-score, annotated by StrokeStatus/Sex/Hypertension/Diabetes/ScanDate)