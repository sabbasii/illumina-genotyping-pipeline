#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(limma); library(ggplot2)
})

# ------------ paths ------------
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
OUT_DIR      <- Sys.getenv("OUT_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(OUT_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")
OUTDIR   <- file.path(OUT_DIR, "expr", "dge_unadjusted") # this builds a path string
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE) # dir.create(path, ) creates the directory

cat("[I] expr:", expr_csv, "\n[I] meta:", meta_csv, "\n[I] out :", OUTDIR, "\n\n")
stopifnot(file.exists(expr_csv), file.exists(meta_csv))

# ------------ helpers ------------
pick_first <- function(df, exact, regex=NULL, fallback=NA_character_) {
  idx <- integer(0)
  if (!is.null(exact)) idx <- c(idx, which(names(df) == exact))
  if (!is.null(regex)) idx <- c(idx, grep(regex, names(df)))
  idx <- unique(idx)
  if (length(idx) >= 1) df[[idx[1]]] else rep_len(fallback, nrow(df))
}

# ------------ load ------------
expr <- readr::read_csv(expr_csv, guess_max = 200000, show_col_types = FALSE)
meta <- readr::read_csv(meta_csv, show_col_types = FALSE)

# ------------ locate ID & symbol columns ------------
feature_id_col <- names(expr)[1]  # first column = ProbeID/GeneID
symbol_candidates <- c("Gene Symbol","GeneSymbol","Gene_Symbol","Symbol","HGNC symbol","HGNC Symbol")
gs_idx <- which(names(expr) %in% symbol_candidates)
if (length(gs_idx) == 0) gs_idx <- grep("symbol", names(expr), ignore.case = TRUE)
gene_symbol_col <- if (length(gs_idx) >= 1) names(expr)[gs_idx[1]] else NA_character_

# ------------ canonicalize metadata ------------
meta$UASG <- pick_first(meta, "UASG", "^UASG(\\.|\\.{3}\\d+)")
if (!"StrokeStatus" %in% names(meta)) {
  FD <- pick_first(meta, "Final Diagnosis", "^Final Diagnosis(\\.|\\.{3}\\d+)")
  meta$StrokeStatus <- ifelse(!is.na(FD) & tolower(trimws(FD))=="control", "Control", "Non-control")
} # this creates a new column named StrokeStatus
meta <- meta %>% filter(!is.na(StrokeStatus), !is.na(UASG))

# Set Group with clean, syntactic level names
raw_levels <- c("Control","Non-control")
meta$Group <- factor(meta$StrokeStatus, levels = raw_levels)
levels(meta$Group) <- make.names(levels(meta$Group))      # "Non.control"
levels(meta$Group) <- gsub("\\.", "_", levels(meta$Group)) # "Non_control"
cat("[I] Group levels:", paste(levels(meta$Group), collapse = ", "), "\n")

# ------------ build expression matrix ------------
sample_idx <- which(startsWith(colnames(expr), "UASG-"))
if (length(sample_idx) == 0) stop("No columns starting with 'UASG-' in the expression file.")
X <- as.data.frame(expr[, sample_idx, drop = FALSE])
X[] <- lapply(X, function(z) suppressWarnings(as.numeric(z)))

common_ids <- intersect(colnames(X), meta$UASG)
if (length(common_ids) < 2) stop("Fewer than 2 overlapping UASG IDs between expression and metadata.")
meta_sub <- meta %>% filter(UASG %in% common_ids)
X <- X[, meta_sub$UASG, drop = FALSE]

# ------------ log2 if needed; filter zero-variance ------------
xmax <- suppressWarnings(max(X, na.rm = TRUE))
did_log <- FALSE
if (is.finite(xmax) && xmax > 100) { X <- log2(X + 1); did_log <- TRUE; cat("[I] Applied log2(x+1).\n") }

var_ok <- apply(X, 1, function(r){ v <- stats::var(r, na.rm = TRUE); is.finite(v) && v > 0 })
Xf <- as.matrix(X[var_ok, , drop = FALSE])
cat("[I] Kept non-zero-variance features:", sum(var_ok), "of", length(var_ok), "\n")

# ------------ limma design & fit ------------
design <- model.matrix(~ 0 + Group, data = meta_sub)  # no intercept
colnames(design) <- levels(meta_sub$Group)            # "Control","Non_control"
cat("[I] Design columns:", paste(colnames(design), collapse = ", "), "\n")

fit <- lmFit(Xf, design)
contrast.matrix <- makeContrasts(Non_control_vs_Control = Non_control - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ------------ results ------------
feature_ids  <- expr[[feature_id_col]][var_ok]
gene_symbols <- if (!is.na(gene_symbol_col)) expr[[gene_symbol_col]][var_ok] else rep(NA_character_, sum(var_ok))

tbl <- topTable(fit2, coef = "Non_control_vs_Control", number = Inf, sort.by = "P")
tbl$FeatureID  <- feature_ids[match(rownames(tbl), rownames(Xf))]
tbl$GeneSymbol <- gene_symbols[match(rownames(tbl), rownames(Xf))]
tbl <- tbl %>% relocate(FeatureID, GeneSymbol, .before = 1)

out_all  <- file.path(OUTDIR, "dge_unadjusted_all.csv")
out_fdr  <- file.path(OUTDIR, "dge_unadjusted_FDR_lt_0.05.csv")
out_p001 <- file.path(OUTDIR, "dge_unadjusted_p_lt_0.001.csv")
write_csv(tbl, out_all)
write_csv(tbl %>% filter(adj.P.Val < 0.05), out_fdr)
write_csv(tbl %>% filter(P.Value   < 0.001), out_p001)
cat("[I] Wrote:\n - ", out_all, "\n - ", out_fdr, "\n - ", out_p001, "\n", sep = "")

# ------------ volcano ------------
vol <- tbl %>% mutate(sig = ifelse(adj.P.Val < 0.05, "FDR<0.05", "NS"),
                      neglog10p = -log10(P.Value))
p <- ggplot(vol, aes(x = logFC, y = neglog10p, color = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("NS"="#bdbdbd","FDR<0.05"="#E24A33")) +
  labs(
    title = "Volcano: Non_control vs Control (unadjusted)",
    subtitle = if (did_log) "log2(x+1), limma eBayes" else "limma eBayes",
    x = "log2 Fold Change (Non_control / Control)",
    y = "-log10(p-value)"
  ) + theme_minimal(base_size = 12)
ggsave(file.path(OUTDIR, "volcano_unadjusted.png"), p, width = 7, height = 5, dpi = 150)

# ------------ MA plot (FIX: use AveExpr, not Amean) ------------
ma_df <- tbl %>% mutate(sig = ifelse(adj.P.Val < 0.05, "FDR<0.05", "NS"))
pma <- ggplot(ma_df, aes(x = AveExpr, y = logFC, color = sig)) +
  geom_point(alpha = 0.7, size = 1.4) +
  scale_color_manual(values = c("NS"="#bdbdbd","FDR<0.05"="#1f77b4")) +
  labs(title = "MA plot: Non_control vs Control (unadjusted)",
       x = "Average log-expression (AveExpr)", y = "log2 Fold Change") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUTDIR, "MAplot_unadjusted.png"), pma, width = 7, height = 5, dpi = 150)

# Reproducibility
writeLines(capture.output(sessionInfo()), con = file.path(OUTDIR, "sessionInfo.txt"))
cat("\n[I] Done.\n")


# ------ Run it ------
## Rscript scripts/expr/05_dge_limma_unadjusted.R
