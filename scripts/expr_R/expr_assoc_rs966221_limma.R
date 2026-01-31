#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
})

# ---------- paths from environment (set by scripts/00_config.sh) ----------
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
PLINK_DIR    <- Sys.getenv("PLINK_DIR")
PHENO_DIR    <- Sys.getenv("PHENO_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(PLINK_DIR), nzchar(PHENO_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.csv")  # mRS row removed
geno_tsv <- file.path(PLINK_DIR,    "tmp", "rs966221_genotypes.tsv") # header: IID ADD_G DOM_G
map_tsv  <- file.path(PHENO_DIR,    "iid_to_uasg.tsv")              # header: IID UASG

out_dir  <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- load inputs ----------
message("[load] ", expr_csv)
expr_df <- read_csv(expr_csv, show_col_types = FALSE)

message("[load] ", geno_tsv)
geno_df <- read_tsv(geno_tsv, show_col_types = FALSE)

message("[load] ", map_tsv)
map_df  <- read_tsv(map_tsv, show_col_types = FALSE)

# ---------- basic checks / normalize headers ----------
stopifnot(ncol(expr_df) >= 9)        # expect ≥8th col to be first sample
stopifnot(all(c("IID","UASG") %in% names(map_df)))
if (!"IID" %in% names(geno_df)) {
  if ("Sample" %in% names(geno_df)) names(geno_df)[names(geno_df)=="Sample"] <- "IID"
  if ("SAMPLE" %in% names(geno_df)) names(geno_df)[names(geno_df)=="SAMPLE"] <- "IID"
}
stopifnot(all(c("IID","ADD_G","DOM_G") %in% names(geno_df)))

map_df  <- map_df  |> mutate(IID = trimws(IID), UASG = trimws(UASG))
geno_df <- geno_df |> mutate(IID = trimws(IID))

# ---------- build UASG ↔ genotype table ----------
geno_df <- geno_df |> left_join(map_df, by = "IID")
if (!"UASG" %in% names(geno_df) || anyNA(geno_df$UASG)) {
  n_miss <- sum(is.na(geno_df$UASG))
  stop(sprintf("Failed to map %d genotype rows to UASG via iid_to_uasg.tsv. Check IID values.", n_miss))
}

# ---------- detect Gene Symbol column from expr_df ----------
colnames(expr_df)[1] <- "ProbeID"

lower_names <- tolower(names(expr_df))
gs_idx <- which(grepl("(gene.*symbol|hgnc.*symbol|^symbol$)", lower_names))
GeneSymbol_col <- if (length(gs_idx) > 0) names(expr_df)[gs_idx[1]] else NA_character_

if (!is.na(GeneSymbol_col)) {
  annot_df <- expr_df |>
    select(ProbeID, !!GeneSymbol_col) |>
    rename(GeneSymbol = !!GeneSymbol_col) |>
    mutate(GeneSymbol = as.character(GeneSymbol))
  message("[info] Using GeneSymbol column: ", GeneSymbol_col)
} else {
  annot_df <- expr_df |>
    select(ProbeID) |>
    mutate(GeneSymbol = NA_character_)
  message("[warn] No gene symbol-like column detected; GeneSymbol will be NA. (You can later post-annotate to HGNC/Ensembl.)")
}

# ---------- prepare expression matrix ----------
# Structure: 1=ProbeID, 2..7 annotations, 8..end = UASG-#### samples
sample_cols <- colnames(expr_df)[8:ncol(expr_df)]           # UASG sample headers
stopifnot(length(sample_cols) > 0)

# Coerce ONLY sample columns to numeric
expr_num <- expr_df
for (cn in sample_cols) {
  expr_num[[cn]] <- suppressWarnings(as.numeric(expr_num[[cn]]))
}
all_na_cols <- sample_cols[sapply(sample_cols, function(cn) all(is.na(expr_num[[cn]])))]
if (length(all_na_cols) == length(sample_cols)) {
  stop("All selected sample columns are non-numeric. Check that sample columns start at col 8 and values are numeric.")
}
if (length(all_na_cols) > 0) {
  warning("Dropping non-numeric sample columns: ", paste(all_na_cols, collapse = ", "))
  sample_cols <- setdiff(sample_cols, all_na_cols)
}

expr_mat <- as.matrix(expr_num[, sample_cols, drop = FALSE])
rownames(expr_mat) <- expr_num$ProbeID

# Optional log2 guard for raw-like scales
rng <- range(expr_mat, finite = TRUE)
if (is.finite(rng[2]) && rng[2] > 100) {
  message("[info] Expression appears unlogged; applying log2(x+1).")
  expr_mat <- log2(expr_mat + 1)
}

# ---------- align samples by UASG ----------
uasg_overlap <- intersect(colnames(expr_mat), unique(geno_df$UASG))
if (length(uasg_overlap) < 20) {
  warning("Small overlap between expression and genotype UASGs: ", length(uasg_overlap))
}
expr_mat_sub <- expr_mat[, uasg_overlap, drop = FALSE]
geno_sub <- geno_df |>
  filter(UASG %in% uasg_overlap) |>
  select(UASG, ADD_G, DOM_G) |>
  distinct() |>
  arrange(UASG)
expr_mat_sub <- expr_mat_sub[, geno_sub$UASG, drop = FALSE]
stopifnot(identical(colnames(expr_mat_sub), geno_sub$UASG))

# ---------- design matrices ----------
ADD <- as.numeric(geno_sub$ADD_G)          # 0/1/2 dosage of risk allele G
DOM <- as.numeric(geno_sub$DOM_G)          # 0 vs 1 (AG/GG)

keep_add <- !is.na(ADD)
keep_dom <- !is.na(DOM)
X_add <- expr_mat_sub[, keep_add, drop = FALSE];  ADD <- ADD[keep_add]
X_dom <- expr_mat_sub[, keep_dom, drop = FALSE];  DOM <- DOM[keep_dom]

design_ADD <- model.matrix(~ ADD)   # Intercept + ADD
design_DOM <- model.matrix(~ DOM)   # Intercept + DOM

# ---------- fit limma + save (with GeneSymbol merged) ----------
fit_and_save <- function(X, design, prefix, annot_df) {
  fit <- lmFit(X, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "P")  # coef 2 = ADD/DOM
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt$ProbeID <- rownames(tt)
  tt <- tt |>
    dplyr::left_join(annot_df, by = "ProbeID") |>
    dplyr::relocate(ProbeID, GeneSymbol, .before = 1)

  # write full table
  full_out <- file.path(out_dir, paste0("limma_", prefix, "_results.csv"))
  write.csv(tt, full_out, row.names = FALSE)

  # relaxed-FDR and nominal p-value exports
  fdr10 <- dplyr::filter(tt, adj.P.Val < 0.10)
  p001  <- dplyr::filter(tt, P.Value   < 0.001)
  p01   <- dplyr::filter(tt, P.Value   < 0.01)
  top100 <- tt[order(tt$P.Value), ][1:min(100, nrow(tt)), ]

  write.csv(fdr10, file.path(out_dir, paste0("hits_", prefix, "_FDR_lt_0.10.csv")), row.names = FALSE)
  write.csv(p001,  file.path(out_dir, paste0("hits_", prefix, "_P_lt_1e-3.csv")),  row.names = FALSE)
  write.csv(p01,   file.path(out_dir, paste0("hits_", prefix, "_P_lt_1e-2.csv")),  row.names = FALSE)
  write.csv(top100,file.path(out_dir, paste0("top100_", prefix, "_by_P.csv")),     row.names = FALSE)

  message("[ok] Wrote: ", full_out)
  invisible(list(full=tt, fdr10=fdr10, p001=p001, p01=p01, top100=top100))
}

res_add <- fit_and_save(X_add, design_ADD, "ADD_G", annot_df)
res_dom <- fit_and_save(X_dom, design_DOM, "DOM_G", annot_df)

# ---------- run summary (includes nominal p-value counts) ----------
summary_path <- file.path(out_dir, "run_summary.txt")
sink(summary_path)
cat("Expression–Genotype association (limma)\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Expr samples (total):", ncol(expr_mat), "\n")
cat("Genotype UASGs (total):", nrow(geno_df), "\n")
cat("Overlapped UASGs (ADD):", ncol(X_add), "\n")
cat("Overlapped UASGs (DOM):", ncol(X_dom), "\n")
cat("Genes tested (rows):", nrow(expr_mat_sub), "\n\n")
cat("Model ADD_G: Intercept + ADD (0/1/2; risk allele = G)\n")
cat("Model DOM_G: Intercept + DOM (1 if AG/GG else 0)\n\n")

cat("ADD_G FDR<0.10 hits:", nrow(res_add$fdr10), "\n")
cat("ADD_G p<0.01 (nominal):", nrow(res_add$p01), " | p<0.001:", nrow(res_add$p001), "\n")
cat("DOM_G FDR<0.10 hits:", nrow(res_dom$fdr10), "\n")
cat("DOM_G p<0.01 (nominal):", nrow(res_dom$p01), " | p<0.001:", nrow(res_dom$p001), "\n")
sink()
message("[ok] Wrote: ", summary_path)
message("[done] Association run complete.")
