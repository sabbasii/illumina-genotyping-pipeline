#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
})

# ================================
# Paths from environment
# ================================
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
PLINK_DIR    <- Sys.getenv("PLINK_DIR")
PHENO_DIR    <- Sys.getenv("PHENO_DIR")
stopifnot(nzchar(EXPR_OUT_DIR), nzchar(PLINK_DIR), nzchar(PHENO_DIR))

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")  # samples start at col 8
geno_tsv <- file.path(PLINK_DIR,    "tmp", "rs966221_genotypes.tsv")  # columns: IID ADD_G DOM_G
map_tsv  <- file.path(PHENO_DIR,    "iid_to_uasg.tsv")                # columns: IID UASG (with header)
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")              # has "Age At Onset" + "Final Diagnosis"

out_dir  <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================
# Load expression (samples from col 8)
# ================================
message("[load] ", expr_csv)
expr_df <- readr::read_csv(expr_csv, show_col_types = FALSE)
stopifnot(ncol(expr_df) >= 9)

colnames(expr_df)[1] <- "ProbeID"
sample_cols <- colnames(expr_df)[8:ncol(expr_df)]
stopifnot(length(sample_cols) > 0)

# Detect GeneSymbol column (best-effort)
lower_names <- tolower(names(expr_df))
gs_idx <- which(grepl("(gene.*symbol|hgnc.*symbol|^symbol$)", lower_names))
GeneSymbol_col <- if (length(gs_idx) > 0) names(expr_df)[gs_idx[1]] else NA_character_

if (!is.na(GeneSymbol_col)) {
  annot_df <- expr_df |>
    dplyr::select(ProbeID, !!GeneSymbol_col) |>
    dplyr::rename(GeneSymbol = !!GeneSymbol_col) |>
    dplyr::mutate(GeneSymbol = as.character(GeneSymbol))
  message("[info] Using GeneSymbol column: ", GeneSymbol_col)
} else {
  annot_df <- expr_df |>
    dplyr::select(ProbeID) |>
    dplyr::mutate(GeneSymbol = NA_character_)
  message("[warn] No gene symbol-like column detected; GeneSymbol will be NA.")
}

# Coerce sample columns to numeric
expr_num <- expr_df
for (cn in sample_cols) expr_num[[cn]] <- suppressWarnings(as.numeric(expr_num[[cn]]))

all_na_cols <- sample_cols[sapply(sample_cols, function(cn) all(is.na(expr_num[[cn]])))]
if (length(all_na_cols) == length(sample_cols)) {
  stop("All selected sample columns are non-numeric. Check that samples start at col 8 and values are numeric.")
}
if (length(all_na_cols) > 0) {
  warning("Dropping non-numeric sample columns: ", paste(all_na_cols, collapse = ", "))
  sample_cols <- setdiff(sample_cols, all_na_cols)
}

expr_mat <- as.matrix(expr_num[, sample_cols, drop = FALSE])
rownames(expr_mat) <- expr_num$ProbeID

# Optional log2 guard for raw-like input
rng <- range(expr_mat, finite = TRUE)
if (is.finite(rng[2]) && rng[2] > 100) {
  message("[info] Expression appears unlogged; applying log2(x+1).")
  expr_mat <- log2(expr_mat + 1)
}

# ================================
# Load genotypes + IID→UASG map
# ================================
message("[load] ", geno_tsv)
geno_df <- readr::read_tsv(geno_tsv, show_col_types = FALSE)
if (!"IID" %in% names(geno_df)) {
  if ("Sample" %in% names(geno_df)) names(geno_df)[names(geno_df)=="Sample"] <- "IID"
  if ("SAMPLE" %in% names(geno_df)) names(geno_df)[names(geno_df)=="SAMPLE"] <- "IID"
}
stopifnot(all(c("IID","ADD_G","DOM_G") %in% names(geno_df)))
geno_df <- geno_df |> dplyr::mutate(IID = trimws(IID))

message("[load] ", map_tsv)
map_df <- readr::read_tsv(map_tsv, show_col_types = FALSE) |>
  dplyr::mutate(IID = trimws(IID), UASG = trimws(UASG))
stopifnot(all(c("IID","UASG") %in% names(map_df)))

# Attach UASG to genotypes
geno_df <- geno_df |> dplyr::left_join(map_df, by = "IID")
if (!"UASG" %in% names(geno_df) || anyNA(geno_df$UASG)) {
  n_miss <- sum(is.na(geno_df$UASG))
  stop(sprintf("Failed to map %d genotype rows to UASG via iid_to_uasg.tsv. Check IID values.", n_miss))
}

# ================================
# Load meta; repair dup names; pick FIRST occurrences
# ================================
message("[load] ", meta_csv)
meta_raw <- read.csv(meta_csv, check.names = FALSE)

orig_names <- names(meta_raw)
names(meta_raw) <- make.unique(orig_names, sep = "__dup")  # UASG, UASG__dup1, ...

# helper: first repaired name for an original header
first_by_original <- function(original_names_vec, repaired_names_vec, target_exact) {
  idx <- which(original_names_vec == target_exact)
  if (length(idx) == 0) return(NA_character_) else return(repaired_names_vec[idx[1]])
}

# UASG: allow variants like UASG...1 / UASG.1 in ORIGINAL names
uasg_idx <- grep("^UASG(\\b|[.]{3}[0-9]+|\\.[0-9]+)?$", orig_names)
uasg_name <- if (length(uasg_idx)) names(meta_raw)[uasg_idx[1]] else NA_character_
sex_name  <- first_by_original(orig_names, names(meta_raw), "Sex")
age_name  <- first_by_original(orig_names, names(meta_raw), "Age At Onset")
diag_name <- first_by_original(orig_names, names(meta_raw), "Final Diagnosis")

needed <- c(uasg_name, sex_name, age_name, diag_name)
if (anyNA(needed)) {
  stop("[ERR] Could not locate required columns. Found -> ",
       "UASG=", uasg_name, " | Sex=", sex_name,
       " | Age At Onset=", age_name, " | Final Diagnosis=", diag_name,
       ". Check header spellings in meta_selected.csv")
}
message("[info] Using columns -> UASG=", uasg_name,
        " | Sex=", sex_name,
        " | Age At Onset=", age_name,
        " | Final Diagnosis=", diag_name)

# Build covariate table (Sex factor; Age centered; Final Diagnosis collapsed)
meta <- meta_raw |>
  dplyr::transmute(
    UASG = trimws(.data[[uasg_name]]),
    Sex  = factor(.data[[sex_name]]),
    AgeAtOnset_raw = suppressWarnings(as.numeric(as.character(.data[[age_name]]))),
    FinalDiagnosis_raw = as.character(.data[[diag_name]])
  ) |>
  dplyr::mutate(
    FinalDiagnosis2 = ifelse(FinalDiagnosis_raw == "Control", "Control", "Stroke"),
    FinalDiagnosis2 = factor(FinalDiagnosis2),
    FinalDiagnosis2 = stats::relevel(FinalDiagnosis2, ref = "Control"),
    AgeCentered = AgeAtOnset_raw - mean(AgeAtOnset_raw, na.rm = TRUE)
  ) |>
  dplyr::select(UASG, Sex, AgeCentered, FinalDiagnosis2)

# ================================
# Align expr, geno, covariates
# ================================
uasg_overlap <- Reduce(intersect, list(colnames(expr_mat), unique(geno_df$UASG), meta$UASG))
if (length(uasg_overlap) < 20) {
  warning("Small overlap across expr/genotype/covariates: ", length(uasg_overlap))
}

# order by common UASG
expr_mat_sub <- expr_mat[, uasg_overlap, drop = FALSE]
geno_sub <- geno_df |>
  dplyr::filter(UASG %in% uasg_overlap) |>
  dplyr::select(UASG, ADD_G, DOM_G) |>
  dplyr::distinct() |>
  dplyr::arrange(match(UASG, uasg_overlap))
cov_sub <- meta |>
  dplyr::filter(UASG %in% uasg_overlap) |>
  dplyr::arrange(match(UASG, uasg_overlap))

stopifnot(identical(colnames(expr_mat_sub), geno_sub$UASG))
stopifnot(identical(geno_sub$UASG, cov_sub$UASG))

# Genotype vectors
ADD <- as.numeric(geno_sub$ADD_G)     # 0/1/2
DOM <- as.numeric(geno_sub$DOM_G)     # 0/1

# Remove samples with NA in needed covariates/genotype
keep_add <- !is.na(ADD) & !is.na(cov_sub$Sex) & !is.na(cov_sub$AgeCentered) & !is.na(cov_sub$FinalDiagnosis2)
keep_dom <- !is.na(DOM) & !is.na(cov_sub$Sex) & !is.na(cov_sub$AgeCentered) & !is.na(cov_sub$FinalDiagnosis2)

X_add <- expr_mat_sub[, keep_add, drop = FALSE]
X_dom <- expr_mat_sub[, keep_dom, drop = FALSE]

cov_add <- cov_sub[keep_add, , drop = FALSE]
cov_dom <- cov_sub[keep_dom, , drop = FALSE]

ADD_k <- ADD[keep_add]
DOM_k <- DOM[keep_dom]

# ================================
# Design matrices (adjusted models)
# Expression ~ 1 + Genotype + Sex + AgeCentered + FinalDiagnosis2
# ================================
design_ADD <- model.matrix(~ ADD_k + cov_add$Sex + cov_add$AgeCentered + cov_add$FinalDiagnosis2)
colnames(design_ADD)[2] <- "ADD"

design_DOM <- model.matrix(~ DOM_k + cov_dom$Sex + cov_dom$AgeCentered + cov_dom$FinalDiagnosis2)
colnames(design_DOM)[2] <- "DOM"

# ================================
# Fit limma + save
# ================================
fit_and_save <- function(X, design, coef_name, prefix, annot_df) {
  fit <- lmFit(X, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  tt  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt$ProbeID <- rownames(tt)
  tt <- tt |>
    dplyr::left_join(annot_df, by = "ProbeID") |>
    dplyr::relocate(ProbeID, GeneSymbol, .before = 1)
  out <- file.path(out_dir, paste0("limma_", prefix, "_adjcov_results.csv"))
  write.csv(tt, out, row.names = FALSE)
  message("[ok] Wrote: ", out)
  invisible(tt)
}

add_res <- fit_and_save(X_add, design_ADD, "ADD", "ADD_G", annot_df)
dom_res <- fit_and_save(X_dom, design_DOM, "DOM", "DOM_G", annot_df)

# ================================
# Run summary
# ================================
summary_path <- file.path(out_dir, "run_summary_adjcov.txt")
sink(summary_path)
cat("Expression–Genotype association (limma) — adjusted for Sex, Age At Onset (centered), Final Diagnosis (Control vs Stroke)\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Expr samples (total):", ncol(expr_mat), "\n")
cat("Genotype UASGs (total):", nrow(geno_df), "\n")
cat("Overlap (expr ∩ geno ∩ cov):", length(uasg_overlap), "\n")
cat("ADD model samples used:", ncol(X_add), "\n")
cat("DOM model samples used:", ncol(X_dom), "\n")
cat("Genes tested (rows):", nrow(expr_mat_sub), "\n\n")
cat("Model ADD_G: Expression ~ 1 + ADD + Sex + AgeCentered + FinalDiagnosis2\n")
cat("Model DOM_G: Expression ~ 1 + DOM + Sex + AgeCentered + FinalDiagnosis2\n\n")
cat("ADD_G FDR<0.05 hits:", sum(add_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("DOM_G FDR<0.05 hits:", sum(dom_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("ADD_G FDR<0.10 hits:", sum(add_res$adj.P.Val < 0.10, na.rm = TRUE), "\n")
cat("DOM_G FDR<0.10 hits:", sum(dom_res$adj.P.Val < 0.10, na.rm = TRUE), "\n")
cat("ADD_G p<0.01 (nominal):", sum(add_res$P.Value < 0.01, na.rm = TRUE), " | p<0.001:", sum(add_res$P.Value < 0.001, na.rm = TRUE), "\n")
cat("DOM_G p<0.01 (nominal):", sum(dom_res$P.Value < 0.01, na.rm = TRUE), " | p<0.001:", sum(dom_res$P.Value < 0.001, na.rm = TRUE), "\n")
sink()
message("[ok] Wrote: ", summary_path)
message("[done] Adjusted association run complete.")

# ------ Run it ------
## Rscript scripts/expr_assoc_rs966221_limma_adjcov.R

# ------ Outputs ($EXPR_OUT_DIR/assoc_rs966221/) ------
# /home/sima/git_projects/illumina-genotyping-pipeline/output/genotype_run1/expr/explore/assoc_rs966221/
  ## limma_ADD_G_adjcov_results.csv
  ## limma_DOM_G_adjcov_results.csv
  ## run_summary_adjcov.txt