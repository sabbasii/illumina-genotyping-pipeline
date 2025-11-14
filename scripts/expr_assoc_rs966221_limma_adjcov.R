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

colnames(expr_df)[1] <- "ProbeID"  # this is changing the first col name into ProbeID
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
# loop over each column name in sample_cols and, for every one of those columns in expr_num, convert its values to numeric with as.numeric(), silencing any warnings that would normally appear during the conversion.
for (cn in sample_cols) expr_num[[cn]] <- suppressWarnings(as.numeric(expr_num[[cn]]))

# create a 'character vector' of column names that are completely NA. all(is.na(expr_num[[cn]])) → TRUE if every value in that column is NA.
all_na_cols <- sample_cols[sapply(sample_cols, function(cn) all(is.na(expr_num[[cn]])))]
# If all sample columns are entirely NA, the script stops with an error because there is no usable numeric data to analyze.
if (length(all_na_cols) == length(sample_cols)) {
  stop("All selected sample columns are non-numeric. Check that samples start at col 8 and values are numeric.")
}
# If only some sample columns are all NA, the script warns you about those columns and drops them, keeping only the valid numeric sample columns for further analysis.
if (length(all_na_cols) > 0) {
  warning("Dropping non-numeric sample columns: ", paste(all_na_cols, collapse = ", "))
  sample_cols <- setdiff(sample_cols, all_na_cols)
}

# as.matrix(...) → converts that data.frame/tibble subset into a matrix.
# takes only the columns listed in sample_cols from expr_num, and keeps them as a 2D table
expr_mat <- as.matrix(expr_num[, sample_cols, drop = FALSE])
# expr_mat is a numeric matrix of expression values (rows = probes/genes, columns = samples), ready for matrix-based functions (e.g., limma/linear modeling).
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
# Check if the genotype table is missing an IID column; if so, look for a column named Sample or SAMPLE and rename whichever exists to IID so downstream code has a standard ID column name.
if (!"IID" %in% names(geno_df)) {
  if ("Sample" %in% names(geno_df)) names(geno_df)[names(geno_df)=="Sample"] <- "IID"
  if ("SAMPLE" %in% names(geno_df)) names(geno_df)[names(geno_df)=="SAMPLE"] <- "IID"
}
stopifnot(all(c("IID","ADD_G","DOM_G") %in% names(geno_df)))
geno_df <- geno_df |> dplyr::mutate(IID = trimws(IID))

message("[load] ", map_tsv)
map_df <- readr::read_tsv(map_tsv, show_col_types = FALSE) |>
  dplyr::mutate(IID = trimws(IID), UASG = trimws(UASG)) # trims white space and overwrites the old values with the trimmed  ones in the same columns.
stopifnot(all(c("IID","UASG") %in% names(map_df)))

# Attach UASG to genotypes
geno_df <- geno_df |> dplyr::left_join(map_df, by = "IID")
if (!"UASG" %in% names(geno_df) || anyNA(geno_df$UASG)) {
  n_miss <- sum(is.na(geno_df$UASG))
  stop(sprintf("Failed to map %d genotype rows to UASG via iid_to_uasg.tsv. Check IID values.", n_miss))
}
# ==========sideNote==========
# A 'left join' keeps all rows from the left table (geno_df) 
# and adds matching columns from the right table (map_df) based on the key (IID).
# Rows in geno_df with no match in map_df get NA for the new columns.

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
    # Recode diagnosis into two groups: "Control" stays "Control", everything else becomes "Stroke".
    FinalDiagnosis2 = ifelse(FinalDiagnosis_raw == "Control", "Control", "Stroke"),
    # Turn that new variable into a categorical factor
    FinalDiagnosis2 = factor(FinalDiagnosis2),
    # Set "Control" as the reference level (baseline) for modeling
    FinalDiagnosis2 = stats::relevel(FinalDiagnosis2, ref = "Control"),
    # subtract the 'mean age' so the age variable is centered around zero for better model stability and interpretation.
    AgeCentered = AgeAtOnset_raw - mean(AgeAtOnset_raw, na.rm = TRUE)
  ) |>
  dplyr::select(UASG, Sex, AgeCentered, FinalDiagnosis2)

# ==========sideNote==========
# transmute creates only these new columns (drops everything else)
# .data → “the current data frame inside this mutate/transmute call”
# (here: meta_raw)


# ================================
# Align expr, geno, covariates
# ================================
uasg_overlap <- Reduce(
  intersect,
  list(colnames(expr_mat), unique(geno_df$UASG), meta$UASG)
)

if (length(uasg_overlap) < 20) {
  warning("Small overlap across expr/genotype/covariates: ",
          length(uasg_overlap))
}
# uasg_overlap = vector of UASG IDs present
# in expression, genotype, and metadata.

# Keep only the columns (samples) in expr_mat
# whose UASG IDs are in uasg_overlap.
expr_mat_sub <- expr_mat[, uasg_overlap, drop = FALSE]

geno_sub <- geno_df |>
  dplyr::filter(UASG %in% uasg_overlap) |>
  dplyr::select(UASG, ADD_G, DOM_G) |>
  dplyr::distinct() |>
  dplyr::arrange(match(UASG, uasg_overlap))
# ==========sideNote==========
# ADD_G → additive genotype (0/1/2)
# DOM_G → dominant genotype (0/1)
# distinct():
# Drop any duplicate rows
# (e.g., if the same UASG appeared twice with the same genotypes).

cov_sub <- meta |>
  # Keep only rows whose UASG is in uasg_overlap
  dplyr::filter(UASG %in% uasg_overlap) |>
  # Sort rows so UASG matches the order in uasg_overlap
  dplyr::arrange(match(UASG, uasg_overlap))

stopifnot(identical(colnames(expr_mat_sub), geno_sub$UASG))
stopifnot(identical(geno_sub$UASG, cov_sub$UASG))

# Genotype vectors
# Pull genotype columns from geno_sub and convert them to numeric vectors
ADD <- as.numeric(geno_sub$ADD_G)     # 0/1/2
DOM <- as.numeric(geno_sub$DOM_G)     # 0/1

# Remove samples with NA in needed covariates/genotype
keep_add <- !is.na(ADD) & !is.na(cov_sub$Sex) & !is.na(cov_sub$AgeCentered) & !is.na(cov_sub$FinalDiagnosis2)
keep_dom <- !is.na(DOM) & !is.na(cov_sub$Sex) & !is.na(cov_sub$AgeCentered) & !is.na(cov_sub$FinalDiagnosis2)

# Subset the expression matrix to samples with complete data:
# X_add: samples valid for the additive model.
# X_dom: samples valid for the dominant model.
X_add <- expr_mat_sub[, keep_add, drop = FALSE] # dim = (21448, 223)
X_dom <- expr_mat_sub[, keep_dom, drop = FALSE] # dim = (21448, 223)

# Subset the covariate table
# X_add should matche cov_add sample-wise, and X_dom should matche cov_dom
cov_add <- cov_sub[keep_add, , drop = FALSE] # dim = (223, 4)
cov_dom <- cov_sub[keep_dom, , drop = FALSE] # dim = (223, 4)

ADD_k <- ADD[keep_add] # length = 223
DOM_k <- DOM[keep_dom] # length = 223

# ================================
# Design matrices (adjusted models)
# Expression ~ 1 + Genotype + Sex + AgeCentered + FinalDiagnosis2
# ================================
design_ADD <- model.matrix(~ ADD_k + cov_add$Sex + cov_add$AgeCentered + cov_add$FinalDiagnosis2)
colnames(design_ADD)[2] <- "ADD" # rename "ADD_k" to "ADD"

design_DOM <- model.matrix(~ DOM_k + cov_dom$Sex + cov_dom$AgeCentered + cov_dom$FinalDiagnosis2)
colnames(design_DOM)[2] <- "DOM" # rename "DOM_k" to "DOM"

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
# ==========sideNote==========
# coef_name tells limma which predictor you want to test.
# design matrix has one column per predictor in the model
# (e.g., ADD, SexMale, AgeCentered, etc.).
# Each column corresponds to a coefficient (a β value) in the regression.
# When you run:
# topTable(fit, coef = "ADD")
# you are asking:
# ➡️ “Show me the genes whose expression is associated with the ADD genotype effect.”
# Coefficient = a predictor’s effect in the model.
# coef_name = which effect (column) you want limma to return results for.

# lmFit : fit the linear model for each gene using the design matrix.
# eBayes : apply empirical Bayes shrinkage (stabilizes variances across genes).
# robust = TRUE, trend = TRUE make it more robust and allow mean–variance trend.
# Extract the per-gene results table (tt) for the coefficient coef_name.
# number = Inf → keep all genes, not just the top 10.
# Sort by p-value.
# Recompute FDR-adjusted p-values using Benjamini–Hochberg
#and store in adj.P.Val.
# Add a ProbeID column from the row names (each row = one probe).
# Join with your annotation table annot_df to attach gene symbols
#(and maybe other info).
# Move ProbeID and GeneSymbol to the front of the table.
# Build an output path like:
# limma_ADD_adjcov_results.csv or limma_DOM_adjcov_results.csv in out_dir


# ================================
add_res <- fit_and_save(X_add, design_ADD, "ADD", "ADD_G", annot_df)
dom_res <- fit_and_save(X_dom, design_DOM, "DOM", "DOM_G", annot_df)

# ================================
# Run summary
# ================================
summary_path <- file.path(out_dir, "run_summary_adjcov.txt")
sink(summary_path)

cat(
  "Expression–Genotype association (limma) — adjusted for Sex, ",
  "Age At Onset (centered), Final Diagnosis (Control vs Stroke)\n",
  sep = ""
)
cat("Date:", format(Sys.time()), "\n\n")
cat("Expr samples (total):", ncol(expr_mat), "\n")
cat("Genotype UASGs (total):", nrow(geno_df), "\n")
cat("Overlap (expr ∩ geno ∩ cov):", length(uasg_overlap), "\n")
cat("ADD model samples used:", ncol(X_add), "\n")
cat("DOM model samples used:", ncol(X_dom), "\n")
cat("Genes tested (rows):", nrow(expr_mat_sub), "\n\n")
cat("Model ADD_G: Expression ~ 1 + ADD + Sex + AgeCentered + FinalDiagnosis2\n")
# ADD
cat(
  "ADD_G FDR<0.05 hits:",
  sum(add_res$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)
cat(
  "ADD_G FDR<0.10 hits:",
  sum(add_res$adj.P.Val < 0.10, na.rm = TRUE),
  "\n"
)
cat(
  "ADD_G p<0.05:",  sum(add_res$P.Value < 0.05,  na.rm = TRUE),
  " | p<0.01:",     sum(add_res$P.Value < 0.01,  na.rm = TRUE),
  " | p<0.001:",    sum(add_res$P.Value < 0.001, na.rm = TRUE),
  "\n\n"
)
# DOM
cat("Model DOM_G: Expression ~ 1 + DOM + Sex + AgeCentered + FinalDiagnosis2\n")
cat(
  "DOM_G FDR<0.05 hits:",
  sum(dom_res$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)
cat(
  "DOM_G FDR<0.10 hits:",
  sum(dom_res$adj.P.Val < 0.10, na.rm = TRUE),
  "\n"
)
cat(
  "DOM_G p<0.05:",  sum(dom_res$P.Value < 0.05,  na.rm = TRUE),
  " | p<0.01:",     sum(dom_res$P.Value < 0.01,  na.rm = TRUE),
  " | p<0.001:",    sum(dom_res$P.Value < 0.001, na.rm = TRUE),
  "\n"
)
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