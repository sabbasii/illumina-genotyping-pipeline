#!/usr/bin/env Rscript
# Expression–Genotype association for rs966221 with covariates
# Adjusts for Age, Sex, StrokeStatus (no Plate/Run).
# Uses expr_selected_clean.nomrs.csv (UASGs start at column 8).
# In meta_selected.csv (which may have duplicate UASG columns), uses the FIRST UASG-like column.

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

expr_csv <- file.path(EXPR_OUT_DIR, "expr_selected_clean.nomrs.csv")     # NOTE: nomrs
geno_tsv <- file.path(PLINK_DIR,    "tmp", "rs966221_genotypes.tsv")     # header: IID ADD_G DOM_G
map_tsv  <- file.path(PHENO_DIR,    "iid_to_uasg.tsv")                   # header: IID UASG
meta_csv <- file.path(EXPR_OUT_DIR, "meta_selected.csv")

out_dir  <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- load expression ----------
message("[load] ", expr_csv)
expr_df <- read_csv(expr_csv, show_col_types = FALSE)
stopifnot(ncol(expr_df) >= 9)
colnames(expr_df)[1] <- "ProbeID"

# Detect GeneSymbol-like column
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
  message("[warn] No gene symbol-like column detected; GeneSymbol will be NA.")
}

# Samples start at column 8
sample_cols <- colnames(expr_df)[8:ncol(expr_df)]
stopifnot(length(sample_cols) > 0)

# Coerce sample columns to numeric
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

# Optional log2 guard
rng <- range(expr_mat, finite = TRUE)
if (is.finite(rng[2]) && rng[2] > 100) {
  message("[info] Expression appears unlogged; applying log2(x+1).")
  expr_mat <- log2(expr_mat + 1)
}

# ---------- load genotypes + IID→UASG map ----------
message("[load] ", geno_tsv)
geno_df <- read_tsv(geno_tsv, show_col_types = FALSE)
if (!"IID" %in% names(geno_df)) {
  if ("Sample" %in% names(geno_df)) names(geno_df)[names(geno_df)=="Sample"] <- "IID"
  if ("SAMPLE" %in% names(geno_df)) names(geno_df)[names(geno_df)=="SAMPLE"] <- "IID"
}
stopifnot(all(c("IID","ADD_G","DOM_G") %in% names(geno_df)))
geno_df <- geno_df |> mutate(IID = trimws(IID))

message("[load] ", map_tsv)
map_df  <- read_tsv(map_tsv, show_col_types = FALSE) |>
  mutate(IID = trimws(IID), UASG = trimws(UASG))
stopifnot(all(c("IID","UASG") %in% names(map_df)))

# Map to UASG
geno_df <- geno_df |> left_join(map_df, by = "IID")
if (!"UASG" %in% names(geno_df) || anyNA(geno_df$UASG)) {
  n_miss <- sum(is.na(geno_df$UASG))
  stop(sprintf("Failed to map %d genotype rows to UASG via iid_to_uasg.tsv. Check IID values.", n_miss))
}

# ---------- load metadata; pick FIRST UASG-like column ----------
if (!file.exists(meta_csv)) stop("[ERR] Missing metadata: ", meta_csv)
meta_raw <- read.csv(meta_csv, check.names = FALSE)

uasg_idx <- grep("^UASG(\\b|[.]{3}[0-9]+|\\.[0-9]+)?$", names(meta_raw))
if (length(uasg_idx) == 0) stop("[ERR] No UASG-like column found in meta_selected.csv")
uasg_first <- uasg_idx[1]
message("[info] Using UASG column from meta index ", uasg_first, " (", names(meta_raw)[uasg_first], ")")

need_cov <- c("Sex","Age","StrokeStatus")
missing_cov <- setdiff(need_cov, names(meta_raw))
if (length(missing_cov)) stop("[ERR] Missing columns in meta_selected.csv: ", paste(missing_cov, collapse=", "))

meta <- meta_raw |>
  transmute(
    UASG = trimws(.data[[uasg_first]]),
    Sex  = .data[["Sex"]],
    Age  = .data[["Age"]],
    StrokeStatus = .data[["StrokeStatus"]]
  )
meta$Sex <- factor(meta$Sex)
meta$StrokeStatus <- factor(meta$StrokeStatus)
suppressWarnings(meta$Age <- as.numeric(as.character(meta$Age)))

# ---------- align across expr, geno, cov ----------
uasg_overlap <- Reduce(intersect, list(colnames(expr_mat), unique(geno_df$UASG), meta$UASG))
if (length(uasg_overlap) < 20) {
  warning("Small overlap across expr/genotype/covariates: ", length(uasg_overlap))
}

# Order everything by the same UASG vector
expr_mat_sub <- expr_mat[, uasg_overlap, drop = FALSE]
geno_sub <- geno_df |>
  filter(UASG %in% uasg_overlap) |>
  select(UASG, ADD_G, DOM_G) |>
  distinct() |>
  arrange(match(UASG, uasg_overlap))
cov_sub  <- meta |>
  filter(UASG %in% uasg_overlap) |>
  arrange(match(UASG, uasg_overlap))

stopifnot(identical(colnames(expr_mat_sub), geno_sub$UASG))
stopifnot(identical(geno_sub$UASG, cov_sub$UASG))

# Genotype vectors
ADD <- as.numeric(geno_sub$ADD_G)     # 0/1/2 dosage
DOM <- as.numeric(geno_sub$DOM_G)     # 0/1 (AG/GG vs AA)

# Remove samples with any NA in covariates or genotype (per-model)
keep_add <- !is.na(ADD) & !is.na(cov_sub$Age) & !is.na(cov_sub$Sex) & !is.na(cov_sub$StrokeStatus)
keep_dom <- !is.na(DOM) & !is.na(cov_sub$Age) & !is.na(cov_sub$Sex) & !is.na(cov_sub$StrokeStatus)

X_add <- expr_mat_sub[, keep_add, drop = FALSE];  ADD_k <- ADD[keep_add]
X_dom <- expr_mat_sub[, keep_dom, drop = FALSE];  DOM_k <- DOM[keep_dom]
cov_add <- cov_sub[keep_add, , drop = FALSE]
cov_dom <- cov_sub[keep_dom, , drop = FALSE]

# ---------- design matrices (adjusted) ----------
# Expression ~ 1 + Genotype + Sex + Age + StrokeStatus
design_ADD <- model.matrix(~ ADD_k + cov_add$Sex + cov_add$Age + cov_add$StrokeStatus)
colnames(design_ADD)[2] <- "ADD"   # rename for readability

design_DOM <- model.matrix(~ DOM_k + cov_dom$Sex + cov_dom$Age + cov_dom$StrokeStatus)
colnames(design_DOM)[2] <- "DOM"

# ---------- fit limma + save (with GeneSymbol merged) ----------
fit_and_save <- function(X, design, coef_name, prefix, annot_df) {
  fit <- lmFit(X, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  tt  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  tt$ProbeID <- rownames(tt)
  tt <- tt |>
    left_join(annot_df, by = "ProbeID") |>
    relocate(ProbeID, GeneSymbol, .before = 1)
  out <- file.path(out_dir, paste0("limma_", prefix, "_adjcov_results.csv"))
  write.csv(tt, out, row.names = FALSE)
  message("[ok] Wrote: ", out)
  invisible(tt)
}

add_res <- fit_and_save(X_add, design_ADD, "ADD", "ADD_G", annot_df)
dom_res <- fit_and_save(X_dom, design_DOM, "DOM", "DOM_G", annot_df)

# ---------- run summary ----------
summary_path <- file.path(out_dir, "run_summary_adjcov.txt")
sink(summary_path)
cat("Expression–Genotype association (limma) — adjusted for Age, Sex, StrokeStatus\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Expr samples (total):", ncol(expr_mat), "\n")
cat("Genotype UASGs (total):", nrow(geno_df), "\n")
cat("Overlap (expr ∩ geno ∩ cov):", length(uasg_overlap), "\n")
cat("ADD model samples used:", ncol(X_add), "\n")
cat("DOM model samples used:", ncol(X_dom), "\n")
cat("Genes tested (rows):", nrow(expr_mat_sub), "\n\n")
cat("Model ADD_G: Expression ~ 1 + ADD + Sex + Age + StrokeStatus\n")
cat("Model DOM_G: Expression ~ 1 + DOM + Sex + Age + StrokeStatus\n\n")
cat("ADD_G FDR<0.05 hits:", sum(add_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("DOM_G FDR<0.05 hits:", sum(dom_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
sink()
message("[ok] Wrote: ", summary_path)
message("[done] Adjusted association run complete.")
