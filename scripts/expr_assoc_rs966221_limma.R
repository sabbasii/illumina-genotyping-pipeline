#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
  library(readr)
  library(dplyr)
})

# ---- paths from environment (00_config.sh) ----
EXPR_OUT_DIR <- Sys.getenv("EXPR_OUT_DIR")
PLINK_DIR    <- Sys.getenv("PLINK_DIR")
PHENO_DIR    <- Sys.getenv("PHENO_DIR")

stopifnot(nzchar(EXPR_OUT_DIR), nzchar(PLINK_DIR), nzchar(PHENO_DIR))

expr_csv   <- file.path(EXPR_OUT_DIR, "expr_selected.csv")
meta_csv   <- file.path(EXPR_OUT_DIR, "meta_selected.csv")   # not used in 1st pass
geno_tsv   <- file.path(PLINK_DIR,    "tmp", "rs966221_genotypes.tsv")
keep_list  <- file.path(PHENO_DIR,    "iid_selected.keep")

out_dir    <- file.path(EXPR_OUT_DIR, "assoc_rs966221")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load data ----
expr_df <- read_csv(expr_csv, show_col_types = FALSE)
geno_df <- read_tsv(geno_tsv, show_col_types = FALSE, col_types = "ccii") # IID,GT,ADD_G,DOM_G
keep_ids <- read_tsv(keep_list, col_names = c("FID","IID"), show_col_types = FALSE)

# ---- prepare expression matrix ----
# expr_selected.csv: rows = probes/genes, first column = probe IDs, remaining columns = UASG-#### samples
stopifnot(ncol(expr_df) >= 2)
colnames(expr_df)[1] <- "ProbeID"

# Samples as in header:
sample_cols <- setdiff(colnames(expr_df), "ProbeID")

# Ensure matrix
expr_mat <- as.matrix(expr_df[, sample_cols])
rownames(expr_mat) <- expr_df$ProbeID
mode(expr_mat) <- "numeric"

# Optional log2 guard (if values look like raw intensities)
rng <- range(expr_mat, finite = TRUE)
if (rng[2] > 100) {
  message("[info] Expression range suggests non-log values; applying log2(x+1).")
  expr_mat <- log2(expr_mat + 1)
}

# ---- link sample IDs: UASG (expression) -> IID (genotypes) ----
# Your IID format is SentrixBarcode_A_SentrixPosition_A in VCF/PLINK world,
# but the genotype table you generated (rs966221_genotypes.tsv) should already use IID consistent with PLINK.
# For expression, headers are UASG-####. We need a UASG <-> IID map.
# We derive UASG from keep_list by merging with SampleSheet later; for now, use UASG present in expression
# and rely on rs966221_genotypes.tsv having a UASG column if you added it. If not, add a quick mapping step:

# Fast path: if geno_df already includes UASG column, use it; else, try to parse from IID via a map file if you have one.
if (!"UASG" %in% colnames(geno_df)) {
  # Try to read a previously made IID<->UASG map if present:
  map_path <- file.path(PHENO_DIR, "iid_to_uasg.tsv")
  if (file.exists(map_path)) {
    map_df <- read_tsv(map_path, show_col_types = FALSE, col_types = "ccc")
    geno_df <- geno_df %>% left_join(map_df, by = "IID")
  }
}

if (!"UASG" %in% colnames(geno_df)) {
  stop("Genotype table lacks UASG and no iid_to_uasg.tsv found. Please provide PHENO_DIR/iid_to_uasg.tsv with columns IID,UASG.")
}

# Intersect on UASG
expr_uasg <- sample_cols
geno_uasg <- unique(geno_df$UASG)
uasg_overlap <- intersect(expr_uasg, geno_uasg)

if (length(uasg_overlap) < 50) {
  warning("Small overlap between expression and genotype sample IDs: ", length(uasg_overlap))
}

# subset/align
expr_mat_sub <- expr_mat[, uasg_overlap, drop = FALSE]
geno_sub <- geno_df %>% filter(UASG %in% uasg_overlap) %>%
  select(UASG, ADD_G, DOM_G) %>%
  distinct()

# Order alignment
geno_sub <- geno_sub %>% arrange(UASG)
expr_mat_sub <- expr_mat_sub[, geno_sub$UASG, drop = FALSE]
stopifnot(identical(colnames(expr_mat_sub), geno_sub$UASG))

# ---- design matrices ----
# Primary: ADD_G (numeric 0/1/2)
ADD <- as.numeric(geno_sub$ADD_G)
if (any(is.na(ADD))) {
  keep <- !is.na(ADD)
  expr_mat_sub <- expr_mat_sub[, keep, drop = FALSE]
  ADD <- ADD[keep]
}

design_ADD <- model.matrix(~ ADD)  # Intercept + ADD

# Sensitivity: DOM_G (0/1)
DOM <- as.numeric(geno_sub$DOM_G)
if (any(is.na(DOM))) {
  keep2 <- !is.na(DOM)
  # Keep alignment with ADD analysis, but run separately
  expr_mat_dom <- expr_mat_sub[, keep2, drop = FALSE]
  DOM <- DOM[keep2]
} else {
  expr_mat_dom <- expr_mat_sub
}
design_DOM <- model.matrix(~ DOM)

# ---- fit limma (microarray-like matrix; voom not required) ----
fit_and_save <- function(X, design, prefix) {
  fit <- lmFit(X, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")  # coef 2 = ADD or DOM
  tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
  out_path <- file.path(out_dir, paste0("limma_", prefix, "_results.csv"))
  write.csv(tt, out_path, row.names = TRUE)
  message("[ok] Wrote: ", out_path)
}

fit_and_save(expr_mat_sub, design_ADD, "ADD_G")
fit_and_save(expr_mat_dom, design_DOM, "DOM_G")

# ---- QC summaries ----
sink(file.path(out_dir, "run_summary.txt"))
cat("Expression–Genotype association (limma)\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Samples overlapped (UASG):", length(colnames(expr_mat_sub)), "\n")
cat("Genes tested:", nrow(expr_mat_sub), "\n")
cat("Design ADD_G: Intercept + ADD (0/1/2; risk allele = G)\n")
cat("Design DOM_G: Intercept + DOM (1 if AG/GG else 0)\n\n")

# Simple hit summaries
add_res <- read.csv(file.path(out_dir, "limma_ADD_G_results.csv"))
dom_res <- read.csv(file.path(out_dir, "limma_DOM_G_results.csv"))
cat("ADD_G FDR<0.05 hits:", sum(add_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("DOM_G FDR<0.05 hits:", sum(dom_res$adj.P.Val < 0.05, na.rm = TRUE), "\n")
sink()

message("[done] Association run complete.")
