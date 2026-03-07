#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# scripts/eqtl/05_inspect_results/curated_gene_list_analysis/08_test_clinical_outcome_by_signature.R
#
# Purpose
#   Compute a transcriptional signature score per sample:
#     signature_score = mean(zscore(expression of signature genes))
#   then test whether the score predicts clinical outcomes:
#     mRS_90d  ~ signature_score (+ covariates)
#     dmRS_90d ~ signature_score (+ covariates)
#     good_outcome (mRS_90d<=2) ~ signature_score (+ covariates)
#
#   Optional: also include genotype for a SNP (additive 0/1/2) to test:
#     outcome ~ signature_score + genotype (+ covariates)
#
# Inputs
#   - output/eqtl/GE_overlap.txt
#   - metadata/clinical_data.csv
#   - Signature gene list (default):
#       input_data/target_lists/stroke_inflammation_signature_human.txt
#   - Optional genotype:
#       output/eqtl/SNP_overlap.txt (if --snp is provided)
#
# Filtering
#   --sex male|female
#   --diagnosis <value[,value2,...]>  (case-insensitive exact match)
#
# Covariates
#   Default: age + sex + ancestry (included only if present AND have >=2 levels after filtering)
#   Disable with: --covariates false
#
# Output
#   output/eqtl/results/inspect/pde/clinical_by_signature/<LABEL>/
#     results_dx-..._sex-..._cov-..._snp-....tsv
#     latest.tsv
#
# How to run
#   source scripts/00_config.sh
#
#   # signature only
#   Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/08_test_clinical_outcome_by_signature.R \
#     --diagnosis ischemic_stroke
#
#   # signature + genotype
#   Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/08_test_clinical_outcome_by_signature.R \
#     --diagnosis ischemic_stroke \
#     --snp rs12212126
#
# Optional arguments
#   --gene-set <path>
#   --sex male|female
#   --diagnosis ischemic_stroke,control
#   --covariates false
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) stop(paste0("Missing value after ", flag))
  args[w[1] + 1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (!nzchar(REPO_ROOT)) stop("REPO_ROOT not set. Source scripts/00_config.sh first.")

DEFAULT_GENE_SET <- file.path("input_data", "target_lists", "stroke_inflammation_signature_human.txt")

gene_set_path <- get_arg("--gene-set", DEFAULT_GENE_SET)
if (!grepl("^/", gene_set_path)) gene_set_path <- file.path(REPO_ROOT, gene_set_path)

sex_filter    <- tolower(trimws(get_arg("--sex", "")))
diag_filter_s <- trimws(get_arg("--diagnosis", ""))
add_cov       <- tolower(get_arg("--covariates", "true")) %in% c("true","t","1","yes","y")
snp_in        <- trimws(get_arg("--snp", ""))  # optional

if (nzchar(sex_filter) && !(sex_filter %in% c("male", "female"))) stop("--sex must be 'male' or 'female' (or omit).")

diag_keep <- character()
if (nzchar(diag_filter_s)) {
  diag_keep <- tolower(trimws(unlist(strsplit(diag_filter_s, ",", fixed = TRUE))))
  diag_keep <- diag_keep[nzchar(diag_keep)]
}

GE_file       <- file.path(REPO_ROOT, "output/eqtl/GE_overlap.txt")
clinical_file <- file.path(REPO_ROOT, "metadata/clinical_data.csv")
SNP_file      <- file.path(REPO_ROOT, "output/eqtl/SNP_overlap.txt")

for (f in c(GE_file, clinical_file, gene_set_path)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}
if (nzchar(snp_in) && !file.exists(SNP_file)) stop("Missing file (needed for --snp): ", SNP_file)

# ------------------------------------------------------------
# Load signature genes
# ------------------------------------------------------------
sig_genes <- fread(gene_set_path, header = FALSE)[[1]]
sig_genes <- unique(toupper(trimws(sig_genes)))
sig_genes <- sig_genes[nzchar(sig_genes)]
if (!length(sig_genes)) stop("Signature gene list is empty: ", gene_set_path)

# ------------------------------------------------------------
# Load expression matrix (genes x samples)
# ------------------------------------------------------------
GE <- fread(GE_file)
if (!("geneid" %in% names(GE))) stop("GE_overlap.txt must have first column named 'geneid'")

GE[, geneU := toupper(geneid)]
GEc <- GE[geneU %in% sig_genes]
if (nrow(GEc) < 3) stop("Too few signature genes found in GE_overlap.txt (found ", nrow(GEc), ").")

ge_samples <- setdiff(names(GE), c("geneid", "geneU"))

GEc2 <- GEc[, c("geneid", ge_samples), with = FALSE]
expr_mat <- as.matrix(GEc2[, -"geneid"])
rownames(expr_mat) <- GEc2$geneid
mode(expr_mat) <- "numeric"

# Drop genes with too few finite values or zero variance
ok_gene <- apply(expr_mat, 1, function(x) {
  x <- x[is.finite(x)]
  length(x) >= 10 && stats::sd(x) > 0
})
expr_mat <- expr_mat[ok_gene, , drop = FALSE]
if (nrow(expr_mat) < 3) stop("Too few usable signature genes after QC (need >=3).")

# Z-score each gene across samples (rows=genes)
z_gene <- t(scale(t(expr_mat)))  # gene-wise z
# Score per sample = mean z across genes
sig_score <- colMeans(z_gene, na.rm = TRUE)

sig_dt <- data.table(sample = names(sig_score), signature_score = as.numeric(sig_score))

# ------------------------------------------------------------
# Load clinical data + harmonize sample IDs
# ------------------------------------------------------------
clin <- fread(clinical_file)
need_cols <- c("sample_name","sex","ancestry","age","diagnosis","mRS_90d","dmRS_90d")
for (nm in need_cols) if (!(nm %in% names(clin))) stop("clinical_data.csv missing required column: ", nm)

clin[, sample_raw  := as.character(sample_name)]
clin[, sample_dash := gsub("_", "-", sample_raw)]
clin[, sample_usg  := gsub("-", "_", sample_raw)]

# decide which mapping gives best overlap with expression sample IDs
m_raw  <- sum(clin$sample_raw   %chin% sig_dt$sample)
m_dash <- sum(clin$sample_dash  %chin% sig_dt$sample)
m_usg  <- sum(clin$sample_usg   %chin% sig_dt$sample)

if (max(m_raw, m_dash, m_usg) == 0) {
  stop(
    "No overlap between clinical sample_name and GE sample columns.\n",
    "Example clinical: ", clin$sample_name[1], "\n",
    "Example GE sample: ", sig_dt$sample[1]
  )
}

if (m_dash >= m_raw && m_dash >= m_usg) {
  clin[, sample := sample_dash]
} else if (m_raw >= m_usg) {
  clin[, sample := sample_raw]
} else {
  clin[, sample := sample_usg]
}

clin[, sex_l := tolower(trimws(as.character(sex)))]
clin[, diagnosis_l := tolower(trimws(as.character(diagnosis)))]
clin[, ancestry_l := tolower(trimws(as.character(ancestry)))]

# Merge signature score into clinical
df <- merge(sig_dt, clin, by = "sample", all.x = TRUE)

# Filters
if (nzchar(sex_filter)) df <- df[sex_l == sex_filter]
if (length(diag_keep))  df <- df[diagnosis_l %in% diag_keep]

# Drop samples with no clinical match (optional but makes counts honest)
df <- df[!is.na(sample_name)]

if (nrow(df) < 30) stop("Too few samples after filtering (n=", nrow(df), ").")

# Outcomes + covariates
df[, mRS_90d  := suppressWarnings(as.numeric(mRS_90d))]
df[, dmRS_90d := suppressWarnings(as.numeric(dmRS_90d))]
df[, age      := suppressWarnings(as.numeric(age))]
df[, good_outcome := as.integer(is.finite(mRS_90d) & mRS_90d <= 2)]

# Optional genotype merge
has_geno <- FALSE
if (nzchar(snp_in)) {
  SNP <- fread(SNP_file)
  if (!("snpid" %in% names(SNP))) stop("SNP_overlap.txt must have first column named 'snpid'")
  setnames(SNP, "snpid", "SNP")

  row <- SNP[SNP == snp_in]
  if (nrow(row) == 0) stop("SNP not found in SNP_overlap.txt: ", snp_in)

  sample_cols <- setdiff(names(SNP), "SNP")
  g <- suppressWarnings(as.numeric(unlist(row[, ..sample_cols])))

  geno <- data.table(sample = sample_cols, genotype = g)
  geno <- geno[is.finite(genotype) & genotype %in% c(0,1,2)]

  # Merge to df (keep signature-scored samples)
  df <- merge(df, geno, by = "sample", all.x = TRUE)
  df[, genotype := suppressWarnings(as.numeric(genotype))]
  has_geno <- TRUE
}

# ------------------------------------------------------------
# Covariates: keep only if informative after filtering
# ------------------------------------------------------------
keep_cov <- function(v, dat) {
  if (!(v %in% names(dat))) return(FALSE)
  if (v == "age") {
    x <- dat[[v]]
    x <- x[is.finite(x)]
    return(length(x) >= 10 && stats::sd(x) > 0)
  }
  x <- tolower(trimws(as.character(dat[[v]])))
  x <- x[nzchar(x)]
  return(length(unique(x)) >= 2)
}

cov_terms <- character()
if (add_cov) {
  if (keep_cov("age", df))      cov_terms <- c(cov_terms, "age")
  if (keep_cov("sex", df))      cov_terms <- c(cov_terms, "sex")
  if (keep_cov("ancestry", df)) cov_terms <- c(cov_terms, "ancestry")
}

rhs_sig <- paste(c("signature_score", cov_terms), collapse = " + ")
rhs_sig_geno <- paste(c("signature_score", "genotype", cov_terms), collapse = " + ")

# ------------------------------------------------------------
# Safe model fits
# ------------------------------------------------------------
safe_lm <- function(formula, dat) {
  mf <- model.frame(formula, data = dat, na.action = na.omit)
  if (nrow(mf) < 30) return(list(n = nrow(mf), beta = NA_real_, p = NA_real_))
  fit <- lm(formula, data = mf)
  sm <- summary(fit)$coefficients
  if (!("signature_score" %in% rownames(sm))) return(list(n=nrow(mf), beta=NA_real_, p=NA_real_))
  list(n = nrow(mf), beta = sm["signature_score","Estimate"], p = sm["signature_score","Pr(>|t|)"])
}

safe_glm <- function(formula, dat) {
  mf <- model.frame(formula, data = dat, na.action = na.omit)
  if (nrow(mf) < 30) return(list(n = nrow(mf), beta = NA_real_, p = NA_real_))
  fit <- glm(formula, data = mf, family = binomial())
  sm <- summary(fit)$coefficients
  if (!("signature_score" %in% rownames(sm))) return(list(n=nrow(mf), beta=NA_real_, p=NA_real_))
  list(n = nrow(mf), beta = sm["signature_score","Estimate"], p = sm["signature_score","Pr(>|z|)"])
}

# Build formulas
f_mrs_sig  <- as.formula(paste("mRS_90d ~", rhs_sig))
f_dmrs_sig <- as.formula(paste("dmRS_90d ~", rhs_sig))
f_good_sig <- as.formula(paste("good_outcome ~", rhs_sig))

# Fit signature-only
r_mrs_sig  <- safe_lm(f_mrs_sig, df)
r_dmrs_sig <- safe_lm(f_dmrs_sig, df)
r_good_sig <- safe_glm(f_good_sig, df)

# Fit signature+genotype if requested and enough genotype data
r_mrs_sg <- r_dmrs_sg <- r_good_sg <- list(n=NA_integer_, beta=NA_real_, p=NA_real_)
sg_ok <- FALSE
if (has_geno) {
  if (sum(is.finite(df$genotype) & df$genotype %in% c(0,1,2)) >= 30) {
    sg_ok <- TRUE
    f_mrs_sg  <- as.formula(paste("mRS_90d ~", rhs_sig_geno))
    f_dmrs_sg <- as.formula(paste("dmRS_90d ~", rhs_sig_geno))
    f_good_sg <- as.formula(paste("good_outcome ~", rhs_sig_geno))
    r_mrs_sg  <- safe_lm(f_mrs_sg, df)
    r_dmrs_sg <- safe_lm(f_dmrs_sg, df)
    r_good_sg <- safe_glm(f_good_sg, df)
  }
}

# ------------------------------------------------------------
# Console printout
# ------------------------------------------------------------
n_all  <- nrow(df)
n_mrs  <- sum(is.finite(df$mRS_90d))
n_dmrs <- sum(is.finite(df$dmRS_90d))
n_good <- sum(!is.na(df$good_outcome))
cov_used <- if (length(cov_terms)) paste(cov_terms, collapse = ",") else "none"

fmt_int <- function(x) ifelse(is.na(x), "NA", formatC(x, format="d", big.mark=","))
fmt_num <- function(x, d=3) ifelse(is.na(x), "NA", formatC(x, format="f", digits=d))
fmt_p   <- function(p) ifelse(is.na(p), "NA", formatC(p, format="e", digits=2))

cat("\n=== Clinical outcome ~ signature score ===\n")
cat(sprintf("Gene set: %s\n", basename(gene_set_path)))
cat(sprintf("Genes: found=%s | used_after_QC=%s\n", fmt_int(nrow(GEc)), fmt_int(nrow(expr_mat))))
cat(sprintf("Filters: sex=%s | diagnosis=%s | covariates=%s | snp=%s\n",
            if (nzchar(sex_filter)) sex_filter else "all",
            if (length(diag_keep)) paste(diag_keep, collapse=",") else "all",
            cov_used,
            if (has_geno) snp_in else "none"))
cat(sprintf("Samples: total=%s | mRS=%s | dmRS=%s | good_outcome=%s\n",
            fmt_int(n_all), fmt_int(n_mrs), fmt_int(n_dmrs), fmt_int(n_good)))

ss <- df[, .(
  n = .N,
  mean = mean(signature_score, na.rm=TRUE),
  sd   = sd(signature_score, na.rm=TRUE),
  min  = min(signature_score, na.rm=TRUE),
  max  = max(signature_score, na.rm=TRUE)
)]
cat(sprintf("Signature score: mean=%s  sd=%s  range=[%s, %s]\n",
            fmt_num(ss$mean,3), fmt_num(ss$sd,3), fmt_num(ss$min,3), fmt_num(ss$max,3)))

cat("\nModel results (beta for signature_score):\n")
cat(sprintf("  mRS_90d   : beta=%s  p=%s  n=%s\n", fmt_num(r_mrs_sig$beta,3),  fmt_p(r_mrs_sig$p),  fmt_int(r_mrs_sig$n)))
cat(sprintf("  dmRS_90d  : beta=%s  p=%s  n=%s\n", fmt_num(r_dmrs_sig$beta,3), fmt_p(r_dmrs_sig$p), fmt_int(r_dmrs_sig$n)))
cat(sprintf("  good<=2   : logodds=%s  p=%s  n=%s\n", fmt_num(r_good_sig$beta,3), fmt_p(r_good_sig$p), fmt_int(r_good_sig$n)))

if (has_geno) {
  if (sg_ok) {
    cat("\nSignature + genotype model (beta for signature_score):\n")
    cat(sprintf("  mRS_90d   : beta=%s  p=%s  n=%s\n", fmt_num(r_mrs_sg$beta,3),  fmt_p(r_mrs_sg$p),  fmt_int(r_mrs_sg$n)))
    cat(sprintf("  dmRS_90d  : beta=%s  p=%s  n=%s\n", fmt_num(r_dmrs_sg$beta,3), fmt_p(r_dmrs_sg$p), fmt_int(r_dmrs_sg$n)))
    cat(sprintf("  good<=2   : logodds=%s  p=%s  n=%s\n", fmt_num(r_good_sg$beta,3), fmt_p(r_good_sg$p), fmt_int(r_good_sg$n)))
  } else {
    cat("\nSignature + genotype: skipped (not enough genotype data after filtering).\n")
  }
}

# ------------------------------------------------------------
# Output table
# ------------------------------------------------------------
LABEL <- "stroke_signature"
if (basename(gene_set_path) != "stroke_inflammation_signature_human.txt") {
  lbl <- basename(gene_set_path)
  lbl <- gsub("\\.[^.]+$", "", lbl)
  lbl <- gsub("[^A-Za-z0-9_\\-]+", "_", lbl)
  LABEL <- paste0("signature_", lbl)
}

out <- data.table(
  label = LABEL,
  gene_set = gene_set_path,
  snp = if (has_geno) snp_in else "",
  sex_filter = if (nzchar(sex_filter)) sex_filter else "all",
  diagnosis_filter = if (length(diag_keep)) paste(diag_keep, collapse = ",") else "all",
  covariates = cov_used,

  n_samples = n_all,
  n_mRS_90d = n_mrs,
  n_dmRS_90d = n_dmrs,
  n_good_outcome = n_good,

  mRS_sig_beta = r_mrs_sig$beta,
  mRS_sig_p    = r_mrs_sig$p,
  mRS_sig_n    = r_mrs_sig$n,

  dmRS_sig_beta = r_dmrs_sig$beta,
  dmRS_sig_p    = r_dmrs_sig$p,
  dmRS_sig_n    = r_dmrs_sig$n,

  good_sig_logodds_beta = r_good_sig$beta,
  good_sig_p            = r_good_sig$p,
  good_sig_n            = r_good_sig$n,

  mRS_sig_geno_beta = if (sg_ok) r_mrs_sg$beta else NA_real_,
  mRS_sig_geno_p    = if (sg_ok) r_mrs_sg$p else NA_real_,
  mRS_sig_geno_n    = if (sg_ok) r_mrs_sg$n else NA_integer_,

  dmRS_sig_geno_beta = if (sg_ok) r_dmrs_sg$beta else NA_real_,
  dmRS_sig_geno_p    = if (sg_ok) r_dmrs_sg$p else NA_real_,
  dmRS_sig_geno_n    = if (sg_ok) r_dmrs_sg$n else NA_integer_,

  good_sig_geno_logodds_beta = if (sg_ok) r_good_sg$beta else NA_real_,
  good_sig_geno_p            = if (sg_ok) r_good_sg$p else NA_real_,
  good_sig_geno_n            = if (sg_ok) r_good_sg$n else NA_integer_
)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
out_dir <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/clinical_by_signature")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Simple filename logic
outfile <- if (has_geno) {
  file.path(out_dir, "stroke_signature_expr_geno_mRS.tsv")
} else {
  file.path(out_dir, "stroke_signature_expr_mRS.tsv")
}

# Write results
fwrite(out, outfile, sep = "\t")

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------
cat("\n=== Results (signature_score effect) ===\n")
cat(sprintf("mRS_90d:  beta=%0.3f  p=%s  n=%d\n",
            out$mRS_sig_beta, format(out$mRS_sig_p, scientific=TRUE, digits=2), out$mRS_sig_n))
cat(sprintf("dmRS_90d: beta=%0.3f  p=%s  n=%d\n",
            out$dmRS_sig_beta, format(out$dmRS_sig_p, scientific=TRUE, digits=2), out$dmRS_sig_n))
cat(sprintf("good<=2:  logodds=%0.3f  p=%s  n=%d\n",
            out$good_sig_logodds_beta, format(out$good_sig_p, scientific=TRUE, digits=2), out$good_sig_n))

if (has_geno) {
  cat("\n(signature + genotype model)\n")
  cat(sprintf("mRS_90d:  beta=%s  p=%s  n=%s\n",
              format(out$mRS_sig_geno_beta, digits=3),
              format(out$mRS_sig_geno_p, scientific=TRUE, digits=2),
              format(out$mRS_sig_geno_n)))
  cat(sprintf("dmRS_90d: beta=%s  p=%s  n=%s\n",
              format(out$dmRS_sig_geno_beta, digits=3),
              format(out$dmRS_sig_geno_p, scientific=TRUE, digits=2),
              format(out$dmRS_sig_geno_n)))
  cat(sprintf("good<=2:  logodds=%s  p=%s  n=%s\n",
              format(out$good_sig_geno_logodds_beta, digits=3),
              format(out$good_sig_geno_p, scientific=TRUE, digits=2),
              format(out$good_sig_geno_n)))
}

cat("\nSaved file:\n")
cat("  ", outfile, "\n\n", sep = "")