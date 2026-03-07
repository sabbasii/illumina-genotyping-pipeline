#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# scripts/eqtl/05_inspect_results/curated_gene_list_analysis/06_test_clinical_outcome_by_snp.R
#
# Purpose
#   Test whether genotype at a SNP is associated with clinical stroke outcomes
#   using metadata/clinical_data.csv, merged to genotypes from output/eqtl/SNP_overlap.txt.
#
# Outcomes
#   1) mRS_90d  (linear):   mRS_90d  ~ genotype (+ optional covariates)
#   2) dmRS_90d (linear):   dmRS_90d ~ genotype (+ optional covariates)
#   3) good_outcome (logistic): I(mRS_90d <= 2) ~ genotype (+ optional covariates)
#
# Models
#   A) Additive genotype: genotype coded 0 / 1 / 2
#   B) Genotype-group:    genotype_f (factor with levels 0/1/2; overall 3-group test)
#
# Filtering
#   --sex male|female
#   --diagnosis <value[,value2,...]>   (case-insensitive; exact match after trimming)
#
# Covariates
#   Default: age + sex + ancestry (included only if present)
#   Disable with: --covariates false
#   NOTE: Any covariate with <2 levels after filtering is automatically dropped
#         (prevents contrasts errors when you subset to one sex, etc.)
#
# SNP selection (choose one)
#   --snp rsXXXX
#   --pde-gene PDE5A   (uses top eQTL SNP for that PDE gene from eqtl_all_pde.tsv)
#
# Output
#   All runs are saved under ONE folder:
#     output/eqtl/results/inspect/pde/clinical_by_snp/<BASE_LABEL>/
#
#   Each run writes a TSV with a suffix describing the filters:
#     <BASE_LABEL>__sex-female__dx-ischemic_stroke__cov-true.tsv
#
# Console
#   Prints genotype x diagnosis counts (no file is saved for this table).
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

snp_in        <- get_arg("--snp", "")
pde_gene      <- toupper(get_arg("--pde-gene", ""))
add_cov       <- tolower(get_arg("--covariates", "true")) %in% c("true", "t", "1", "yes", "y")
sex_filter    <- tolower(trimws(get_arg("--sex", "")))
diag_filter_s <- trimws(get_arg("--diagnosis", ""))

if (!nzchar(snp_in) && !nzchar(pde_gene)) stop("Provide either --snp rsXXXX or --pde-gene PDE5A")
if (nzchar(snp_in) && nzchar(pde_gene)) stop("Provide only one of --snp or --pde-gene (not both).")
if (nzchar(sex_filter) && !(sex_filter %in% c("male", "female"))) stop("--sex must be 'male' or 'female' (or omit it).")

diag_keep <- character()
if (nzchar(diag_filter_s)) {
  diag_keep <- tolower(trimws(unlist(strsplit(diag_filter_s, ",", fixed = TRUE))))
  diag_keep <- diag_keep[nzchar(diag_keep)]
}

clinical_file <- file.path(REPO_ROOT, "metadata/clinical_data.csv")
snp_file      <- file.path(REPO_ROOT, "output/eqtl/SNP_overlap.txt")
eqtl_pde_file <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/eqtl_all_pde.tsv")

if (!file.exists(clinical_file)) stop("Missing file: ", clinical_file)
if (!file.exists(snp_file)) stop("Missing file: ", snp_file)

# ------------------------------------------------------------
# Determine SNP to test
# ------------------------------------------------------------
top_snp <- snp_in
top_p   <- NA_real_

if (!nzchar(top_snp)) {
  if (!file.exists(eqtl_pde_file)) stop("Missing file needed for --pde-gene: ", eqtl_pde_file)
  hits <- fread(eqtl_pde_file)
  if (!all(c("SNP", "gene") %in% names(hits))) stop("eqtl_all_pde.tsv must include columns: SNP, gene")

  p_col <- NULL
  for (cand in c("pvalue", "p-value", "p.value", "PValue", "P", "p")) {
    if (cand %in% names(hits)) { p_col <- cand; break }
  }
  if (is.null(p_col)) stop("No p-value column found in: ", eqtl_pde_file)

  hits[, geneU := toupper(gene)]
  hits[, P := suppressWarnings(as.numeric(get(p_col)))]
  hits <- hits[is.finite(P) & P > 0 & P <= 1]

  sub <- hits[geneU == pde_gene]
  if (nrow(sub) == 0) stop("No hits found for PDE gene: ", pde_gene)

  setorder(sub, P)
  top_snp <- sub$SNP[1]
  top_p   <- sub$P[1]
}

BASE_LABEL <- if (nzchar(pde_gene)) paste0(pde_gene, "_", top_snp) else top_snp

# ------------------------------------------------------------
# Read genotype row
# ------------------------------------------------------------
SNP <- fread(snp_file)
if (!("snpid" %in% names(SNP))) stop("SNP_overlap.txt must have first column named 'snpid'")
setnames(SNP, "snpid", "SNP")

row <- SNP[SNP == top_snp]
if (nrow(row) == 0) stop("SNP not found in SNP_overlap.txt: ", top_snp)

sample_cols <- setdiff(names(SNP), "SNP")
g <- suppressWarnings(as.numeric(unlist(row[, ..sample_cols])))

geno <- data.table(sample = sample_cols, genotype = g)
geno <- geno[is.finite(genotype) & genotype %in% c(0, 1, 2)]
if (nrow(geno) < 20) stop("Too few samples with genotype 0/1/2 after filtering: ", nrow(geno))

# ------------------------------------------------------------
# Read clinical + harmonize IDs
# ------------------------------------------------------------
clin <- fread(clinical_file)
for (nm in c("sample_name","sex","ancestry","age","diagnosis","mRS_90d","dmRS_90d")) {
  if (!(nm %in% names(clin))) stop("clinical_data.csv missing required column: ", nm)
}

clin[, sample_raw  := as.character(sample_name)]
clin[, sample_dash := gsub("_", "-", sample_raw)]
clin[, sample_usg  := gsub("-", "_", sample_raw)]

m_raw  <- sum(clin$sample_raw  %chin% geno$sample)
m_dash <- sum(clin$sample_dash %chin% geno$sample)
m_usg  <- sum(clin$sample_usg  %chin% geno$sample)

if (max(m_raw, m_dash, m_usg) == 0) {
  stop(
    "No overlap between clinical sample_name and SNP sample columns.\n",
    "Example clinical: ", clin$sample_name[1], "\n",
    "Example SNP col:  ", geno$sample[1]
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

# ------------------------------------------------------------
# Merge + filters
# ------------------------------------------------------------
df <- merge(geno, clin, by = "sample", all.x = TRUE)
df[, genotype := suppressWarnings(as.numeric(genotype))]
df <- df[is.finite(genotype) & genotype %in% c(0, 1, 2)]

if (nzchar(sex_filter)) df <- df[sex_l == sex_filter]
if (length(diag_keep))  df <- df[diagnosis_l %in% diag_keep]

if (nrow(df) < 20) stop("Too few samples after filtering (n=", nrow(df), ").")

# Require at least 2 genotype groups
if (uniqueN(df$genotype) < 2) stop("Need >=2 genotype groups after filtering to run models.")

df[, genotype_f := factor(genotype, levels = c(0, 1, 2))]

# Outcomes
df[, mRS_90d  := suppressWarnings(as.numeric(mRS_90d))]
df[, dmRS_90d := suppressWarnings(as.numeric(dmRS_90d))]
df[, age      := suppressWarnings(as.numeric(age))]
df[, good_outcome := as.integer(is.finite(mRS_90d) & mRS_90d <= 2)]

# ------------------------------------------------------------
# Build covariates (drop any with <2 levels after filtering)
# ------------------------------------------------------------
keep_cov <- function(v, dat) {
  if (!(v %in% names(dat))) return(FALSE)
  if (v == "age") {
    x <- dat[[v]]
    x <- x[is.finite(x)]
    return(length(x) >= 3 && stats::sd(x) > 0)
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

rhs_add <- paste(c("genotype", cov_terms), collapse = " + ")
rhs_fac <- paste(c("genotype_f", cov_terms), collapse = " + ")

# ------------------------------------------------------------
# Safe model fitting helpers
# ------------------------------------------------------------
safe_lm <- function(formula, dat) {
  mf <- model.frame(formula, data = dat, na.action = na.omit)
  if (nrow(mf) < 30) return(list(n = nrow(mf), beta = NA_real_, p = NA_real_))
  fit <- lm(formula, data = mf)
  sm <- summary(fit)$coefficients

  if ("genotype" %in% rownames(sm)) {
    return(list(n = nrow(mf), beta = sm["genotype","Estimate"], p = sm["genotype","Pr(>|t|)"]))
  }

  a <- anova(fit)
  rn <- grep("^genotype_f", rownames(a), value = TRUE)
  p <- if (length(rn)) a[rn[1], "Pr(>F)"] else NA_real_
  list(n = nrow(mf), beta = NA_real_, p = p)
}

safe_glm <- function(formula, dat) {
  mf <- model.frame(formula, data = dat, na.action = na.omit)
  if (nrow(mf) < 30) return(list(n = nrow(mf), beta = NA_real_, p = NA_real_))
  fit <- glm(formula, data = mf, family = binomial())
  sm <- summary(fit)$coefficients

  if ("genotype" %in% rownames(sm)) {
    return(list(n = nrow(mf), beta = sm["genotype","Estimate"], p = sm["genotype","Pr(>|z|)"]))
  }

  d <- drop1(fit, test = "Chisq")
  rn <- grep("^genotype_f", rownames(d), value = TRUE)
  p <- if (length(rn)) d[rn[1], "Pr(>Chi)"] else NA_real_
  list(n = nrow(mf), beta = NA_real_, p = p)
}

# ------------------------------------------------------------
# Fit models
# ------------------------------------------------------------
f_mrs_add  <- as.formula(paste("mRS_90d ~", rhs_add))
f_dmrs_add <- as.formula(paste("dmRS_90d ~", rhs_add))
f_good_add <- as.formula(paste("good_outcome ~", rhs_add))

f_mrs_fac  <- as.formula(paste("mRS_90d ~", rhs_fac))
f_dmrs_fac <- as.formula(paste("dmRS_90d ~", rhs_fac))
f_good_fac <- as.formula(paste("good_outcome ~", rhs_fac))

res_mrs_add  <- safe_lm(f_mrs_add, df)
res_dmrs_add <- safe_lm(f_dmrs_add, df)
res_good_add <- safe_glm(f_good_add, df)

res_mrs_fac  <- safe_lm(f_mrs_fac, df)
res_dmrs_fac <- safe_lm(f_dmrs_fac, df)
res_good_fac <- safe_glm(f_good_fac, df)

# ------------------------------------------------------------
# Print genotype x diagnosis distribution (no file saved)
# ------------------------------------------------------------
cat("\nGenotype x diagnosis counts (filtered dataset):\n")
print(df[, .N, by = .(diagnosis = diagnosis_l, genotype)][order(diagnosis, genotype)])

# ------------------------------------------------------------
# Summaries + output
# ------------------------------------------------------------
geno_counts <- df[, .N, by = genotype][order(genotype)]
geno_counts_str <- paste0(geno_counts$genotype, "=", geno_counts$N, collapse = ", ")

n_genotyped <- nrow(df)
n_mrs       <- sum(is.finite(df$mRS_90d))
n_dmrs      <- sum(is.finite(df$dmRS_90d))
n_good      <- sum(!is.na(df$good_outcome))

cov_used <- if (length(cov_terms)) paste(cov_terms, collapse = ",") else "none"

out <- data.table(
  base_label = BASE_LABEL,
  snp = top_snp,
  top_snp_p = top_p,

  sex_filter = if (nzchar(sex_filter)) sex_filter else "all",
  diagnosis_filter = if (length(diag_keep)) paste(diag_keep, collapse = ",") else "all",
  covariates = cov_used,

  n_genotyped = n_genotyped,
  genotype_counts = geno_counts_str,
  n_mRS_90d = n_mrs,
  n_dmRS_90d = n_dmrs,
  n_good_outcome = n_good,

  mRS_add_beta = res_mrs_add$beta,
  mRS_add_p    = res_mrs_add$p,
  mRS_add_n    = res_mrs_add$n,

  dmRS_add_beta = res_dmrs_add$beta,
  dmRS_add_p    = res_dmrs_add$p,
  dmRS_add_n    = res_dmrs_add$n,

  good_add_logodds_beta = res_good_add$beta,
  good_add_p            = res_good_add$p,
  good_add_n            = res_good_add$n,

  mRS_factor_p  = res_mrs_fac$p,
  mRS_factor_n  = res_mrs_fac$n,
  dmRS_factor_p = res_dmrs_fac$p,
  dmRS_factor_n = res_dmrs_fac$n,
  good_factor_p = res_good_fac$p,
  good_factor_n = res_good_fac$n
)

# One folder per SNP/PDE-label only (no per-filter subfolders)
out_dir <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde", "clinical_by_snp", BASE_LABEL)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# File suffix to avoid overwrites
suffix_parts <- c(
  paste0("sex-", if (nzchar(sex_filter)) sex_filter else "all"),
  paste0("dx-", if (length(diag_keep)) paste(diag_keep, collapse = "-") else "all"),
  paste0("cov-", if (add_cov) "true" else "false")
)
suffix <- paste(suffix_parts, collapse = "__")

outfile <- file.path(out_dir, paste0(BASE_LABEL, "__", suffix, ".tsv"))
fwrite(out, outfile, sep = "\t")

cat("\nBASE_LABEL: ", BASE_LABEL, "\n", sep = "")
cat("Saved: ", outfile, "\n", sep = "")
cat("Covariates used: ", cov_used, "\n", sep = "")
cat("Genotypes: ", geno_counts_str, "\n", sep = "")