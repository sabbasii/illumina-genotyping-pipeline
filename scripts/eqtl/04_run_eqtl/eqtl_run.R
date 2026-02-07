#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(MatrixEQTL)
})

# ------------------------------------------------------------
# eqtl_run.R
#
# Runs MatrixEQTL in two modes:
#   1) ALL eQTL (single output file)
#   2) CIS + TRANS eQTL (separate output files)
#
# Update: supports --base-dir so the same script can run on either:
#   - output/genotype_run1/eqtl
#   - output/genotype_run1/eqtl_allSNPs
#
# Run examples:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/04_run_eqtl/eqtl_run.R
#
#   # Run on all-SNP inputs
#   Rscript scripts/eqtl/04_run_eqtl/eqtl_run.R --base-dir output/genotype_run1/eqtl_allSNPs
# ------------------------------------------------------------

# ---------------------------
# CLI args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste0("Missing value after ", flag))
  args[hit + 1]
}

has_flag <- function(flag) flag %in% args

base_rel <- get_arg("--base-dir", "output/genotype_run1/eqtl")
NO_COV <- has_flag("--no-covariates")

# ---------------------------
# config (REPO_ROOT)
# ---------------------------
REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

base.dir <- file.path(REPO_ROOT, base_rel)

# Inputs (prepared earlier in pipeline)
SNP_file_name        <- file.path(base.dir, "SNP_overlap.txt")
expression_file_name <- file.path(base.dir, "GE_overlap.txt")

covariates_file_name <- file.path(base.dir, "Covariates_numeric.txt")

snpsloc_file_name <- file.path(base.dir, "snpsloc.txt")
geneloc_file_name <- file.path(base.dir, "geneloc.txt")

out.dir <- file.path(base.dir, "results")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Choose model:
#   modelLINEAR (additive; recommended default)
#   modelANOVA
#   modelLINEAR_CROSS
useModel <- modelLINEAR

# Output thresholds (smaller = fewer rows written)
pvOutputThreshold       <- 1e-3   # ALL
pvOutputThreshold.cis   <- 1e-5   # CIS
pvOutputThreshold.trans <- 1e-5   # TRANS

# CIS distance (bp): common choices 1e6 (1Mb) or 5e5 (500kb)
cisDist <- 1e6

errorCovariance <- numeric()

# ---------------------------
# helpers
# ---------------------------
load_sliced <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  sd <- SlicedData$new()
  sd$fileDelimiter <- "\t"
  sd$fileOmitCharacters <- "NA"
  sd$fileSkipRows <- 1
  sd$fileSkipColumns <- 1
  sd$fileSliceSize <- 2000
  sd$LoadFile(path)
  sd
}

load_covariates <- function(path) {
  if (NO_COV) return(SlicedData$new())
  if (!file.exists(path)) {
    message("NOTE: Covariates_numeric.txt not found, running with NO covariates: ", path)
    return(SlicedData$new())
  }
  load_sliced(path)
}

# ---------------------------
# sanity checks
# ---------------------------
need_files <- c(SNP_file_name, expression_file_name, snpsloc_file_name, geneloc_file_name)
missing <- need_files[!file.exists(need_files)]
if (length(missing) > 0) {
  stop("Missing required input(s):\n  ", paste(missing, collapse = "\n  "))
}

# ---------------------------
# load data
# ---------------------------
snps <- load_sliced(SNP_file_name)
gene <- load_sliced(expression_file_name)
cvrt <- load_covariates(covariates_file_name)

# Locations (required for cis/trans)
snpspos <- read.table(snpsloc_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genepos <- read.table(geneloc_file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# MatrixEQTL expects:
# snpspos: 3 cols -> snpid, chr, pos
# genepos: 4 cols -> geneid, chr, left, right
stopifnot(all(c("snpid", "chr", "pos") %in% names(snpspos)))
stopifnot(all(c("geneid", "chr", "left", "right") %in% names(genepos)))

# ---------------------------
# 1) ALL eQTL (single output)
# ---------------------------
out_all <- file.path(out.dir, "eqtl_all.tsv")

me_all <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = out_all,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

# ---------------------------
# 2) CIS + TRANS eQTL
# ---------------------------
out_cis   <- file.path(out.dir, "eqtl_cis.tsv")
out_trans <- file.path(out.dir, "eqtl_trans.tsv")

me_ct <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = out_trans,                 # TRANS
  pvOutputThreshold = pvOutputThreshold.trans,
  output_file_name.cis = out_cis,               # CIS
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

# ---------------------------
# run summary
# ---------------------------
summary_path <- file.path(out.dir, "eqtl_run_summary.txt")

model_name <- if (identical(useModel, modelLINEAR)) "modelLINEAR" else
  if (identical(useModel, modelANOVA)) "modelANOVA" else
    if (identical(useModel, modelLINEAR_CROSS)) "modelLINEAR_CROSS" else "unknown"

cat(
  "Matrix eQTL run summary\n",
  "----------------------\n",
  "Base dir: ", base.dir, "\n",
  "Model: ", model_name, "\n",
  "cisDist (bp): ", cisDist, "\n",
  "Covariates used: ", if (NO_COV) "NO (--no-covariates)" else if (file.exists(covariates_file_name)) "YES" else "NO (file missing)", "\n",
  "\nInputs:\n",
  "  SNP_overlap: ", SNP_file_name, "\n",
  "  GE_overlap:  ", expression_file_name, "\n",
  "  snpsloc:     ", snpsloc_file_name, "\n",
  "  geneloc:     ", geneloc_file_name, "\n",
  if (!NO_COV) paste0("  covariates:  ", covariates_file_name, "\n") else "",
  "\nOutputs:\n",
  "  ALL:   ", out_all, "\n",
  "  CIS:   ", out_cis, "\n",
  "  TRANS: ", out_trans, "\n",
  "\nThresholds:\n",
  "  ALL pvOutputThreshold:       ", pvOutputThreshold, "\n",
  "  CIS pvOutputThreshold.cis:   ", pvOutputThreshold.cis, "\n",
  "  TRANS pvOutputThreshold:     ", pvOutputThreshold.trans, "\n",
  file = summary_path, sep = ""
)

cat(
  "Done.\nWrote:\n  ",
  out_all, "\n  ",
  out_cis, "\n  ",
  out_trans, "\n  ",
  summary_path, "\n",
  sep = ""
)

# ------------------------------------------------------------
# OUTPUTS (Default paths; can be changed with --base-dir)
# ------------------------------------------------------------
#   $REPO_ROOT/output/genotype_run1/eqtl/results/
#     - eqtl_all.tsv
#     - eqtl_cis.tsv
#     - eqtl_trans.tsv
#     - eqtl_run_summary.txt
# ------------------------------------------------------------
