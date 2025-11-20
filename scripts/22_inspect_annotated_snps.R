#!/usr/bin/env Rscript
# inspect_genotype_csv.R
#
# Part 3 of the gene annotation workflow:
#   1) glist_hg19_gene_annotation_prepare.sh
#   2) vcf_gene_annotation_export.sh
#   3) inspect_annotated_snps.R  <-- this script (QC / exploration)
#
# This script inspects the SNP–GENE–GENOTYPE CSV produced by:
#   scripts/vcf_gene_annotation_export.sh
#
# Default input (in output/genotype_run1/vcf/):
#   cohort_snps_genes_genotypes.csv
#
# Columns are expected to be:
#   ID, CHROM, POS, REF, ALT, GENE, <sample1>, <sample2>, ...
#
# Usage:
#   cd output/genotype_run1/vcf
#   Rscript inspect_genotype_csv.R
#
#   # or specify a different CSV path:
#   Rscript inspect_genotype_csv.R path/to/your_file.csv
#
# You can also source this in an R console / VS Code for interactive work.

suppressPackageStartupMessages({
  library(data.table)   # fast fread / data.table
})

# ==============================
# 1. Paths and basic setup
# ==============================

args <- commandArgs(trailingOnly = TRUE)

# Default file name (matches vcf_gene_annotation_export.sh output)
default_csv <- "cohort_snps_genes_genotypes.csv"

csv_path <- if (length(args) >= 1) args[[1]] else default_csv

cat("[INFO] Input CSV:", csv_path, "\n")

if (!file.exists(csv_path)) {
  stop("[ERROR] File not found: ", csv_path,
       "\n       Make sure you ran scripts/vcf_gene_annotation_export.sh first,",
       "\n       or pass the CSV path explicitly to this script.")
}

# fread() is fast and memory-efficient; it auto-detects types.
geno_dt <- fread(
  csv_path,
  sep = ",",
  header = TRUE,
  data.table = TRUE
)

# ==============================
# 2. Basic structure checks
# ==============================

cat("\n[INFO] Basic structure\n")
cat("  Rows (SNPs):                 ", nrow(geno_dt), "\n")
cat("  Columns (meta + samples):    ", ncol(geno_dt), "\n\n")

cat("[INFO] Column names (first 20):\n")
print(names(geno_dt)[1:min(20, ncol(geno_dt))])
cat("\n")

cat("[INFO] Glimpse of meta columns (first 10 rows):\n")

# We expect these meta columns from the export script:
meta_cols <- c("ID", "CHROM", "POS", "REF", "ALT", "GENE")
missing_meta <- setdiff(meta_cols, names(geno_dt))

if (length(missing_meta) > 0) {
  cat("[WARN] The following expected meta columns are missing:\n")
  print(missing_meta)
  cat("       Check that the CSV was produced by vcf_gene_annotation_export.sh.\n\n")
} else {
  print(head(geno_dt[, ..meta_cols], 10))
  cat("\n")
}

cat("[INFO] Structure (first few columns):\n")
str(geno_dt[, 1:min(15, ncol(geno_dt)), with = FALSE])
cat("\n")

# ==============================
# 3. Separate meta vs genotype columns
# ==============================

all_cols  <- names(geno_dt)
geno_cols <- setdiff(all_cols, meta_cols)

cat("[INFO] Meta columns:\n")
print(meta_cols)
cat("[INFO] Number of genotype columns (samples):", length(geno_cols), "\n\n")

meta_dt <- geno_dt[, ..meta_cols]
gt_dt   <- if (length(geno_cols) > 0) geno_dt[, ..geno_cols] else NULL

# ==============================
# 4. Inspect GENE annotation
# ==============================

if (!"GENE" %in% names(meta_dt)) {
  cat("[WARN] No GENE column detected. Skipping gene-level summary.\n\n")
} else {
  cat("[INFO] GENE column summary:\n")
  cat("  Example values:\n")
  print(head(meta_dt$GENE, 20))
  cat("\n")

  n_total     <- nrow(meta_dt)
  n_with_gene <- sum(meta_dt$GENE != ".", na.rm = TRUE)
  n_without   <- sum(meta_dt$GENE == ".", na.rm = TRUE)

  cat("[INFO] Gene annotation counts:\n")
  cat("  Total SNPs          :", n_total,      "\n")
  cat("  SNPs with a GENE    :", n_with_gene,  "\n")
  cat("  SNPs without GENE . :", n_without,    "\n")
  if (n_total > 0) {
    cat("  Proportion with GENE:", sprintf("%.3f", n_with_gene / n_total), "\n")
  }
  cat("\n")

  uniq_genes <- sort(unique(meta_dt$GENE[meta_dt$GENE != "."]))
  cat("[INFO] Number of unique genes with at least one SNP:", length(uniq_genes), "\n")
  cat("       First 20 genes:\n")
  print(head(uniq_genes, 20))
  cat("\n")
}

# ==============================
# 5. Inspect genotypes for a few samples
# ==============================

if (length(geno_cols) > 0) {
  cat("[INFO] Quick genotype check for the first 3 samples:\n")

  sample_cols <- geno_cols[1:min(3, length(geno_cols))]
  # Show ID, CHROM, POS + a few sample columns
  print(head(geno_dt[, c("ID", "CHROM", "POS", sample_cols), with = FALSE], 10))
  cat("\n")

  for (s in sample_cols) {
    cat("[INFO] Genotype distribution for sample:", s, "\n")
    print(table(geno_dt[[s]], useNA = "ifany"))
    cat("\n")
  }
} else {
  cat("[WARN] No genotype columns detected. Check that your CSV includes sample columns.\n\n")
}

# ==============================
# 6. Example: extract all SNPs for one gene
# ==============================

gene_of_interest <- "FAM87B"  # change as needed

if ("GENE" %in% names(meta_dt) && gene_of_interest %in% meta_dt$GENE) {
  cat("[INFO] Example: SNPs annotated to gene:", gene_of_interest, "\n")

  snps_gene <- geno_dt[GENE == gene_of_interest]
  cat("  Number of SNPs for this gene:", nrow(snps_gene), "\n")
  cat("  Preview (meta cols + first 3 samples):\n")

  first_samples <- geno_cols[1:min(3, length(geno_cols))]
  cols_to_show  <- c(meta_cols, first_samples)
  cols_to_show  <- intersect(cols_to_show, names(snps_gene))

  print(snps_gene[, ..cols_to_show])
  cat("\n")
} else {
  cat("[INFO] Gene", gene_of_interest,
      "not found in GENE column (or GENE column missing).\n\n")
}

cat("[DONE] CSV inspection complete.\n")
