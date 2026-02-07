#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
})

# ------------------------------------------------------------
# make_geneloc_from_gtf.R
#
# Build Matrix-eQTL geneloc.txt from an Ensembl GTF, filtered to
# the gene IDs actually present in GE.txt (gene symbols).
#
# Inputs:
#   - output/genotype_run1/eqtl/GE.txt            (geneid column)
#   - reference/annotation/Homo_sapiens.GRCh37.87.gtf.gz
#
# Output:
#   - output/genotype_run1/eqtl/geneloc.txt
#     columns: geneid, chr, left, right
#     chr format: chr1..chr22, chrX, chrY, chrMT
#
# Run:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/01_build_inputs/make_geneloc_from_gtf.R
# ------------------------------------------------------------

# ---- paths (REPO_ROOT) ----
REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

ge_path  <- file.path(REPO_ROOT, "output/genotype_run1/eqtl/GE.txt")
gtf_path <- file.path(REPO_ROOT, "reference/annotation/Homo_sapiens.GRCh37.87.gtf.gz")

out_dir  <- file.path(REPO_ROOT, "output/genotype_run1/eqtl")
out_path <- file.path(out_dir, "geneloc.txt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ge_path)) stop("Missing GE.txt: ", ge_path)
if (!file.exists(gtf_path)) stop("Missing GTF: ", gtf_path)

# ---- load gene IDs from GE.txt (authoritative) ----
ge <- fread(ge_path, sep = "\t", header = TRUE, data.table = FALSE, quote = "")
if (!("geneid" %in% names(ge))) stop("GE.txt must have a 'geneid' column.")

symbols <- unique(ge$geneid)
symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

# ---- import GTF (genes only) ----
gtf <- import(gtf_path)                 # rtracklayer reads .gtf.gz
gtf_genes <- gtf[gtf$type == "gene"]

gene_name <- mcols(gtf_genes)$gene_name
chr_raw   <- as.character(seqnames(gtf_genes))
left      <- start(gtf_genes)
right     <- end(gtf_genes)

# Keep standard chromosomes; add "chr" prefix to match snpsloc (chr1, chrX, chrMT)
keep_chr <- chr_raw %in% c(as.character(1:22), "X", "Y", "MT")

df <- data.frame(
  geneid = gene_name,
  chr    = paste0("chr", chr_raw),
  left   = left,
  right  = right,
  stringsAsFactors = FALSE
)

df <- df[keep_chr & !is.na(df$geneid) & nzchar(df$geneid), , drop = FALSE]

# ---- filter to genes present in GE + de-duplicate (keep widest span per geneid) ----
df <- df[df$geneid %in% symbols, , drop = FALSE]
df$span <- df$right - df$left

# If a gene appears multiple times, keep the entry with the largest span
df <- df[order(df$geneid, -df$span), ]
df <- df[!duplicated(df$geneid), c("geneid", "chr", "left", "right")]

# ---- write output ----
fwrite(df, out_path, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:", out_path, "\n")
cat("Mapped genes:", nrow(df), "of", length(symbols), "unique geneid in GE.txt\n")

# Run
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/01_build_inputs/make_geneloc_from_gtf.R