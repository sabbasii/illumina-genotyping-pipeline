#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
})

# ---- paths (REPO_ROOT) ----
REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

expr_path <- file.path(REPO_ROOT, "input_data/expr/expr_table.tsv")
gtf_path  <- file.path(REPO_ROOT, "reference/annotation/Homo_sapiens.GRCh37.87.gtf.gz")
out_dir   <- file.path(REPO_ROOT, "output/genotype_run1/eqtl")
out_path  <- file.path(out_dir, "geneloc.txt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load expression symbols ----
expr <- fread(expr_path, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("PROBEID", "SYMBOL", "GENENAME") %in% names(expr))) {
  stop("expr_table.tsv must have columns: PROBEID, SYMBOL, GENENAME")
}
symbols <- unique(expr$SYMBOL)
symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

# ---- import GTF (genes only) ----
gtf <- import(gtf_path)                 # Ensembl GTF; rtracklayer reads .gtf.gz fine
gtf_genes <- gtf[gtf$type == "gene"]

gene_name <- mcols(gtf_genes)$gene_name
chr_raw   <- as.character(seqnames(gtf_genes))
left      <- start(gtf_genes)
right     <- end(gtf_genes)

# Keep standard chromosomes and add "chr" prefix to match your example (chr1, chrX, chrMT)
keep_chr <- chr_raw %in% c(as.character(1:22), "X", "Y", "MT")
df <- data.frame(
  geneid = gene_name,
  chr    = paste0("chr", chr_raw),
  left   = left,
  right  = right,
  stringsAsFactors = FALSE
)
df <- df[keep_chr & !is.na(df$geneid) & nzchar(df$geneid), , drop = FALSE]

# ---- filter to genes present in expression + de-duplicate (keep widest span per geneid) ----
df <- df[df$geneid %in% symbols, , drop = FALSE]
df$span <- df$right - df$left

# If a gene appears multiple times, keep the entry with the largest span
df <- df[order(df$geneid, -df$span), ]
df <- df[!duplicated(df$geneid), c("geneid", "chr", "left", "right")]

# ---- write output ----
fwrite(df, out_path, sep = "\t", quote = FALSE)

cat("Wrote:", out_path, "\n")
cat("Mapped genes:", nrow(df), "of", length(symbols), "unique SYMBOLs\n")

# ---- Run ----
# source scripts/00_config.sh
# Rscript scripts/eQTL/make_geneloc_from_gtf.R