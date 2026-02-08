#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
})

# ------------------------------------------------------------
# make_geneloc_from_gtf.R
#
# Build Matrix-eQTL geneloc.txt from an Ensembl GTF, filtered to
# the gene IDs present in GE.txt (gene symbols).
#
# Inputs:
#   - output/eqtl/GE.txt
#   - reference/annotation/Homo_sapiens.GRCh37.87.gtf.gz
#
# Outputs (under output/eqtl/):
#   1) geneloc.txt
#      columns: geneid, chr, left, right
#      (kept stable for Matrix-eQTL usage)
#
#   2) geneloc_extended.tsv
#      columns: geneid, ensembl_gene_id, chr, left, right, strand, tss, gene_biotype
#      (richer annotation for inspection scripts; does not affect eqtl engine)
#
# Run:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/01_build_inputs/make_geneloc_from_gtf.R
# ------------------------------------------------------------

# ---- paths (REPO_ROOT) ----
REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

ge_path  <- file.path(REPO_ROOT, "output/eqtl/GE.txt")
gtf_path <- file.path(REPO_ROOT, "reference/annotation/Homo_sapiens.GRCh37.87.gtf.gz")

out_dir        <- file.path(REPO_ROOT, "output/eqtl")
out_path       <- file.path(out_dir, "geneloc.txt")
out_path_ext   <- file.path(out_dir, "geneloc_extended.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ge_path)) stop("Missing GE.txt: ", ge_path)
if (!file.exists(gtf_path)) stop("Missing GTF: ", gtf_path)

# ---- load gene IDs from GE.txt (authoritative) ----
ge <- fread(ge_path, sep = "\t", header = TRUE, data.table = FALSE, quote = "")
if (!("geneid" %in% names(ge))) stop("GE.txt must have a 'geneid' column.")

symbols <- unique(ge$geneid)
symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

# ---- import GTF (genes only) ----
gtf <- import(gtf_path)  # reads .gtf.gz
gtf_genes <- gtf[gtf$type == "gene"]

# Pull metadata
gene_name <- as.character(mcols(gtf_genes)$gene_name)
ensg      <- as.character(mcols(gtf_genes)$gene_id)
biotype   <- as.character(mcols(gtf_genes)$gene_biotype)
if (all(is.na(biotype))) biotype <- as.character(mcols(gtf_genes)$gene_type)  # some GTFs use gene_type

chr_raw <- as.character(seqnames(gtf_genes))
left    <- start(gtf_genes)
right   <- end(gtf_genes)
strandv <- as.character(strand(gtf_genes))  # "+" / "-" / "*"

# Keep standard chromosomes; add "chr" prefix (chr1..chr22, chrX, chrY, chrMT)
keep_chr <- chr_raw %in% c(as.character(1:22), "X", "Y", "MT")

dt <- data.table(
  geneid = gene_name,
  ensembl_gene_id = ensg,
  chr_raw = chr_raw,
  chr = paste0("chr", chr_raw),
  left = as.integer(left),
  right = as.integer(right),
  strand = strandv,
  gene_biotype = biotype
)

dt <- dt[keep_chr]
dt <- dt[!is.na(geneid) & nzchar(geneid)]

# ---- filter to genes present in GE ----
dt <- dt[geneid %in% symbols]

# ---- compute TSS (uses strand; falls back to left if strand missing) ----
# '+' => TSS = left, '-' => TSS = right, otherwise left
dt[, tss := left]
dt[strand %in% c("-", "-1"), tss := right]

# ---- de-duplicate (keep widest span per geneid) ----
dt[, span := right - left]
setorder(dt, geneid, -span)
dt <- dt[!duplicated(geneid)]

# ---- write minimal geneloc.txt (stable) ----
dt_min <- dt[, .(geneid, chr, left, right)]
fwrite(dt_min, out_path, sep = "\t", quote = FALSE, na = "NA")

# ---- write extended file (for inspection) ----
dt_ext <- dt[, .(geneid, ensembl_gene_id, chr, left, right, strand, tss, gene_biotype)]
fwrite(dt_ext, out_path_ext, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:\n")
cat("  ", out_path, "\n", sep = "")
cat("  ", out_path_ext, "\n", sep = "")
cat("Mapped genes:", nrow(dt_min), "of", length(symbols), "unique geneid in GE.txt\n")

# Run
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/01_build_inputs/make_geneloc_from_gtf.R
