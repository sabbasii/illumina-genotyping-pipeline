#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 05_test_signature_score_by_top_pde_snp.R
#
# Purpose
#   Test whether genotype shifts the stroke inflammation
#   signature score.
#
# Signature score:
#   mean(zscore(expression of signature genes))
#
# Models tested:
#   1) Collapsed:  G2 vs G01
#   2) Additive:   signature_score ~ genotype
#
# Example run
# source scripts/00_config.sh
#
# Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/05_test_signature_score_by_top_pde_snp.R \
#   --gene PDE7A
#
# Output
#   <PDE>_<TOP_SNP>_signature_score_test.tsv
# ============================================================

args <- commandArgs(trailingOnly=TRUE)

get_arg <- function(flag, default=NULL){
  w <- which(args==flag)
  if(length(w)==0) return(default)
  args[w+1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if(REPO_ROOT=="") stop("Source scripts/00_config.sh first")

gene_in <- get_arg("--gene","")
if(gene_in=="") stop("Provide --gene")

PDE_GENE <- toupper(gene_in)

gene_set_file <- file.path(
REPO_ROOT,
"input_data/target_lists/stroke_inflammation_signature_human.txt"
)

eqtl_file <- file.path(
REPO_ROOT,
"output/eqtl/results/inspect/pde/eqtl_all_pde.tsv"
)

SNP_file <- file.path(REPO_ROOT,"output/eqtl/SNP_overlap.txt")
GE_file  <- file.path(REPO_ROOT,"output/eqtl/GE_overlap.txt")

# ------------------------------------------------------------
# find top SNP
# ------------------------------------------------------------

hits <- fread(eqtl_file)
hits[,geneU:=toupper(gene)]

p_col <- intersect(names(hits),
c("pvalue","p-value","p.value","PValue","P","p"))[1]

hits[,P:=as.numeric(get(p_col))]
hits <- hits[is.finite(P)]

sub <- hits[geneU==PDE_GENE]

setorder(sub,P)

top_snp <- sub$SNP[1]

# ------------------------------------------------------------
# load matrices
# ------------------------------------------------------------

SNP <- fread(SNP_file)
GE  <- fread(GE_file)

setnames(SNP,"snpid","SNP")

snp_row <- SNP[SNP==top_snp]

# ------------------------------------------------------------
# signature genes
# ------------------------------------------------------------

sig_genes <- fread(gene_set_file,header=FALSE)[[1]]
sig_genes <- unique(toupper(sig_genes))

GE[,geneU:=toupper(geneid)]
GEc <- GE[geneU %in% sig_genes]

# ------------------------------------------------------------
# sample overlap
# ------------------------------------------------------------

snp_samples <- names(SNP)[-1]
ge_samples <- names(GE)[!(names(GE) %in% c("geneid","geneU"))]

samples <- intersect(snp_samples,ge_samples)

# genotype

g <- as.numeric(unlist(snp_row[,..samples]))

df <- data.table(
sample=samples,
genotype=g
)

df <- df[genotype %in% c(0,1,2)]

# ------------------------------------------------------------
# expression matrix
# ------------------------------------------------------------

GEc2 <- GEc[,c("geneid",samples),with=FALSE]

expr <- as.matrix(GEc2[,-1])
rownames(expr) <- GEc2$geneid

mode(expr) <- "numeric"

# z-score genes

expr_z <- t(scale(t(expr)))

# signature score per sample

signature_score <- colMeans(expr_z,na.rm=TRUE)

df$signature_score <- signature_score[df$sample]

# collapsed groups

df[,group:=ifelse(genotype %in% c(0,1),"G01","G2")]

# ------------------------------------------------------------
# statistical tests
# ------------------------------------------------------------

# additive
m_add <- lm(signature_score ~ genotype, data=df)

# collapsed
m_group <- lm(signature_score ~ group, data=df)

res <- data.table(
PDE_gene = PDE_GENE,
top_snp = top_snp,

additive_beta = coef(m_add)["genotype"],
additive_p = summary(m_add)$coefficients["genotype","Pr(>|t|)"],

group_beta = coef(m_group)["groupG2"],
group_p = summary(m_group)$coefficients["groupG2","Pr(>|t|)"],

n = nrow(df)
)

# ------------------------------------------------------------
# output
# ------------------------------------------------------------

out_dir <- file.path(
REPO_ROOT,
"output/eqtl/results/inspect/pde",
"overall_expression_by_top_pde_snp",
PDE_GENE
)

dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

outfile <- file.path(
out_dir,
paste0(PDE_GENE,"_",top_snp,"_signature_score_test.tsv")
)

fwrite(res,outfile,sep="\t")

print(res)

cat("\nOutput written to:\n",outfile,"\n")