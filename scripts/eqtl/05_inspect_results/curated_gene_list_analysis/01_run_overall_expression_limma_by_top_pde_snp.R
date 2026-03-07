#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ============================================================
# 01_run_overall_expression_limma_by_top_pde_snp.R
#
# Purpose
#   For a given PDE gene (e.g. PDE3A), identify its top eQTL SNP
#   and test whether expression of a curated gene set differs
#   between genotype groups using limma.
#
# Grouping logic
#   G01 = genotype 0 or 1
#   G2  = genotype 2
#
# Contrast tested
#   G2 - G01
#
# How to run
#
#   source scripts/00_config.sh
#
#   Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/01_run_overall_expression_limma_by_top_pde_snp.R \
#       --gene PDE3A
#
# Optional arguments
#
#   --gene-set input_data/target_lists/stroke_inflammation_signature_human.txt
#   --top 12
#
# Output
#
#   output/eqtl/results/inspect/pde/overall_expression_by_top_pde_snp/<PDE_GENE>/
#
#       <PDE_GENE>_<TOP_SNP>_overall_expression_limma.tsv
# ============================================================

DEFAULT_GENE_SET <- file.path(
  "input_data",
  "target_lists",
  "stroke_inflammation_signature_human.txt"
)

DEFAULT_TOP_N <- 12L

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default=NULL){
  w <- which(args == flag)
  if(length(w)==0) return(default)
  if(w[1]==length(args)) stop(paste0("Missing value after ",flag))
  args[w[1]+1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if(REPO_ROOT=="") stop("REPO_ROOT not set. Source scripts/00_config.sh first.")

gene_in <- get_arg("--gene","")
if(!nzchar(gene_in)) stop("Provide --gene (e.g. PDE3A)")

gene_set_path <- get_arg("--gene-set",DEFAULT_GENE_SET)

if(!grepl("^/",gene_set_path))
  gene_set_path <- file.path(REPO_ROOT,gene_set_path)

TOP_N <- suppressWarnings(as.integer(get_arg("--top",DEFAULT_TOP_N)))
if(!is.finite(TOP_N) || TOP_N<=0) TOP_N <- DEFAULT_TOP_N

PDE_GENE <- toupper(gene_in)

# ------------------------------------------------------------
# input paths
# ------------------------------------------------------------

eqtl_pde_file <- file.path(
  REPO_ROOT,
  "output/eqtl/results/inspect/pde/eqtl_all_pde.tsv"
)

SNP_file <- file.path(REPO_ROOT,"output/eqtl/SNP_overlap.txt")
GE_file  <- file.path(REPO_ROOT,"output/eqtl/GE_overlap.txt")

for(f in c(eqtl_pde_file,SNP_file,GE_file,gene_set_path)){
  if(!file.exists(f)) stop("Missing file: ",f)
}

# ------------------------------------------------------------
# identify top SNP
# ------------------------------------------------------------

hits <- fread(eqtl_pde_file)

p_col <- NULL
for(cand in c("pvalue","p-value","p.value","PValue","P","p")){
  if(cand %in% names(hits)){p_col<-cand;break}
}

if(is.null(p_col)) stop("No p-value column found")

hits[,geneU:=toupper(gene)]
hits[,P:=as.numeric(get(p_col))]
hits <- hits[is.finite(P)]

sub <- hits[geneU==PDE_GENE]

if(nrow(sub)==0)
  stop("No PDE hits found for gene: ",PDE_GENE)

setorder(sub,P)

top_snp <- sub$SNP[1]
top_p   <- sub$P[1]

# ------------------------------------------------------------
# output directory
# ------------------------------------------------------------

out_dir <- file.path(
  REPO_ROOT,
  "output/eqtl/results/inspect/pde",
  "overall_expression_by_top_pde_snp",
  PDE_GENE
)

dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

# ------------------------------------------------------------
# load matrices
# ------------------------------------------------------------

SNP <- fread(SNP_file)
GE  <- fread(GE_file)

setnames(SNP,"snpid","SNP")

snp_row <- SNP[SNP==top_snp]

if(nrow(snp_row)==0)
  stop("Top SNP not found in SNP matrix")

gs <- fread(gene_set_path,header=FALSE)[[1]]
gs <- unique(toupper(trimws(gs)))

GE[,geneU:=toupper(geneid)]
GEc <- GE[geneU %in% gs]

if(nrow(GEc)==0)
  stop("None of curated genes found in expression matrix")

# ------------------------------------------------------------
# sample overlap
# ------------------------------------------------------------

snp_samples <- names(SNP)[-1]
ge_samples  <- names(GE)[!(names(GE) %in% c("geneid","geneU"))]

common_samples <- intersect(snp_samples,ge_samples)

g <- as.numeric(unlist(snp_row[,..common_samples]))

df_g <- data.table(
  sample=common_samples,
  genotype=g
)

df_g <- df_g[genotype %in% c(0,1,2)]

df_g[,group:=ifelse(genotype %in% c(0,1),"G01","G2")]

df_g[,group:=factor(group,levels=c("G01","G2"))]

if(length(unique(df_g$group))<2)
  stop("Both genotype groups required")

# ------------------------------------------------------------
# expression matrix
# ------------------------------------------------------------

keep_samples <- intersect(common_samples,df_g$sample)

GEc2 <- GEc[,c("geneid",keep_samples),with=FALSE]

expr_mat <- as.matrix(GEc2[,-1])
rownames(expr_mat) <- GEc2$geneid
mode(expr_mat) <- "numeric"

ok_gene <- apply(expr_mat,1,function(x){
  x<-x[is.finite(x)]
  length(x)>=3 && sd(x)>0
})

expr_mat <- expr_mat[ok_gene,,drop=FALSE]

# ------------------------------------------------------------
# limma model
# ------------------------------------------------------------

df_g2 <- df_g[match(colnames(expr_mat),sample)]

design <- model.matrix(~0+group,data=df_g2)
colnames(design) <- levels(df_g2$group)

fit <- lmFit(expr_mat,design)

contr <- makeContrasts(
  G2vsG01 = G2 - G01,
  levels=design
)

fit2 <- eBayes(
  contrasts.fit(fit,contr),
  trend=TRUE
)

tt <- topTable(
  fit2,
  coef="G2vsG01",
  number=Inf,
  sort.by="P"
)

out <- as.data.table(tt,keep.rownames="geneid")

out[,pde_gene:=PDE_GENE]
out[,top_snp:=top_snp]
out[,top_snp_p:=top_p]
out[,contrast:="G2 - G01"]

# ------------------------------------------------------------
# output file
# ------------------------------------------------------------

outfile <- file.path(
  out_dir,
  paste0(PDE_GENE,"_",top_snp,"_overall_expression_limma.tsv")
)

fwrite(out,outfile,sep="\t")

# ------------------------------------------------------------
# console summary
# ------------------------------------------------------------

cat("PDE gene: ",PDE_GENE,"\n",sep="")
cat("Top SNP:  ",top_snp," (p=",format(top_p,scientific=TRUE,digits=3),")\n",sep="")
cat("Grouping: G01=(0/1) vs G2=(2)\n")

cat("\nTop ",min(TOP_N,nrow(out))," hits:\n",sep="")
print(out[1:min(TOP_N,nrow(out)),.(geneid,logFC,P.Value,adj.P.Val)])

cat("\nOutput written to:\n",outfile,"\n")


# Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/01_run_overall_expression_limma_by_top_pde_snp.R \
#     --gene PDE7A