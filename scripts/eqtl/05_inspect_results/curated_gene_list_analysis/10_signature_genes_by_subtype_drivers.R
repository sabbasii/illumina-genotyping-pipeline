#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 10_signature_genes_by_subtype_drivers.R
# Purpose:
#   Identify which genes in a curated signature differ by stroke subtype (etiology),
#   optionally adjusting for time-to-draw and demographics, and visualize drivers.
#
# Outputs (saved under):
#   output/eqtl/results/inspect/pde/exprset_by_time_and_subtype/
#     gene_subtype_tests.tsv
#     pairwise__<A>_vs_<B>.tsv
#     volcano__<A>_vs_<B>.png
#     topgenes_by_subtype_meanz.png
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  w <- which(args == flag)
  if (!length(w)) return(default)
  if (w[1] == length(args)) stop("Missing value after ", flag)
  args[w[1] + 1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (!nzchar(REPO_ROOT)) stop("REPO_ROOT not set. Source scripts/00_config.sh first.")

gene_set_path <- get_arg("--gene-set", file.path(REPO_ROOT, "input_data/target_lists/stroke_inflammation_signature_human.txt"))
if (!grepl("^/", gene_set_path)) gene_set_path <- file.path(REPO_ROOT, gene_set_path)

diagnosis_keep <- tolower(trimws(get_arg("--diagnosis", "ischemic_stroke")))
minN <- as.integer(get_arg("--min_subtype_n", "15"))
topN <- as.integer(get_arg("--topN", "40"))

GE_file <- file.path(REPO_ROOT, "output/eqtl/GE_overlap.txt")
clinical_file <- file.path(REPO_ROOT, "metadata/clinical_data.csv")
stopifnot(file.exists(GE_file), file.exists(clinical_file), file.exists(gene_set_path))

out_dir <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/exprset_by_time_and_subtype")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ---- Load gene list
genes <- unique(toupper(trimws(fread(gene_set_path, header=FALSE)[[1]])))
genes <- genes[nzchar(genes)]

# ---- Load GE
GE <- fread(GE_file)
stopifnot("geneid" %in% names(GE))
GE[, geneU := toupper(geneid)]
GEc <- GE[geneU %in% genes]
if (nrow(GEc) < 10) stop("Too few curated genes found in GE_overlap.txt: ", nrow(GEc))
sample_cols <- setdiff(names(GEc), c("geneid","geneU"))

# ---- Load clinical
clin <- fread(clinical_file)
need <- c("sample_name","diagnosis","subtype","stroke_to_draw_hours","age","sex","ancestry")
for (nm in need) if (!(nm %in% names(clin))) stop("clinical_data.csv missing: ", nm)

clin[, sample_raw := as.character(sample_name)]
clin[, sample_dash := gsub("_","-", sample_raw)]
clin[, sample_usg  := gsub("-","_", sample_raw)]

# pick best mapping
m_raw  <- sum(clin$sample_raw  %chin% sample_cols)
m_dash <- sum(clin$sample_dash %chin% sample_cols)
m_usg  <- sum(clin$sample_usg  %chin% sample_cols)
if (max(m_raw,m_dash,m_usg) == 0) stop("No overlap between clinical sample_name and GE samples.")
clin[, sample := if (m_dash >= m_raw && m_dash >= m_usg) sample_dash else if (m_raw >= m_usg) sample_raw else sample_usg]

clin[, diagnosis_l := tolower(trimws(as.character(diagnosis)))]
clin[, subtype_l   := tolower(trimws(as.character(subtype)))]
clin[, sex_l       := tolower(trimws(as.character(sex)))]
clin[, ancestry_l  := tolower(trimws(as.character(ancestry)))]
clin[, age := suppressWarnings(as.numeric(age))]
clin[, stroke_to_draw_hours := suppressWarnings(as.numeric(stroke_to_draw_hours))]
clin[, log_hours := log1p(stroke_to_draw_hours)]

# filter diagnosis
clin2 <- clin[diagnosis_l == diagnosis_keep]
if (nrow(clin2) < 50) stop("Too few samples after diagnosis filter: ", nrow(clin2))

# keep subtypes with enough n (focus on main etiologies)
sub_counts <- clin2[nzchar(subtype_l), .N, by=subtype_l][order(-N)]
keep_sub <- sub_counts[N >= minN]$subtype_l
clin2 <- clin2[subtype_l %in% keep_sub]
clin2[, subtype_l := droplevels(factor(subtype_l, levels=keep_sub))]

if (length(levels(clin2$subtype_l)) < 3) stop("Need >=3 subtypes after min_subtype_n filter.")

# ---- Build long expression table for curated genes
# Convert GEc (genes x samples) -> long: gene, sample, expr
long <- melt(
  GEc[, c("geneid", ..sample_cols), with=FALSE],
  id.vars = "geneid",
  variable.name = "sample",
  value.name = "expr"
)

long[, expr := suppressWarnings(as.numeric(expr))]
long <- long[is.finite(expr)]

# merge clinical
dt <- merge(long, clin2, by="sample", all.x=FALSE, all.y=FALSE)
if (nrow(dt) == 0) stop("No overlap after merging expression and clinical.")

# gene-wise z-score to compare effect sizes across genes
dt[, z_expr := as.numeric(scale(expr)), by = geneid]

# ---- Covariates: include only if informative
has_var <- function(v) {
  if (!(v %in% names(dt))) return(FALSE)
  x <- dt[[v]]
  if (is.numeric(x)) return(sum(is.finite(x)) >= 30 && stats::sd(x, na.rm=TRUE) > 0)
  x <- x[!is.na(x)]
  length(unique(x)) >= 2
}
cov_terms <- c()
if (has_var("log_hours")) cov_terms <- c(cov_terms, "log_hours")
if (has_var("age"))       cov_terms <- c(cov_terms, "age")
if (has_var("sex_l"))     cov_terms <- c(cov_terms, "sex_l")
if (has_var("ancestry_l"))cov_terms <- c(cov_terms, "ancestry_l")

rhs <- paste(c("subtype_l", cov_terms), collapse=" + ")
form <- as.formula(paste("z_expr ~", rhs))

# ---- Per-gene global subtype test (ANOVA on subtype term)
genes_u <- unique(dt$geneid)
res <- rbindlist(lapply(genes_u, function(g) {
  d <- dt[geneid == g]
  if (nrow(d) < 30) return(NULL)
  fit <- lm(form, data=d)
  a <- anova(fit)
  if (!("subtype_l" %in% rownames(a))) return(NULL)
  data.table(
    geneid = g,
    n = nrow(model.frame(fit)),
    p_subtype = a["subtype_l","Pr(>F)"]
  )
}), fill=TRUE)

if (!nrow(res)) stop("No genes produced subtype test results.")
res[, fdr_subtype := p.adjust(p_subtype, method="BH")]
setorder(res, fdr_subtype, p_subtype)

fwrite(res, file.path(out_dir, "gene_subtype_tests.tsv"), sep="\t")

# ---- Pairwise contrasts (difference in adjusted means using simple model coefficients)
# We'll do estimated mean difference from lm with subtype as factor (reference changes per contrast).
pairwise <- function(A, B) {
  dt2 <- dt[subtype_l %in% c(A,B)]
  dt2[, subtype_ab := droplevels(factor(subtype_l, levels=c(A,B)))]

  rhs2 <- paste(c("subtype_ab", cov_terms), collapse=" + ")
  form2 <- as.formula(paste("z_expr ~", rhs2))

  out <- rbindlist(lapply(genes_u, function(g) {
    d <- dt2[geneid == g]
    if (nrow(d) < 30) return(NULL)
    fit <- lm(form2, data=d)
    co <- summary(fit)$coefficients
    term <- "subtype_abB"
    # R will name it subtype_ab<levelB>, which is B
    term <- grep("^subtype_ab", rownames(co), value=TRUE)
    if (!length(term)) return(NULL)
    data.table(
      geneid = g,
      n = nrow(model.frame(fit)),
      beta = co[term[1],"Estimate"],
      p = co[term[1],"Pr(>|t|)"]
    )
  }), fill=TRUE)

  out[, fdr := p.adjust(p, method="BH")]
  out[, abs_beta := abs(beta)]
  out[, contrast := paste0(as.character(A), "_vs_", as.character(B))]
  setorder(out, fdr, p, -abs_beta)
  out
}

levs <- levels(clin2$subtype_l)
# focus on main three etiologies if present
# (your file has cardioembolic, small_vessel_disease, large_vessel_disease, cryptogenic)
want <- intersect(c("cardioembolic","small_vessel_disease","large_vessel_disease","cryptogenic"), levs)
if (length(want) >= 3) levs_use <- want else levs_use <- levs

pairs <- combn(levs_use, 2, simplify=FALSE)
for (pr in pairs) {
  A <- pr[[1]]; B <- pr[[2]]
  pw <- pairwise(A,B)
  if (!nrow(pw)) next
  out_tsv <- file.path(out_dir, sprintf("pairwise__%s_vs_%s.tsv", A, B))
  fwrite(pw, out_tsv, sep="\t")

  # Volcano plot (base R)
  png(file.path(out_dir, sprintf("volcano__%s_vs_%s.png", A, B)), width=1200, height=900, res=150)
  x <- pw$beta
  y <- -log10(pw$p)
  plot(x, y, pch=16, cex=0.7,
       xlab=sprintf("Effect (z-expression) %s - %s", B, A),
       ylab="-log10(p)",
       main=sprintf("Signature genes: %s vs %s", A, B))
  abline(h=-log10(0.05), lty=2)
  # label top genes by FDR then |beta|
  top <- pw[order(fdr, -abs(beta))][1:min(12, .N)]
  text(top$beta, -log10(top$p), labels=top$geneid, pos=3, cex=0.8)
  dev.off()
}

# ---- Heatmap-like plot: top genes by subtype mean z_expr
top_genes <- res[order(fdr_subtype, p_subtype)][1:min(topN, .N)]$geneid
means <- dt[geneid %in% top_genes, .(mean_z = mean(z_expr, na.rm=TRUE)), by=.(geneid, subtype_l)]
mat <- dcast(means, geneid ~ subtype_l, value.var="mean_z")
gene_order <- mat$geneid
m <- as.matrix(mat[, -1, with=FALSE])
rownames(m) <- gene_order

# reorder genes by subtype effect spread
spread <- apply(m, 1, function(v) max(v, na.rm=TRUE) - min(v, na.rm=TRUE))
ord <- order(spread, decreasing=TRUE)
m <- m[ord, , drop=FALSE]

png(file.path(out_dir, "topgenes_by_subtype_meanz.png"), width=1400, height=1000, res=150)
par(mar=c(8, 10, 3, 2))
image(t(m[nrow(m):1, , drop=FALSE]), axes=FALSE,
      xlab="", ylab="", main=sprintf("Top %d signature genes: mean z-expression by subtype", min(topN, length(top_genes))))
axis(1, at=seq(0,1,length.out=ncol(m)), labels=colnames(m), las=2)
axis(2, at=seq(0,1,length.out=nrow(m)), labels=rev(rownames(m)), las=2, cex.axis=0.6)
box()
dev.off()

cat("\nSaved driver-gene outputs in:\n  ", out_dir, "\n", sep="")
cat("Key files:\n")
cat("  gene_subtype_tests.tsv (global subtype p/FDR per gene)\n")
cat("  pairwise__A_vs_B.tsv + volcano__A_vs_B.png\n")
cat("  topgenes_by_subtype_meanz.png\n\n")