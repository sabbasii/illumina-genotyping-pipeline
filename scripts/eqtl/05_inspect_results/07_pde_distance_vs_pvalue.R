#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 07_pde_distance_vs_pvalue.R
#
# Scatter: signed SNP distance from gene TSS (kb) vs -log10(p),
# using eqtl_all results, restricted to PDE genes.
# Cis points are highlighted by color.
#
# Inputs (REPO_ROOT/output/eqtl):
#   - results/eqtl_all.tsv
#   - snpsloc.txt              (snpid, chr, pos)
#   - geneloc_extended.tsv     (geneid, chr, tss, strand, ...)
#
# Output:
#   - output/eqtl/results/inspect/pde/pde_dist_vs_logp_all.png
#
# Run:
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/05_inspect_results/07_pde_distance_vs_pvalue.R
#
# Optional:
#   --cis-dist-bp 1000000   # default 1 Mb
#   --pde-regex "^PDE"      # default "^PDE[0-9]"
#   --facet true|false      # default true
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  ix <- which(args == flag)
  if (length(ix) == 0) return(default)
  if (ix[1] == length(args)) return(default)
  args[ix[1] + 1]
}

cis_dist_bp <- as.numeric(get_arg("--cis-dist-bp", "1000000"))
pde_regex   <- get_arg("--pde-regex", "^PDE[0-9]")
facet_on    <- tolower(get_arg("--facet", "true")) %in% c("true","t","1","yes","y")

base_eqtl <- file.path(REPO_ROOT, "output", "eqtl")
res_file  <- file.path(base_eqtl, "results", "eqtl_all.tsv")
snploc_file <- file.path(base_eqtl, "snpsloc.txt")
genext_file <- file.path(base_eqtl, "geneloc_extended.tsv")

if (!file.exists(res_file)) stop("Missing: ", res_file)
if (!file.exists(snploc_file)) stop("Missing: ", snploc_file)
if (!file.exists(genext_file)) stop("Missing: ", genext_file)

out_dir <- file.path(base_eqtl, "results", "inspect", "pde")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(out_dir, "pde_dist_vs_logp_all.png")

# -------- load --------
x <- fread(res_file)
snploc <- fread(snploc_file)
genes <- fread(genext_file)

# standardize snpsloc: snpid chr pos -> SNP CHR POS
setnames(snploc, names(snploc)[1:3], c("SNP", "CHR_SNP", "POS_SNP"))

# gene annotation must include geneid, chr, tss
need_gene <- c("geneid","chr","tss")
miss_gene <- setdiff(need_gene, names(genes))
if (length(miss_gene) > 0) stop("geneloc_extended.tsv missing: ", paste(miss_gene, collapse=", "))

setnames(genes, "geneid", "gene")
setnames(genes, "chr", "CHR_GENE")

# p-value column (MatrixEQTL sometimes uses p-value)
p_col <- if ("p-value" %in% names(x)) "p-value" else if ("p.value" %in% names(x)) "p.value" else NA
if (is.na(p_col)) stop("No p-value column found (expected 'p-value' or 'p.value').")

# beta column optional (for effect-direction styling later)
beta_candidates <- c("beta","Beta","BETA","slope","Slope","effect","Effect")
beta_col <- beta_candidates[beta_candidates %in% names(x)][1]
if (is.na(beta_col)) beta_col <- NULL

# -------- filter to PDE genes --------
if (!("gene" %in% names(x))) stop("Association table missing 'gene' column.")
if (!("SNP" %in% names(x))) stop("Association table missing 'SNP' column.")

x <- x[grepl(pde_regex, gene)]
if (nrow(x) == 0) stop("No PDE rows found with regex: ", pde_regex)

x[, P := suppressWarnings(as.numeric(get(p_col)))]
x <- x[is.finite(P) & P > 0 & P <= 1]
x[, LOGP := -log10(P)]
if (!is.null(beta_col)) x[, beta := suppressWarnings(as.numeric(get(beta_col)))]

# -------- merge SNP + gene positions --------
df <- merge(x, snploc, by = "SNP", all.x = TRUE)
df <- merge(df, genes, by = "gene", all.x = TRUE)

# normalize chr formats (handle chr1 vs 1)
df[, CHR_SNP := gsub("^chr", "", as.character(CHR_SNP))]
df[, CHR_GENE := gsub("^chr", "", as.character(CHR_GENE))]
df[, POS_SNP := suppressWarnings(as.integer(POS_SNP))]
df[, tss := suppressWarnings(as.integer(tss))]

df <- df[is.finite(POS_SNP) & is.finite(tss) & !is.na(CHR_SNP) & !is.na(CHR_GENE)]
df <- df[CHR_SNP == CHR_GENE]  # distance only meaningful on same chromosome

# signed distance (kb)
df[, dist_bp := POS_SNP - tss]
df[, dist_kb := dist_bp / 1000]

# cis flag (default: within cis_dist_bp)
df[, is_cis := abs(dist_bp) <= cis_dist_bp]
df[, class := ifelse(is_cis, "cis", "other")]

# optional: keep a reasonable window for plotting clarity (still “all” results, just readable)
# comment out if you truly want everything
df_plot <- copy(df)
# df_plot <- df_plot[abs(dist_bp) <= (5 * cis_dist_bp)]

# -------- plot --------
p <- ggplot(df_plot, aes(x = dist_kb, y = LOGP, color = class)) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  labs(
    title = "PDE eQTL architecture: distance to TSS vs strength (all results)",
    subtitle = paste0("Cis = |distance| ≤ ", format(cis_dist_bp, big.mark=","), " bp; PDE filter: ", pde_regex),
    x = "Signed distance from TSS (kb)",
    y = expression(-log[10](p))
  ) +
  scale_color_manual(values = c(cis = "firebrick3", other = "grey55")) +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank())

if (facet_on) {
  # If you have many PDE genes, faceting can get huge.
  # You can restrict to top N genes by #associations if needed later.
  p <- p + facet_wrap(~ gene, scales = "free_x")
}

ggsave(out_png, p, width = if (facet_on) 14 else 8.5, height = if (facet_on) 9 else 5.5, dpi = 300)
cat("Saved: ", out_png, "\n", sep = "")
