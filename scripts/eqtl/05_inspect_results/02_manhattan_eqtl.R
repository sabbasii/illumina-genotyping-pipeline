#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# 02_manhattan_eqtl.R
#
# Manhattan plots for MatrixEQTL results.
#
# Modes:
#   - cis : plots all SNP–gene pairs from eqtl_cis.tsv (no collapsing)
#   - all : plots one point per SNP using min p-value across genes
#           (collapsed background) from eqtl_all.tsv
#
# Inputs (under REPO_ROOT/output/eqtl/):
#   - results/eqtl_cis.tsv
#   - results/eqtl_all.tsv
#   - snpsloc.txt   (snpid, chr, pos; chr may be "1" or "chr1")
#
# Outputs (under REPO_ROOT/output/eqtl/results/inspect/manhattan/):
#   - manhattan__cis.png
#   - manhattan__all_minp_per_snp.png
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

source(file.path(REPO_ROOT, "scripts/eqtl/utils/inspect_helpers.R"))
paths <- get_eqtl_paths()

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
mode <- "cis"
ix <- match("--mode", args)
if (!is.na(ix) && ix < length(args)) mode <- args[ix + 1]
mode <- tolower(mode)
if (!(mode %in% c("cis", "all"))) stop("Invalid --mode. Use: cis or all")

# ---- settings ----
line_suggestive <- -log10(1e-5)
line_genomewide <- -log10(5e-8)

# ---- quiet numeric parse (no warnings) ----
as_num_quiet <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  ok <- nzchar(x) & x != "NA"
  out <- rep(NA_real_, length(x))
  if (any(ok)) out[ok] <- as.numeric(x[ok])
  out
}

# ---- load inputs ----
snploc <- read_snpsloc(paths$snpsloc)  # SNP, CHR_NUM, POS (POS is numeric in updated helper)

res_file <- if (mode == "cis") paths$eqtl_cis else paths$eqtl_all
assert_file(res_file)
res <- fread(res_file)

if (!("SNP" %in% names(res))) stop("Association table missing column: SNP")

p_col <- pick_p_col(res)
if (!(p_col %in% names(res))) stop("P-value column not found: ", p_col)

# ---- prepare plot data (warning-free) ----
if (mode == "all") {
  # One point per SNP: keep the smallest valid p-value across genes for each SNP
  tmp <- res[, {
    p <- as_num_quiet(get(p_col))
    p <- p[is.finite(p) & p > 0 & p <= 1]
    .(P = if (length(p)) min(p) else NA_real_)
  }, by = SNP]
  tmp <- tmp[is.finite(P)]

  df <- merge(tmp, snploc, by = "SNP")  # drops SNPs missing from snpsloc (expected)
  df[, LOGP := -log10(P)]

  out_png <- file.path(paths$manhattan_dir, "manhattan__all_minp_per_snp.png")
  title_txt <- "Manhattan (min p per SNP)"
  pt_size <- 0.45
} else {
  # Cis: keep all SNP–gene pairs (no collapsing)
  tmp <- res[, .(SNP, P = as_num_quiet(get(p_col)))]
  tmp <- tmp[is.finite(P) & P > 0 & P <= 1]

  df <- merge(tmp, snploc, by = "SNP")
  df[, LOGP := -log10(P)]

  out_png <- file.path(paths$manhattan_dir, "manhattan__cis.png")
  title_txt <- "Cis eQTL Manhattan"
  pt_size <- 0.9
}

# ---- cumulative position (uses UPDATED helper add_cumpos: overflow-safe) ----
cum <- add_cumpos(df, chr_col = "CHR_NUM", pos_col = "POS")
df2 <- cum$df
axis_df <- cum$axis

# Alternate colors by chromosome (classic look)
df2[, COL := ifelse(CHR_NUM %% 2 == 0, "grey65", "black")]

p <- ggplot(df2, aes(BP_CUM, LOGP)) +
  geom_point(aes(color = COL), size = pt_size) +
  scale_color_identity() +
  geom_hline(yintercept = line_suggestive, color = "blue", linewidth = 0.5) +
  geom_hline(yintercept = line_genomewide, color = "red",  linewidth = 0.5) +
  scale_x_continuous(
    breaks = axis_df$CENTER,
    labels = axis_df$CHR_NUM,
    expand = c(0.01, 0.01)
  ) +
  labs(
    title = title_txt,
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8))
  )

ggsave(
  out_png, p,
  width  = if (mode == "all") 9 else 12,
  height = if (mode == "all") 3.2 else 4.2,
  dpi = 300
)
cat("Saved:", out_png, "\n")

# ---- Run ----
# source scripts/00_config.sh
# Rscript scripts/eqtl/05_inspect_results/02_manhattan_eqtl.R --mode cis
# Rscript scripts/eqtl/05_inspect_results/02_manhattan_eqtl.R --mode all
