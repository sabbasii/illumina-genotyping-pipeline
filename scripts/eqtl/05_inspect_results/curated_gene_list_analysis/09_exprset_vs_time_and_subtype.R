#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# scripts/eqtl/05_inspect_results/curated_gene_list_analysis/09_exprset_vs_time_and_subtype.R
#
# Purpose
#   Using a curated gene list, test whether the overall expression pattern across
#   those genes differs by:
#     1) stroke_to_draw_hours (continuous)
#     2) subtype              (categorical)
#
#   Overall expression is summarized via PCA on gene-wise z-scored expression.
#   Also supports a key extension:
#     - "Residualize out time before PCA" (gene-by-gene regression on time),
#       then re-run PCA and re-test/plot by subtype.
#
# Inputs
#   - output/eqtl/GE_overlap.txt      (genes x samples; first col 'geneid')
#   - metadata/clinical_data.csv      (must include: sample_name, diagnosis, subtype, stroke_to_draw_hours)
#   - curated gene list (1 gene per line)
#
# Filters
#   --gene-set <path>                (default input_data/target_lists/stroke_inflammation_signature_human.txt)
#   --diagnosis ischemic_stroke      (default; comma-separated ok)
#   --sex male|female                (optional)
#   --min_subtype_n 10               (default)
#   --kpcs 5                         (default; MANOVA uses up to kpcs PCs)
#
# Time outliers (Option A)
#   IQR rule on stroke_to_draw_hours: keep within [Q1-1.5*IQR, Q3+1.5*IQR]
#   Prints which samples were removed + writes time_outliers_removed.tsv
#
# Output (flat dir)
#   output/eqtl/results/inspect/pde/exprset_by_time_and_subtype/
#     pca_scores.tsv
#     time_models.tsv
#     time_outliers_removed.tsv
#     subtype_counts_all.tsv
#     subtype_models.tsv
#     plot_pc1_vs_hours.png
#     plot_pc1_by_subtype.png
#     plot_pca_pc1_pc2_by_subtype.png
#     plot_pca_pc1_pc3_by_subtype.png
#     plot_pca_pc2_pc3_by_subtype.png
#     plot_pca_pc3_pc4_by_subtype.png
#
#   + Time-residualized PCA outputs:
#     pca_scores_time_resid.tsv
#     subtype_models_time_resid.tsv
#     plot_time_resid_pca_pc1_pc2_by_subtype.png
#     plot_time_resid_pca_pc1_pc3_by_subtype.png
#     plot_time_resid_pca_pc2_pc3_by_subtype.png
#     plot_time_resid_pca_pc3_pc4_by_subtype.png
#
# How to run
#   source scripts/00_config.sh
#   Rscript scripts/eqtl/05_inspect_results/curated_gene_list_analysis/09_exprset_vs_time_and_subtype.R \
#     --diagnosis ischemic_stroke
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) stop(paste0("Missing value after ", flag))
  args[w[1] + 1]
}

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (!nzchar(REPO_ROOT)) stop("REPO_ROOT not set. Source scripts/00_config.sh first.")

DEFAULT_GENE_SET <- file.path("input_data", "target_lists", "stroke_inflammation_signature_human.txt")
gene_set_path <- get_arg("--gene-set", DEFAULT_GENE_SET)
if (!grepl("^/", gene_set_path)) gene_set_path <- file.path(REPO_ROOT, gene_set_path)

dx_in   <- trimws(get_arg("--diagnosis", "ischemic_stroke"))
sex_in  <- tolower(trimws(get_arg("--sex", "")))
minN    <- suppressWarnings(as.integer(get_arg("--min_subtype_n", "10")))
kpcs    <- suppressWarnings(as.integer(get_arg("--kpcs", "5")))
if (!is.finite(minN) || minN < 2) minN <- 10
if (!is.finite(kpcs) || kpcs < 2) kpcs <- 5
if (nzchar(sex_in) && !(sex_in %in% c("male", "female"))) stop("--sex must be 'male' or 'female' (or omit).")

dx_keep <- character()
if (nzchar(dx_in)) {
  dx_keep <- tolower(trimws(unlist(strsplit(dx_in, ",", fixed = TRUE))))
  dx_keep <- dx_keep[nzchar(dx_keep)]
}

GE_file       <- file.path(REPO_ROOT, "output/eqtl/GE_overlap.txt")
clinical_file <- file.path(REPO_ROOT, "metadata/clinical_data.csv")
for (f in c(GE_file, clinical_file, gene_set_path)) {
  if (!file.exists(f)) stop("Missing file: ", f)
}

fmt_int <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "d", big.mark = ","))
fmt_num <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_p   <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))

# -----------------------
# Load gene list
# -----------------------
genes <- fread(gene_set_path, header = FALSE)[[1]]
genes <- unique(toupper(trimws(genes)))
genes <- genes[nzchar(genes)]
if (!length(genes)) stop("Gene list is empty: ", gene_set_path)

# -----------------------
# Load expression matrix
# -----------------------
GE <- fread(GE_file)
if (!("geneid" %in% names(GE))) stop("GE_overlap.txt must have first column named 'geneid'")

GE[, geneU := toupper(geneid)]
GEc <- GE[geneU %in% genes]
if (nrow(GEc) < 3) stop("Too few curated genes found in GE_overlap.txt (found ", nrow(GEc), ").")

sample_cols <- setdiff(names(GE), c("geneid", "geneU"))
GEc2 <- GEc[, c("geneid", sample_cols), with = FALSE]

expr_mat <- as.matrix(GEc2[, -"geneid"])
rownames(expr_mat) <- GEc2$geneid
mode(expr_mat) <- "numeric"

# Gene QC
ok_gene <- apply(expr_mat, 1, function(x) {
  x <- x[is.finite(x)]
  length(x) >= 10 && stats::sd(x) > 0
})
expr_mat <- expr_mat[ok_gene, , drop = FALSE]
if (nrow(expr_mat) < 3) stop("Too few usable genes after QC (need >=3).")

# Gene-wise z-score (rows=genes)
z_gene <- t(scale(t(expr_mat)))
z_gene[!is.finite(z_gene)] <- NA_real_

mean_z <- colMeans(z_gene, na.rm = TRUE)

# PCA on samples (samples x genes)
pca <- prcomp(t(z_gene), center = FALSE, scale. = FALSE)
pc_scores <- as.data.table(pca$x)
pc_scores[, sample := rownames(pca$x)]
setcolorder(pc_scores, c("sample", setdiff(names(pc_scores), "sample")))

# Ensure we have enough PCs for plotting (PC4) regardless of kpcs
keep_pcs_main <- paste0("PC", seq_len(min(max(kpcs, 4), ncol(pca$x))))
keep_pcs_main <- keep_pcs_main[keep_pcs_main %in% names(pc_scores)]
pc_scores <- pc_scores[, c("sample", keep_pcs_main), with = FALSE]
pc_scores[, mean_z := as.numeric(mean_z[sample])]

# -----------------------
# Load clinical + harmonize sample IDs
# -----------------------
clin <- fread(clinical_file)
need_cols <- c("sample_name","sex","diagnosis","subtype","stroke_to_draw_hours","age","ancestry")
for (nm in need_cols) if (!(nm %in% names(clin))) stop("clinical_data.csv missing required column: ", nm)

clin[, sample_raw  := as.character(sample_name)]
clin[, sample_dash := gsub("_", "-", sample_raw)]
clin[, sample_usg  := gsub("-", "_", sample_raw)]

m_raw  <- sum(clin$sample_raw  %chin% pc_scores$sample)
m_dash <- sum(clin$sample_dash %chin% pc_scores$sample)
m_usg  <- sum(clin$sample_usg  %chin% pc_scores$sample)

if (max(m_raw, m_dash, m_usg) == 0) {
  stop(
    "No overlap between clinical sample_name and expression sample columns.\n",
    "Example clinical: ", clin$sample_name[1], "\n",
    "Example expr sample: ", pc_scores$sample[1]
  )
}

if (m_dash >= m_raw && m_dash >= m_usg) {
  clin[, sample := sample_dash]
} else if (m_raw >= m_usg) {
  clin[, sample := sample_raw]
} else {
  clin[, sample := sample_usg]
}

clin[, sex_l := tolower(trimws(as.character(sex)))]
clin[, diagnosis_l := tolower(trimws(as.character(diagnosis)))]
clin[, subtype_l := tolower(trimws(as.character(subtype)))]
clin[, ancestry_l := tolower(trimws(as.character(ancestry)))]
clin[, stroke_to_draw_hours := suppressWarnings(as.numeric(stroke_to_draw_hours))]

# Filter clinical first (avoid any merge PC-name collisions later)
clin_f <- copy(clin)
if (nzchar(sex_in)) clin_f <- clin_f[sex_l == sex_in]
if (length(dx_keep)) clin_f <- clin_f[diagnosis_l %in% dx_keep]

# Merge: PCs + filtered clinical (PC names exist only in pc_scores)
df <- merge(pc_scores, clin_f, by = "sample", all.x = TRUE)
df <- df[!is.na(sample_name)]  # keep only samples with clinical match
if (nrow(df) < 30) stop("Too few samples after filtering (n=", nrow(df), ").")

# -----------------------
# Output directory
# -----------------------
out_dir <- file.path(REPO_ROOT, "output/eqtl/results/inspect/pde/exprset_by_time_and_subtype")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save PCA scores + key clinical columns
keep_cols <- c("sample", keep_pcs_main, "mean_z", "diagnosis", "subtype", "stroke_to_draw_hours", "sex", "age", "ancestry")
keep_cols <- keep_cols[keep_cols %in% names(df)]
fwrite(df[, ..keep_cols], file.path(out_dir, "pca_scores.tsv"), sep = "\t")

# -----------------------
# Analysis A: PC1 vs stroke_to_draw_hours (with IQR outlier removal)
# -----------------------
df_time <- df[is.finite(stroke_to_draw_hours)]

time_models <- data.table(
  term = character(),
  outcome = character(),
  beta = numeric(),
  p = numeric(),
  n = integer()
)

df_time_clean <- copy(df_time)
time_outliers <- df_time[0]

if (nrow(df_time) >= 30 && ("PC1" %in% names(df_time))) {
  x <- df_time$stroke_to_draw_hours
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qs[[2]] - qs[[1]]
  lo <- qs[[1]] - 1.5 * iqr
  hi <- qs[[2]] + 1.5 * iqr

  df_time_clean <- df_time[stroke_to_draw_hours >= lo & stroke_to_draw_hours <= hi]
  time_outliers <- df_time[!(sample %chin% df_time_clean$sample)]

  cat(sprintf("\nTime outlier filter (IQR): kept %d / %d (lo=%0.3f, hi=%0.3f)\n",
              nrow(df_time_clean), nrow(df_time), lo, hi))
  if (nrow(time_outliers)) {
    cat("Removed stroke_to_draw_hours outlier sample(s):\n")
    print(as.data.frame(time_outliers[, .(sample, stroke_to_draw_hours, subtype, diagnosis)]), row.names = FALSE)
  } else {
    cat("Removed stroke_to_draw_hours outlier sample(s): none\n")
  }

  fit1 <- lm(PC1 ~ stroke_to_draw_hours, data = df_time_clean)
  fit2 <- lm(PC1 ~ log1p(stroke_to_draw_hours), data = df_time_clean)

  sm1 <- summary(fit1)$coefficients
  sm2 <- summary(fit2)$coefficients

  time_models <- rbind(
    time_models,
    data.table(
      term = "stroke_to_draw_hours",
      outcome = "PC1",
      beta = sm1["stroke_to_draw_hours", "Estimate"],
      p = sm1["stroke_to_draw_hours", "Pr(>|t|)"],
      n = nrow(model.frame(fit1))
    ),
    data.table(
      term = "log1p(stroke_to_draw_hours)",
      outcome = "PC1",
      beta = sm2["log1p(stroke_to_draw_hours)", "Estimate"],
      p = sm2["log1p(stroke_to_draw_hours)", "Pr(>|t|)"],
      n = nrow(model.frame(fit2))
    )
  )

  png(file.path(out_dir, "plot_pc1_vs_hours.png"), width = 1200, height = 900, res = 150)
  plot(df_time_clean$stroke_to_draw_hours, df_time_clean$PC1,
       xlab = "stroke_to_draw_hours (IQR-filtered)",
       ylab = "PC1 (overall curated-list expression)",
       main = "PC1 vs stroke_to_draw_hours", pch = 16)
  abline(lm(PC1 ~ stroke_to_draw_hours, data = df_time_clean), lwd = 2)
  dev.off()
} else {
  cat("\nTime models: skipped (not enough samples with stroke_to_draw_hours and PC1).\n")
}

fwrite(time_models, file.path(out_dir, "time_models.tsv"), sep = "\t")
if (nrow(time_outliers)) {
  fwrite(time_outliers[, .(sample, stroke_to_draw_hours, subtype, diagnosis)],
         file.path(out_dir, "time_outliers_removed.tsv"), sep = "\t")
} else {
  fwrite(data.table(sample=character(), stroke_to_draw_hours=numeric(), subtype=character(), diagnosis=character()),
         file.path(out_dir, "time_outliers_removed.tsv"), sep = "\t")
}

# -----------------------
# Analysis B: subtype (counts + ANOVA + MANOVA)
# -----------------------
df_sub <- df[nzchar(subtype_l)]
tab <- df_sub[, .N, by = subtype_l][order(-N)]
keep_subtypes <- tab[N >= minN]$subtype_l
df_sub <- df_sub[subtype_l %in% keep_subtypes]
df_sub[, subtype_l := factor(subtype_l, levels = keep_subtypes)]

subtype_models <- data.table(model=character(), stat=character(), value=numeric(), p=numeric(), n=integer(), n_subtypes=integer())

if (nrow(df_sub) >= 30 && length(levels(df_sub$subtype_l)) >= 2 && ("PC1" %in% names(df_sub))) {
  aov_fit <- aov(PC1 ~ subtype_l, data = df_sub)
  aov_tab <- summary(aov_fit)[[1]]
  subtype_models <- rbind(
    subtype_models,
    data.table(
      model = "ANOVA",
      stat  = "F",
      value = as.numeric(aov_tab["subtype_l","F value"]),
      p     = as.numeric(aov_tab["subtype_l","Pr(>F)"]),
      n     = nrow(model.frame(aov_fit)),
      n_subtypes = length(levels(df_sub$subtype_l))
    )
  )

  pcs_for_man <- intersect(paste0("PC", seq_len(min(kpcs, ncol(pca$x)))), names(df_sub))
  if (length(pcs_for_man) >= 2) {
    Y <- as.matrix(df_sub[, ..pcs_for_man])
    man <- manova(Y ~ subtype_l, data = df_sub)
    man_sum <- summary(man, test = "Pillai")$stats
    subtype_models <- rbind(
      subtype_models,
      data.table(
        model = "MANOVA",
        stat  = "Pillai",
        value = as.numeric(man_sum["subtype_l","approx F"]),
        p     = as.numeric(man_sum["subtype_l","Pr(>F)"]),
        n     = nrow(df_sub),
        n_subtypes = length(levels(df_sub$subtype_l))
      )
    )
  }

  png(file.path(out_dir, "plot_pc1_by_subtype.png"), width = 1600, height = 900, res = 150)
  boxplot(PC1 ~ subtype_l, data = df_sub, las = 2,
          xlab = sprintf("Subtype (kept if n>=%d)", minN),
          ylab = "PC1 (overall curated-list expression)",
          main = sprintf("PC1 by subtype (min_n=%d)", minN))
  dev.off()
} else {
  cat("\nSubtype models: skipped (not enough samples/subtypes after min_subtype_n filter).\n")
}

fwrite(tab, file.path(out_dir, "subtype_counts_all.tsv"), sep = "\t")
fwrite(subtype_models, file.path(out_dir, "subtype_models.tsv"), sep = "\t")

# -----------------------
# Plot: PCA pairs colored by subtype with % variance explained
# -----------------------
plot_pca_pair <- function(d, pcx, pcy, pca_obj, out_png, title_prefix = "PCA of curated gene-set expression by stroke subtype") {
  if (!(pcx %in% names(d)) || !(pcy %in% names(d))) return(invisible(FALSE))
  dd <- d[is.finite(get(pcx)) & is.finite(get(pcy)) & !is.na(subtype_l)]
  if (nrow(dd) < 3) return(invisible(FALSE))

  dd[, subtype_f := droplevels(as.factor(subtype_l))]
  cols <- as.numeric(dd$subtype_f)

  var_expl <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  ix <- as.integer(sub("^PC", "", pcx))
  iy <- as.integer(sub("^PC", "", pcy))
  vx <- if (is.finite(ix) && ix <= length(var_expl)) 100 * var_expl[ix] else NA_real_
  vy <- if (is.finite(iy) && iy <= length(var_expl)) 100 * var_expl[iy] else NA_real_

  xlab <- if (is.finite(vx)) sprintf("%s (%.1f%% variance)", pcx, vx) else pcx
  ylab <- if (is.finite(vy)) sprintf("%s (%.1f%% variance)", pcy, vy) else pcy

  ttl <- title_prefix
  if (is.finite(vx) && is.finite(vy)) {
    ttl <- sprintf("%s\n%s+%s = %.1f%% variance", title_prefix, pcx, pcy, vx + vy)
  }

  png(out_png, width = 1200, height = 900, res = 150)
  plot(dd[[pcx]], dd[[pcy]],
       col = cols, pch = 16,
       xlab = xlab, ylab = ylab,
       main = ttl)
  legend("topright",
         legend = levels(dd$subtype_f),
         col = seq_along(levels(dd$subtype_f)),
         pch = 16, cex = 0.9)
  dev.off()

  invisible(TRUE)
}

if (nrow(df_sub) > 0) {
  plot_pca_pair(df_sub, "PC1", "PC2", pca, file.path(out_dir, "plot_pca_pc1_pc2_by_subtype.png"))
  plot_pca_pair(df_sub, "PC1", "PC3", pca, file.path(out_dir, "plot_pca_pc1_pc3_by_subtype.png"))
  plot_pca_pair(df_sub, "PC2", "PC3", pca, file.path(out_dir, "plot_pca_pc2_pc3_by_subtype.png"))
  plot_pca_pair(df_sub, "PC3", "PC4", pca, file.path(out_dir, "plot_pca_pc3_pc4_by_subtype.png"))
}

# -----------------------
# Residualize time BEFORE PCA (gene-by-gene), then re-test/plot by subtype
#   (Avoids PC1.x/PC1.y by NEVER merging two tables that both contain PC columns.)
# -----------------------
ellipse_points <- function(mu, Sigma, level = 0.95, n = 120) {
  if (any(!is.finite(mu))) return(NULL)
  if (!all(dim(Sigma) == c(2,2))) return(NULL)
  ev <- eigen(Sigma, symmetric = TRUE)
  if (any(!is.finite(ev$values)) || any(ev$values <= 0)) return(NULL)
  r <- sqrt(stats::qchisq(level, df = 2))
  theta <- seq(0, 2*pi, length.out = n)
  circle <- rbind(cos(theta), sin(theta))
  A <- ev$vectors %*% diag(sqrt(ev$values), 2, 2)
  pts <- t(A %*% circle)
  pts[, 1] <- pts[, 1] + mu[1]
  pts[, 2] <- pts[, 2] + mu[2]
  pts
}

plot_pca_pair_with_centroids <- function(d, pcx, pcy, pca_obj, out_png,
                                        title_prefix = "Time-residualized PCA by stroke subtype",
                                        ellipse_level = 0.95) {
  if (!(pcx %in% names(d)) || !(pcy %in% names(d))) return(invisible(FALSE))
  dd <- d[is.finite(get(pcx)) & is.finite(get(pcy)) & !is.na(subtype_l)]
  if (nrow(dd) < 5) return(invisible(FALSE))

  dd[, subtype_f := droplevels(as.factor(subtype_l))]
  cols <- as.numeric(dd$subtype_f)

  var_expl <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  ix <- as.integer(sub("^PC", "", pcx))
  iy <- as.integer(sub("^PC", "", pcy))
  vx <- if (is.finite(ix) && ix <= length(var_expl)) 100 * var_expl[ix] else NA_real_
  vy <- if (is.finite(iy) && iy <= length(var_expl)) 100 * var_expl[iy] else NA_real_

  xlab <- if (is.finite(vx)) sprintf("%s (%.1f%% variance)", pcx, vx) else pcx
  ylab <- if (is.finite(vy)) sprintf("%s (%.1f%% variance)", pcy, vy) else pcy
  ttl <- title_prefix
  if (is.finite(vx) && is.finite(vy)) {
    ttl <- sprintf("%s\n%s+%s = %.1f%% variance", title_prefix, pcx, pcy, vx + vy)
  }

  png(out_png, width = 1200, height = 900, res = 150)
  plot(dd[[pcx]], dd[[pcy]],
       col = cols, pch = 16,
       xlab = xlab, ylab = ylab,
       main = ttl)

  for (lvl in levels(dd$subtype_f)) {
    g <- dd[subtype_f == lvl]
    if (nrow(g) >= 3) {
      mu <- c(mean(g[[pcx]]), mean(g[[pcy]]))
      Sigma <- stats::cov(cbind(g[[pcx]], g[[pcy]]))
      pts <- ellipse_points(mu, Sigma, level = ellipse_level, n = 120)
      col_i <- as.numeric(factor(lvl, levels(dd$subtype_f)))
      if (!is.null(pts)) lines(pts[,1], pts[,2], col = col_i, lwd = 2)
      points(mu[1], mu[2], pch = 4, cex = 1.6, lwd = 2, col = col_i)
      text(mu[1], mu[2], labels = lvl, pos = 3, cex = 0.8, col = col_i)
    }
  }

  legend("topright",
         legend = levels(dd$subtype_f),
         col = seq_along(levels(dd$subtype_f)),
         pch = 16, cex = 0.9)
  dev.off()

  invisible(TRUE)
}

# Only if we have enough time-clean samples
if (nrow(df_time_clean) >= 30) {

  # Use time-clean cohort
  time_samples <- df_time_clean$sample
  common_samps <- intersect(colnames(z_gene), time_samples)

  if (length(common_samps) >= 30) {

    # Keep order aligned to clinical/time data
    df_time_sub <- df_time_clean[sample %chin% common_samps]
    df_time_sub <- df_time_sub[match(common_samps, df_time_sub$sample)]
    stopifnot(all(df_time_sub$sample == common_samps))

    # Residualize each gene's z-scored expression on log1p(time)
    Zt <- z_gene[, common_samps, drop = FALSE]
    Tcov <- log1p(df_time_sub$stroke_to_draw_hours)

    resid_mat <- matrix(NA_real_, nrow = nrow(Zt), ncol = ncol(Zt),
                        dimnames = dimnames(Zt))

    for (i in seq_len(nrow(Zt))) {
      y <- as.numeric(Zt[i, ])
      ok <- is.finite(y) & is.finite(Tcov)
      if (sum(ok) >= 10 && stats::sd(y[ok]) > 0) {
        fit <- stats::lm(y[ok] ~ Tcov[ok])
        r <- rep(NA_real_, length(y))
        r[ok] <- stats::residuals(fit)
        resid_mat[i, ] <- r
      }
    }

    # Re-zscore residuals gene-wise
    resid_z <- t(scale(t(resid_mat)))
    resid_z[!is.finite(resid_z)] <- NA_real_

    ok_gene2 <- apply(resid_z, 1, function(x) {
      x <- x[is.finite(x)]
      length(x) >= 10 && stats::sd(x) > 0
    })
    resid_z <- resid_z[ok_gene2, , drop = FALSE]

    if (nrow(resid_z) >= 3) {

      # PCA on residualized expression
      pca_time_resid <- prcomp(t(resid_z), center = FALSE, scale. = FALSE)

      # keep up to max(kpcs,4) for plots
      keep_pcs_res <- paste0("PC", seq_len(min(max(kpcs, 4), ncol(pca_time_resid$x))))
      keep_pcs_res <- keep_pcs_res[keep_pcs_res %in% colnames(pca_time_resid$x)]

      pcR <- as.data.table(pca_time_resid$x)[, ..keep_pcs_res]
      pcR[, sample := rownames(pca_time_resid$x)]
      setcolorder(pcR, c("sample", keep_pcs_res))

      # Merge ONLY with clinical columns (no PC columns in clin_f), so no PC1.x/PC1.y ever happens
      clin_min <- clin_f[, .(sample, subtype_l, subtype, diagnosis, sex, age, ancestry)]
      dfR <- merge(pcR, clin_min, by = "sample", all.x = TRUE)

      # Apply same subtype minN filter
      dfR_sub <- dfR[nzchar(subtype_l)]
      tabR <- dfR_sub[, .N, by = subtype_l][order(-N)]
      keep_subR <- tabR[N >= minN]$subtype_l
      dfR_sub <- dfR_sub[subtype_l %in% keep_subR]
      dfR_sub[, subtype_l := factor(subtype_l, levels = keep_subR)]

      fwrite(dfR_sub, file.path(out_dir, "pca_scores_time_resid.tsv"), sep = "\t")

      # Subtype models on time-residualized PCs
      subtype_models_R <- data.table(model=character(), stat=character(), value=numeric(), p=numeric(), n=integer(), n_subtypes=integer())

      if (nrow(dfR_sub) >= 30 && length(levels(dfR_sub$subtype_l)) >= 2 && ("PC1" %in% names(dfR_sub))) {

        aovR <- aov(PC1 ~ subtype_l, data = dfR_sub)
        aov_tabR <- summary(aovR)[[1]]
        subtype_models_R <- rbind(
          subtype_models_R,
          data.table(
            model="ANOVA_time_resid",
            stat="F",
            value=as.numeric(aov_tabR["subtype_l","F value"]),
            p=as.numeric(aov_tabR["subtype_l","Pr(>F)"]),
            n=nrow(model.frame(aovR)),
            n_subtypes=length(levels(dfR_sub$subtype_l))
          )
        )

        pcs_for_man_R <- intersect(paste0("PC", seq_len(min(kpcs, ncol(pca_time_resid$x)))), names(dfR_sub))
        if (length(pcs_for_man_R) >= 2) {
          YR <- as.matrix(dfR_sub[, ..pcs_for_man_R])
          manR <- manova(YR ~ subtype_l, data = dfR_sub)
          man_sumR <- summary(manR, test = "Pillai")$stats
          subtype_models_R <- rbind(
            subtype_models_R,
            data.table(
              model="MANOVA_time_resid",
              stat="Pillai",
              value=as.numeric(man_sumR["subtype_l","approx F"]),
              p=as.numeric(man_sumR["subtype_l","Pr(>F)"]),
              n=nrow(dfR_sub),
              n_subtypes=length(levels(dfR_sub$subtype_l))
            )
          )
        }
      }

      fwrite(subtype_models_R, file.path(out_dir, "subtype_models_time_resid.tsv"), sep = "\t")

      # Plots with centroids + ellipses
      plot_pca_pair_with_centroids(dfR_sub, "PC1", "PC2", pca_time_resid,
                                   file.path(out_dir, "plot_time_resid_pca_pc1_pc2_by_subtype.png"))
      plot_pca_pair_with_centroids(dfR_sub, "PC1", "PC3", pca_time_resid,
                                   file.path(out_dir, "plot_time_resid_pca_pc1_pc3_by_subtype.png"))
      plot_pca_pair_with_centroids(dfR_sub, "PC2", "PC3", pca_time_resid,
                                   file.path(out_dir, "plot_time_resid_pca_pc2_pc3_by_subtype.png"))
      plot_pca_pair_with_centroids(dfR_sub, "PC3", "PC4", pca_time_resid,
                                   file.path(out_dir, "plot_time_resid_pca_pc3_pc4_by_subtype.png"))
    }
  }
}

# -----------------------
# Console summary
# -----------------------
cat("\n=== Curated gene-list overall expression: time + subtype ===\n")
cat(sprintf("Gene set: %s\n", basename(gene_set_path)))
cat(sprintf("Genes: requested=%s | found_in_GE=%s | used_after_QC=%s\n",
            fmt_int(length(genes)), fmt_int(nrow(GEc)), fmt_int(nrow(expr_mat))))
cat(sprintf("Filters: diagnosis=%s | sex=%s\n",
            if (length(dx_keep)) paste(dx_keep, collapse=",") else "all",
            if (nzchar(sex_in)) sex_in else "all"))
cat(sprintf("Samples analyzed (after filters + clinical match): n=%s\n", fmt_int(nrow(df))))
cat(sprintf("Time available (pre-filter): n=%s | post-outlier: n=%s\n",
            fmt_int(nrow(df_time)), fmt_int(nrow(df_time_clean))))
cat(sprintf("Subtype non-missing: n=%s | subtypes_kept(minN=%d): %s\n",
            fmt_int(nrow(df_sub)), minN, fmt_int(length(unique(df_sub$subtype_l)))))

if (nrow(time_models)) {
  cat("\nTime models (PC1 ~ time):\n")
  for (i in seq_len(nrow(time_models))) {
    cat(sprintf("  %s: beta=%s  p=%s  n=%s\n",
                time_models$outcome[i],
                fmt_num(time_models$beta[i], 4),
                fmt_p(time_models$p[i]),
                fmt_int(time_models$n[i])))
    cat(sprintf("    term=%s\n", time_models$term[i]))
  }
}

if (nrow(subtype_models)) {
  cat("\nSubtype models:\n")
  for (i in seq_len(nrow(subtype_models))) {
    cat(sprintf("  %s (%s): value=%s  p=%s  n=%s  k=%s\n",
                subtype_models$model[i],
                subtype_models$stat[i],
                fmt_num(subtype_models$value[i], 3),
                fmt_p(subtype_models$p[i]),
                fmt_int(subtype_models$n[i]),
                fmt_int(subtype_models$n_subtypes[i])))
  }
}

cat("\nSaved outputs in:\n  ", out_dir, "\n", sep = "")