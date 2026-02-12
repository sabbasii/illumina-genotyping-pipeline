#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
})

# ------------------------------------------------------------
# 06_target_eqtl_network.R
#
# Build an SNP → gene network for a target gene set.
#
# Two input modes:
#   1) Independent (default):
#      Reads eqtl_{cis|all}.tsv and filters genes by:
#        --gene-regex and/or --gene-list
#
#   2) Reuse filtered hits:
#      Reads a pre-filtered hits file from step 05:
#        --hits-file <path/to/eqtl_<mode>_target_hits.tsv>
#
# The network shows:
#   - Nodes: SNPs and genes
#   - Edges: significant SNP–gene associations
#     edge color: sign(beta) (positive/negative)
#     edge width: -log10(p) or |beta|
#
# Output folder:
#   REPO_ROOT/output/eqtl/results/inspect/targets/<target_name>/<mode>/network/
#
# Files:
#   - eqtl_<mode>_target_network.png
#   - eqtl_<mode>_target_edges.tsv
#   - eqtl_<mode>_target_hub_snps.tsv
#   - eqtl_<mode>_target_hub_genes.tsv
# ------------------------------------------------------------

REPO_ROOT <- Sys.getenv("REPO_ROOT")
if (REPO_ROOT == "") stop("REPO_ROOT is not set. Did you source scripts/00_config.sh?")

source(file.path(REPO_ROOT, "scripts/eqtl/utils/inspect_helpers.R"))
paths <- get_eqtl_paths()

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) return(default)
  args[w[1] + 1]
}

mode <- tolower(get_arg("--mode", "cis"))
if (!(mode %in% c("cis", "all"))) stop("Invalid --mode. Use: cis or all")

target_name <- get_arg("--target-name", "target")
if (!nzchar(target_name)) target_name <- "target"

hits_file <- get_arg("--hits-file", "")
gene_regex <- get_arg("--gene-regex", "")
gene_list_file <- get_arg("--gene-list", "")

USE_FDR <- tolower(get_arg("--use-fdr", "true")) %in% c("true", "t", "1", "yes", "y")
FDR_CUTOFF <- suppressWarnings(as.numeric(get_arg("--fdr", "0.05")))
P_CUTOFF <- suppressWarnings(as.numeric(get_arg("--p", "1e-5")))

MAX_EDGES <- suppressWarnings(as.integer(get_arg("--max-edges", "2500")))
if (!is.finite(MAX_EDGES) || MAX_EDGES <= 0) MAX_EDGES <- 2500L

LABEL_TOP_SNPS  <- suppressWarnings(as.integer(get_arg("--label-snps", "12")))
LABEL_TOP_GENES <- suppressWarnings(as.integer(get_arg("--label-genes", "12")))
if (!is.finite(LABEL_TOP_SNPS)  || LABEL_TOP_SNPS  < 0) LABEL_TOP_SNPS  <- 12L
if (!is.finite(LABEL_TOP_GENES) || LABEL_TOP_GENES < 0) LABEL_TOP_GENES <- 12L

EDGE_WIDTH_MODE <- tolower(get_arg("--edge-width", "logp"))
if (!(EDGE_WIDTH_MODE %in% c("logp", "absbeta"))) stop("Invalid --edge-width. Use: logp or absbeta")

# ---- output dir ----
out_dir <- file.path(paths$inspect_dir, "targets", safe_name(target_name), mode, "network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png    <- file.path(out_dir, paste0("eqtl_", mode, "_target_network.png"))
out_edges  <- file.path(out_dir, paste0("eqtl_", mode, "_target_edges.tsv"))
out_hsnps  <- file.path(out_dir, paste0("eqtl_", mode, "_target_hub_snps.tsv"))
out_hgenes <- file.path(out_dir, paste0("eqtl_", mode, "_target_hub_genes.tsv"))

# ---- load edges ----
e <- NULL

if (nzchar(hits_file)) {
  assert_file(hits_file)
  dt <- fread(hits_file)

  # accept either (SNP,gene,P,FDRv,BETAv) or (SNP,gene,p-value,FDR,beta)
  if (!("SNP" %in% names(dt))) stop("hits-file missing column: SNP")
  if (!("gene" %in% names(dt))) stop("hits-file missing column: gene")

  # p
  p_col <- NULL
  if ("P" %in% names(dt)) p_col <- "P"
  if (is.null(p_col)) p_col <- pick_p_col(dt)

  # fdr
  fdr_col <- NULL
  if ("FDRv" %in% names(dt)) fdr_col <- "FDRv"
  if (is.null(fdr_col) && any(c("FDR","fdr","qvalue","q.value","q-value") %in% names(dt))) fdr_col <- pick_fdr_col(dt)

  # beta
  beta_col <- NULL
  if ("BETAv" %in% names(dt)) beta_col <- "BETAv"
  if (is.null(beta_col)) {
    beta_candidates <- c("beta", "Beta", "BETA", "slope", "Slope", "effect", "Effect")
    if (any(beta_candidates %in% names(dt))) beta_col <- pick_beta_col(dt)
  }

  dt[, p := suppressWarnings(as.numeric(get(p_col)))]
  dt <- dt[is.finite(p) & p > 0 & p <= 1]

  if (!is.null(fdr_col)) {
    dt[, FDR := suppressWarnings(as.numeric(get(fdr_col)))]
  } else {
    dt[, FDR := NA_real_]
  }

  if (!is.null(beta_col)) {
    dt[, beta := suppressWarnings(as.numeric(get(beta_col)))]
  } else {
    dt[, beta := NA_real_]
  }

  e <- dt[, .(SNP, gene, p, FDR, beta)]
} else {
  # independent mode: load full association results then filter
  res_file <- if (mode == "cis") paths$eqtl_cis else paths$eqtl_all
  assert_file(res_file)

  x <- fread(res_file)

  need_min <- c("SNP", "gene")
  miss_min <- setdiff(need_min, names(x))
  if (length(miss_min) > 0) stop("Association table missing columns: ", paste(miss_min, collapse = ", "))

  p_col <- pick_p_col(x)

  fdr_col <- NULL
  if (any(c("FDR", "fdr", "qvalue", "q.value", "q-value") %in% names(x))) {
    fdr_col <- pick_fdr_col(x)
  }

  beta_col <- NULL
  beta_candidates <- c("beta", "Beta", "BETA", "slope", "Slope", "effect", "Effect")
  if (any(beta_candidates %in% names(x))) beta_col <- pick_beta_col(x)

  x[, p := suppressWarnings(as.numeric(get(p_col)))]
  x <- x[is.finite(p) & p > 0 & p <= 1]

  if (!is.null(fdr_col)) {
    x[, FDR := suppressWarnings(as.numeric(get(fdr_col)))]
  } else {
    x[, FDR := NA_real_]
  }

  if (!is.null(beta_col)) {
    x[, beta := suppressWarnings(as.numeric(get(beta_col)))]
  } else {
    x[, beta := NA_real_]
  }

  # gene filter (regex/list required here)
  genes_keep <- NULL
  if (nzchar(gene_list_file)) {
    if (!file.exists(gene_list_file)) stop("Missing gene list: ", gene_list_file)
    gl <- fread(gene_list_file, header = FALSE, sep = "\n", data.table = FALSE)[, 1]
    gl <- trimws(gl)
    gl <- gl[nzchar(gl)]
    genes_keep <- unique(gl)
  }

  if (nzchar(gene_regex)) {
    if (is.null(genes_keep)) {
      x <- x[grepl(gene_regex, gene)]
    } else {
      x <- x[gene %in% genes_keep & grepl(gene_regex, gene)]
    }
  } else if (!is.null(genes_keep)) {
    x <- x[gene %in% genes_keep]
  } else {
    stop("Independent mode requires --gene-regex and/or --gene-list (or provide --hits-file).")
  }

  if (nrow(x) == 0) stop("No rows left after gene filtering.")

  e <- x[, .(SNP, gene, p, FDR, beta)]
}

# ---- significance filter ----
if (USE_FDR) {
  if (all(is.na(e$FDR))) stop("USE_FDR is true but no FDR column was found.")
  e <- e[is.finite(FDR) & FDR <= FDR_CUTOFF]
} else {
  e <- e[is.finite(p) & p <= P_CUTOFF]
}

e <- e[is.finite(p) & is.finite(beta)]
e <- unique(e, by = c("SNP", "gene"))

if (nrow(e) == 0) {
  stop("No edges after significance filtering. Try relaxing thresholds or switching --use-fdr false.")
}

# Keep strongest edges by p for readability
setorder(e, p)
if (nrow(e) > MAX_EDGES) e <- e[1:MAX_EDGES]

# ---- edge attributes ----
e[, sign := fifelse(beta >= 0, "positive", "negative")]
e[, logp := -log10(p)]
e[, absbeta := abs(beta)]

fwrite(e, out_edges, sep = "\t")

# ---- graph ----
nodes <- data.table(
  name = c(unique(e$SNP), unique(e$gene)),
  type = c(rep("SNP", length(unique(e$SNP))),
           rep("GENE", length(unique(e$gene))))
)

g <- graph_from_data_frame(
  d = e[, .(from = SNP, to = gene, p, FDR, beta, sign, logp, absbeta)],
  directed = TRUE,
  vertices = nodes
)

deg_all <- degree(g, mode = "all")
V(g)$degree <- as.integer(deg_all)

hub_snps <- data.table(
  node = V(g)$name[V(g)$type == "SNP"],
  degree = V(g)$degree[V(g)$type == "SNP"]
)[order(-degree)]

hub_genes <- data.table(
  node = V(g)$name[V(g)$type == "GENE"],
  degree = V(g)$degree[V(g)$type == "GENE"]
)[order(-degree)]

best_by_snp <- e[, .(best_p = min(p), max_absbeta = max(absbeta)), by = SNP]
best_by_gene<- e[, .(best_p = min(p), max_absbeta = max(absbeta)), by = gene]

hub_snps <- merge(hub_snps, best_by_snp, by.x="node", by.y="SNP", all.x=TRUE)
hub_genes<- merge(hub_genes, best_by_gene, by.x="node", by.y="gene", all.x=TRUE)

fwrite(hub_snps, out_hsnps, sep = "\t")
fwrite(hub_genes, out_hgenes, sep = "\t")

# ---- layout ----
set.seed(1)
lay <- layout_with_fr(g)

vl <- as.data.table(as_data_frame(g, what = "vertices"))
vl[, x := lay[, 1]]
vl[, y := lay[, 2]]

el <- as.data.table(as_data_frame(g, what = "edges"))

# attach coords
el <- merge(el, vl[, .(name, x, y)], by.x = "from", by.y = "name")
setnames(el, c("x", "y"), c("x_from", "y_from"))
el <- merge(el, vl[, .(name, x, y)], by.x = "to", by.y = "name")
setnames(el, c("x", "y"), c("x_to", "y_to"))

# width scaling
w <- if (EDGE_WIDTH_MODE == "absbeta") el$absbeta else el$logp
w <- w / max(w)
el[, edge_w := 0.25 + 1.8 * w]
el[, edge_a := 0.15 + 0.75 * w]

# labels: top hubs only
lab_snps  <- hub_snps[1:min(LABEL_TOP_SNPS, .N), node]
lab_genes <- hub_genes[1:min(LABEL_TOP_GENES, .N), node]

vl[, label := ""]
vl[type == "SNP"  & name %in% lab_snps,  label := name]
vl[type == "GENE" & name %in% lab_genes, label := name]

vl[, node_size := pmin(10, 2 + sqrt(degree))]

p <- ggplot() +
  geom_segment(
    data = el,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
        color = sign, alpha = edge_a, linewidth = edge_w),
    lineend = "round"
  ) +
  scale_color_manual(values = c(positive = "firebrick3", negative = "dodgerblue3")) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  geom_point(
    data = vl,
    aes(x = x, y = y, shape = type, size = node_size),
    color = "black", fill = "white", stroke = 0.6
  ) +
  scale_shape_manual(values = c(SNP = 21, GENE = 22)) +
  scale_size_identity() +
  geom_text(
    data = vl[label != ""],
    aes(x = x, y = y, label = label),
    size = 3, vjust = -0.8
  ) +
  labs(
    title = "Target SNP → gene eQTL network",
    subtitle = paste0(
      "Target: ", target_name, " | mode: ", mode, " | ",
      if (USE_FDR) paste0("FDR ≤ ", FDR_CUTOFF) else paste0("p ≤ ", format(P_CUTOFF, scientific = TRUE)),
      " | edge width = ", if (EDGE_WIDTH_MODE == "absbeta") "|beta|" else "-log10(p)",
      " | edges shown: ", nrow(e)
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(out_png, p, width = 11, height = 8, dpi = 300, bg = "white")

cat("Saved:\n  ", out_png, "\n  ", out_edges, "\n  ", out_hsnps, "\n  ", out_hgenes, "\n", sep = "")

# ---- Run ----
# source scripts/00_config.sh
#
# Independent (regex):
# Rscript scripts/eqtl/05_inspect_results/06_target_eqtl_network.R \
#   --mode all --target-name pde --gene-regex "^PDE[0-9]" --use-fdr true --fdr 0.05
#
# Independent (gene list):
# Rscript scripts/eqtl/05_inspect_results/06_target_eqtl_network.R \
#   --mode cis --target-name pde --gene-list input_data/target_lists/PDE_genes.txt --use-fdr false --p 1e-5
#
# Reuse step-05 hits file:
# Rscript scripts/eqtl/05_inspect_results/06_target_eqtl_network.R \
#   --mode all --target-name pde \
#   --hits-file output/eqtl/results/inspect/targets/pde/all/eqtl_all_target_hits.tsv
#
# Options:
#   --edge-width logp|absbeta
#   --max-edges 2500
#   --label-snps 12
#   --label-genes 12
