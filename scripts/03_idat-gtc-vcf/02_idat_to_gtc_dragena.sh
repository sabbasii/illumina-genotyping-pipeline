#!/usr/bin/env bash

# Only enable strict mode when executed directly (not when 'source'd into an interactive shell)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  set -euo pipefail
  IFS=$'\n\t'
fi

# export LC_ALL=C  # disabled to avoid prompt error

############################################
# Resolve repo root (portable, robust)
############################################
# We compute the repo root based on this file's location, and only keep an existing
# REPO_ROOT if it actually looks like the repo (i.e., contains scripts/00_config.sh).
if [[ -n "${BASH_VERSION:-}" ]]; then
  _SELF="${BASH_SOURCE[0]}"
elif [[ -n "${ZSH_VERSION:-}" ]]; then
  _SELF="${(%):-%x}"
else
  _SELF="$0"
fi

_CANDIDATE_ROOT="$(cd -- "$(dirname -- "$_SELF")/.." && pwd -P)"

# Override stale/incorrect REPO_ROOT values
if [[ -z "${REPO_ROOT:-}" || ! -f "$REPO_ROOT/scripts/00_config.sh" ]]; then
  REPO_ROOT="$_CANDIDATE_ROOT"
fi

unset _SELF _CANDIDATE_ROOT
export REPO_ROOT


############################################
# Run label and reference build
# `CHRSET_AUTOSOMES` → tells tools (e.g. plink2) which chromosome numbering scheme to expect.
# `SPLIT_PAR` → tells tools how to split pseudo-autosomal regions (PAR) on sex chromosomes, which differs between builds.
############################################
RUN="${RUN:-genotype_run1}"         # change to genotype_run2 etc. when you start a new run
REF_BUILD="${REF_BUILD:-GRCh37}"    # set to GRCh38 if needed
export RUN REF_BUILD

# Build-dependent constants for PLINK2
case "$REF_BUILD" in
  GRCh37) CHRSET_AUTOSOMES=37; SPLIT_PAR=b37 ;;
  GRCh38) CHRSET_AUTOSOMES=38; SPLIT_PAR=b38 ;;
  *) echo "[ERROR] Unknown REF_BUILD: $REF_BUILD" >&2; return 1 2>/dev/null || exit 1 ;;
esac
export CHRSET_AUTOSOMES SPLIT_PAR

############################################
# Inputs (relative to repo)
############################################
INPUT_DIR="$REPO_ROOT/input_data"

# subdirs
IDAT_DIR="$INPUT_DIR/idat"
MANIFEST_DIR="$INPUT_DIR/manifest"
SAMPLE_SHEET_DIR="$INPUT_DIR/sample_sheet"
CLUSTER_DIR="$INPUT_DIR/cluster"

# defaults (override by exporting these before sourcing this script)
: "${BPM_MANIFEST:=$MANIFEST_DIR/GSAMD-24v3-0-EA_20034606_A1.bpm}"
: "${CSV_MANIFEST:=$MANIFEST_DIR/GSAMD-24v3-0-EA_20034606_A1.csv}"
: "${EGT_CLUSTER:=$CLUSTER_DIR/GSA-24v3-0_A1_ClusterFile_custom.egt}"
: "${SAMPLE_SHEET:=$SAMPLE_SHEET_DIR/All_samples_Examine_SNPs_GWAS studies_GJ-P01_infiniumSampleSheet.csv}"
: "${SHEET_PATH:=$SAMPLE_SHEET}"

export INPUT_DIR IDAT_DIR MANIFEST_DIR SAMPLE_SHEET_DIR CLUSTER_DIR \
       BPM_MANIFEST CSV_MANIFEST EGT_CLUSTER SAMPLE_SHEET SHEET_PATH


############################################
# Reference genome
############################################
REF_DIR="$REPO_ROOT/reference/$REF_BUILD"
case "$REF_BUILD" in
  GRCh37) REFERENCE_FASTA="${REFERENCE_FASTA:-$REF_DIR/reference.fa}" ;;
  GRCh38) REFERENCE_FASTA="${REFERENCE_FASTA:-$REF_DIR/reference.fa}" ;;
esac
export REF_DIR REFERENCE_FASTA

############################################
# Outputs
############################################
OUT_DIR="${OUT_DIR:-$REPO_ROOT/output/$RUN}"
GTC_DIR="$OUT_DIR/gtc"
VCF_DIR="${VCF_DIR:-$OUT_DIR/vcf}"
CNV_DIR="$OUT_DIR/cnv"
QC_DIR="${QC_DIR:-$OUT_DIR/qc}"
LOG_DIR="$OUT_DIR/logs"

# PLINK working area and tmp
PLINK_DIR="${PLINK_DIR:-$QC_DIR/plink}"
TMP_DIR="${TMP_DIR:-$PLINK_DIR/plink_tmp}"

# QC sub-areas
QC_SUMMARIES_DIR="${QC_SUMMARIES_DIR:-$QC_DIR/summaries}"        # afreq/hardy/missing/etc.
QC_SEXCHECK_DIR="${QC_SEXCHECK_DIR:-$QC_DIR/sexcheck}"           # chrX pfiles + sexcheck table
QC_SEXCHECK_REPORTS_DIR="${QC_SEXCHECK_REPORTS_DIR:-$QC_SEXCHECK_DIR/reports}"
QC_REPORTS_DIR="${QC_REPORTS_DIR:-$QC_DIR/reports}"              # plots (e.g., PCA, vcfstats)

export OUT_DIR GTC_DIR VCF_DIR CNV_DIR QC_DIR LOG_DIR \
       PLINK_DIR TMP_DIR \
       QC_SUMMARIES_DIR QC_SEXCHECK_DIR QC_SEXCHECK_REPORTS_DIR QC_REPORTS_DIR

############################################
# Metadata (PSAM + sexmap)
############################################
META_DIR="${META_DIR:-$REPO_ROOT/metadata}"
PSAM_SEX="${PSAM_SEX:-$META_DIR/cohort.sex.psam}"
SEXMAP_PATH="${SEXMAP_PATH:-$META_DIR/sexmap.txt}"
RUN_MANIFEST="${RUN_MANIFEST:-$META_DIR/runs/${RUN}.yaml}"
export META_DIR PSAM_SEX SEXMAP_PATH RUN_MANIFEST

############################################
# Expression / Microarray integration
############################################
# input files (microarray)
: "${EXP_INDIR:=$INPUT_DIR/expression_microarray}"
: "${EXP_TGA_CSV:=$EXP_INDIR/TGA-cohort-BackUp.csv}"
: "${EXP_TRANSPOSE:=$EXP_INDIR/transpose_numbers.csv}"

# output area for expression exploration + phenotype derivation
: "${EXPR_OUT_DIR:=$OUT_DIR/expr/explore}"
: "${EXPR_LISTS_DIR:=$EXPR_OUT_DIR/lists}"
: "${EXPR_META_OVERLAP:=$EXPR_OUT_DIR/meta_overlap.csv}"
: "${EXPR_KEEP_LIST:=$EXPR_LISTS_DIR/samples_overlap.txt}"
: "${EXPR_ONLY_SNP:=$EXPR_LISTS_DIR/samples_only_in_snp.txt}"
: "${EXPR_ONLY_EXPR:=$EXPR_LISTS_DIR/samples_only_in_expression.txt}"

# optional phenotype / PSAM outputs
: "${PHENO_DIR:=$META_DIR/pheno}"
: "${PHENO_PSAM:=$PHENO_DIR/cohort.pheno.psam}"
: "${PHENO_MAP:=$PHENO_DIR/pheno_map.tsv}"

export EXP_INDIR EXP_TGA_CSV EXP_TRANSPOSE \
       EXPR_OUT_DIR EXPR_LISTS_DIR EXPR_META_OVERLAP \
       EXPR_KEEP_LIST EXPR_ONLY_SNP EXPR_ONLY_EXPR \
       PHENO_DIR PHENO_PSAM PHENO_MAP

############################################
# Threads
############################################
THREADS="${THREADS:-16}"
export THREADS

############################################
# Key VCF filenames produced by upstream steps
############################################
VCF_RAW="$VCF_DIR/cohort.gtc.$REF_BUILD.raw.vcf.gz"
VCF_NORM="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"
VCF_QC="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.vcf.gz"
VCF_SNP="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.snps.vcf.gz"
TAGS_BCF="$VCF_DIR/tmp.tags.bcf"
export VCF_RAW VCF_NORM VCF_QC VCF_SNP TAGS_BCF

############################################
# Helpers
############################################
check_tools_exist() {
  local ok=1
  command -v bcftools >/dev/null 2>&1 || { echo "[ERR] bcftools not found"; ok=0; }
  command -v plink2   >/dev/null 2>&1 || { echo "[ERR] plink2 not found";   ok=0; }
  return $ok
}

activate_env() {
  # Optional: enforce a specific conda env
  # if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then echo "[WARN] not in 'array-pipeline' env"; fi
  command -v bcftools >/dev/null 2>&1 || return 1

  # discover gtc2vcf plugin if needed
  if [ -n "${BCFTOOLS_PLUGINS:-}" ] && [ -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then return 0; fi

  local candidates=()
  if [ -n "${CONDA_PREFIX:-}" ]; then
    candidates+=("$CONDA_PREFIX/libexec/bcftools" "$CONDA_PREFIX/lib/bcftools/plugins" "$CONDA_PREFIX/share/bcftools/plugins")
  fi
  local bcftools_bin bc_root
  bcftools_bin="$(command -v bcftools || true)"
  if [ -n "$bcftools_bin" ]; then
    bc_root="$(cd -- "$(dirname -- "$bcftools_bin")/.." && pwd -P 2>/dev/null || true)"
    [ -n "$bc_root" ] && candidates+=("$bc_root/libexec/bcftools" "$bc_root/lib/bcftools/plugins" "$bc_root/share/bcftools/plugins")
  fi
  candidates+=("/usr/local/libexec/bcftools" "/usr/local/lib/bcftools/plugins" "/usr/local/share/bcftools/plugins" \
               "/usr/libexec/bcftools" "/usr/lib/bcftools/plugins" "/usr/share/bcftools/plugins")

  local p
  for p in "${candidates[@]}"; do
    if [ -f "$p/gtc2vcf.so" ]; then export BCFTOOLS_PLUGINS="$p"; return 0; fi
  done

  bcftools +gtc2vcf -h >/dev/null 2>&1 && return 0
  echo "[ERROR] bcftools plugin 'gtc2vcf' not found." >&2
  return 1
}

ensure_dirs() {
  mkdir -p \
    "$GTC_DIR" "$VCF_DIR" "$CNV_DIR" "$QC_DIR" "$LOG_DIR" \
    "$PLINK_DIR" "$TMP_DIR" \
    "$QC_SUMMARIES_DIR" "$QC_SEXCHECK_DIR" "$QC_SEXCHECK_REPORTS_DIR" "$QC_REPORTS_DIR" \
    "$VCF_DIR/tmp" "$QC_DIR/tmp" "$META_DIR" \
    "$EXPR_OUT_DIR" "$EXPR_LISTS_DIR" "$PHENO_DIR"
}

############################################
# IMPORTANT: How to use this script
############################################
# This file defines environment variables used by all downstream scripts.
# It MUST be sourced, not executed.
#
# Correct:
#   source scripts/00_config.sh
#   . scripts/00_config.sh
#
# Incorrect (variables will NOT persist):
#   bash scripts/00_config.sh
#   ./scripts/00_config.sh
#
# After sourcing, variables like REPO_ROOT, RUN, OUT_DIR, REF_BUILD, etc.
# will be available in your current shell and to any scripts you run next.
############################################

# End of config script
