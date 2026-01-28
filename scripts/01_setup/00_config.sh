#!/usr/bin/env bash
# scripts/00_config.sh
# Stable, shell-friendly config for illumina-genotyping-pipeline.
#
# IMPORTANT:
#   - This file MUST be sourced (not executed) to export variables into your shell:
#       source scripts/00_config.sh
#       . scripts/00_config.sh
#   - Works when sourced from bash or zsh (your terminal), and when sourced by bash scripts.

############################################
# Resolve repo root (works for bash + zsh)
############################################
if [ -z "${REPO_ROOT:-}" ]; then
  if [ -n "${BASH_VERSION:-}" ] && [ -n "${BASH_SOURCE[0]:-}" ]; then
    _SCRIPT="${BASH_SOURCE[0]}"
  elif [ -n "${ZSH_VERSION:-}" ]; then
    # When sourced in zsh, (%):-%x is the current file
    _SCRIPT="${(%):-%x}"
  else
    _SCRIPT="$0"
  fi

  _SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
  REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
  unset _SCRIPT _SCRIPT_DIR
fi
export REPO_ROOT

# Optional hard-assertion of expected repo root (customize as needed)
_expected_root="/home/sima/git_projects/illumina-genotyping-pipeline"
if [ "$REPO_ROOT" != "$_expected_root" ]; then
  echo "[WARN] REPO_ROOT resolved as: $REPO_ROOT" >&2
  echo "       Expected:         $_expected_root" >&2
  echo "       (If you moved the repo, ignore this. Otherwise check path resolution.)" >&2
fi
unset _expected_root

############################################
# Run label and reference build
############################################
RUN="${RUN:-genotype_run1}"         # change to genotype_run2 etc. when you start a new run
REF_BUILD="${REF_BUILD:-GRCh37}"    # GRCh37 or GRCh38
export RUN REF_BUILD

# Build-dependent constants for PLINK2
case "$REF_BUILD" in
  GRCh37) CHRSET_AUTOSOMES=37; SPLIT_PAR=b37 ;;
  GRCh38) CHRSET_AUTOSOMES=38; SPLIT_PAR=b38 ;;
  *)
    echo "[ERROR] Unknown REF_BUILD: $REF_BUILD (expected GRCh37 or GRCh38)" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac
export CHRSET_AUTOSOMES SPLIT_PAR

############################################
# Inputs (relative to repo)
############################################
INPUT_DIR="$REPO_ROOT/input_data"

IDAT_DIR="$INPUT_DIR/idat"
MANIFEST_DIR="$INPUT_DIR/manifest"
SAMPLE_SHEET_DIR="$INPUT_DIR/sample_sheet"
CLUSTER_DIR="$INPUT_DIR/cluster"

# defaults (override by exporting these before sourcing this script)
: "${BPM_MANIFEST:=$MANIFEST_DIR/GSAMD-24v3-0-EA_20034606_A1.bpm}"
: "${CSV_MANIFEST:=$MANIFEST_DIR/GSAMD-24v3-0-EA_20034606_A1.csv}"
: "${EGT_CLUSTER:=$CLUSTER_DIR/GSA-24v3-0_A1_ClusterFile_custom.egt}"

# NOTE: file name contains a space; ALWAYS quote "$SAMPLE_SHEET" when used.
: "${SAMPLE_SHEET:=$SAMPLE_SHEET_DIR/All_samples_Examine_SNPs_GWAS studies_GJ-P01_infiniumSampleSheet.csv}"
: "${SHEET_PATH:=$SAMPLE_SHEET}"

export INPUT_DIR IDAT_DIR MANIFEST_DIR SAMPLE_SHEET_DIR CLUSTER_DIR \
       BPM_MANIFEST CSV_MANIFEST EGT_CLUSTER SAMPLE_SHEET SHEET_PATH

############################################
# Reference genome
############################################
REF_DIR="$REPO_ROOT/reference/$REF_BUILD"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REF_DIR/reference.fa}"
export REF_DIR REFERENCE_FASTA

############################################
# Outputs (matches your repo layout)
############################################
OUT_DIR="${OUT_DIR:-$REPO_ROOT/output/$RUN}"
GTC_DIR="${GTC_DIR:-$OUT_DIR/gtc}"
VCF_DIR="${VCF_DIR:-$OUT_DIR/vcf}"
CNV_DIR="${CNV_DIR:-$OUT_DIR/cnv}"
QC_DIR="${QC_DIR:-$OUT_DIR/qc}"
LOG_DIR="${LOG_DIR:-$OUT_DIR/logs}"
TMP_RUN_DIR="${TMP_RUN_DIR:-$OUT_DIR/tmp}"

# PLINK working area and tmp
PLINK_DIR="${PLINK_DIR:-$QC_DIR/plink}"
PLINK_TMP_DIR="${PLINK_TMP_DIR:-$PLINK_DIR/plink_tmp}"

# QC sub-areas
QC_SUMMARIES_DIR="${QC_SUMMARIES_DIR:-$QC_DIR/summaries}"          # afreq/hardy/missing/etc.
QC_SEXCHECK_DIR="${QC_SEXCHECK_DIR:-$QC_DIR/sexcheck}"             # chrX pfiles + sexcheck table
QC_SEXCHECK_REPORTS_DIR="${QC_SEXCHECK_REPORTS_DIR:-$QC_SEXCHECK_DIR/reports}"
QC_REPORTS_DIR="${QC_REPORTS_DIR:-$QC_DIR/reports}"                # plots (e.g., PCA, vcfstats)

export OUT_DIR GTC_DIR VCF_DIR CNV_DIR QC_DIR LOG_DIR TMP_RUN_DIR \
       PLINK_DIR PLINK_TMP_DIR \
       QC_SUMMARIES_DIR QC_SEXCHECK_DIR QC_SEXCHECK_REPORTS_DIR QC_REPORTS_DIR

############################################
# Metadata (tracked)
############################################
META_DIR="${META_DIR:-$REPO_ROOT/metadata}"
PSAM_SEX="${PSAM_SEX:-$META_DIR/cohort.sex.psam}"
SEXMAP_PATH="${SEXMAP_PATH:-$META_DIR/sexmap.txt}"
SEX_OVERRIDES="${SEX_OVERRIDES:-$META_DIR/overrides/sex_overrides.txt}"
RUN_MANIFEST="${RUN_MANIFEST:-$META_DIR/runs/${RUN}.yaml}"

export META_DIR PSAM_SEX SEXMAP_PATH SEX_OVERRIDES RUN_MANIFEST

############################################
# Expression / Microarray integration
############################################
: "${EXP_INDIR:=$INPUT_DIR/expression_microarray}"
: "${EXPR_META_W:=$EXP_INDIR/expression_metadata_wide.csv}"
: "${EXPR_MATRIX_T:=$EXP_INDIR/expression_matrix_transposed.csv}"

# Base output dir for expression prep outputs (scripts 03/04)
: "${EXPR_OUT_DIR:=$OUT_DIR/expr/prep}"

# Structured subdirs (keep outputs tidy)
: "${EXPR_DIR_MATRICES:=$EXPR_OUT_DIR/matrices}"
: "${EXPR_DIR_METADATA:=$EXPR_OUT_DIR/metadata}"
: "${EXPR_DIR_CLINICAL:=$EXPR_OUT_DIR/clinical}"
: "${EXPR_DIR_LISTS:=$EXPR_OUT_DIR/lists}"
: "${EXPR_DIR_REPORTS:=$EXPR_OUT_DIR/reports}"

# Backward-compatible alias (optional): older code may still use EXPR_LISTS_DIR
: "${EXPR_LISTS_DIR:=$EXPR_DIR_LISTS}"

# Common list outputs
: "${EXPR_KEEP_LIST:=$EXPR_DIR_LISTS/samples_overlap.txt}"
: "${EXPR_ONLY_SNP:=$EXPR_DIR_LISTS/samples_only_in_snp.txt}"
: "${EXPR_ONLY_EXPR:=$EXPR_DIR_LISTS/samples_only_in_expression.txt}"

# (Optional legacy var; keep only if something still writes/reads it)
: "${EXPR_META_OVERLAP:=$EXPR_DIR_METADATA/meta_overlap.csv}"

: "${PHENO_DIR:=$META_DIR/pheno}"
: "${PHENO_PSAM:=$PHENO_DIR/cohort.pheno.psam}"
: "${PHENO_MAP:=$PHENO_DIR/pheno_map.tsv}"

export EXP_INDIR EXPR_META_W EXPR_MATRIX_T \
       EXPR_OUT_DIR \
       EXPR_DIR_MATRICES EXPR_DIR_METADATA EXPR_DIR_CLINICAL EXPR_DIR_LISTS EXPR_DIR_REPORTS \
       EXPR_LISTS_DIR EXPR_KEEP_LIST EXPR_ONLY_SNP EXPR_ONLY_EXPR EXPR_META_OVERLAP \
       PHENO_DIR PHENO_PSAM PHENO_MAP


############################################
# Threads
############################################
THREADS="${THREADS:-16}"
export THREADS

############################################
# Key VCF filenames (produced by upstream steps)
############################################
VCF_RAW="$VCF_DIR/cohort.gtc.$REF_BUILD.raw.vcf.gz"
VCF_NORM="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"
VCF_QC="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.vcf.gz"
VCF_SNP="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.snps.vcf.gz"
TAGS_BCF="$VCF_DIR/tmp/tags.bcf"
export VCF_RAW VCF_NORM VCF_QC VCF_SNP TAGS_BCF

############################################
# Helpers
############################################
check_tools_exist() {
  local missing=0
  command -v bcftools >/dev/null 2>&1 || { echo "[ERR] bcftools not found"; missing=1; }
  command -v plink2   >/dev/null 2>&1 || { echo "[ERR] plink2 not found";   missing=1; }
  return "$missing"
}

activate_env() {
  # Optional: enforce a specific conda env
  # if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then echo "[WARN] not in 'array-pipeline' env"; fi

  command -v bcftools >/dev/null 2>&1 || return 1

  # If already set and valid, keep it
  if [ -n "${BCFTOOLS_PLUGINS:-}" ] && [ -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then
    return 0
  fi

  # Candidate plugin directories (conda + bcftools-relative + system)
  local candidates=()
  if [ -n "${CONDA_PREFIX:-}" ]; then
    candidates+=("$CONDA_PREFIX/libexec/bcftools" "$CONDA_PREFIX/lib/bcftools/plugins" "$CONDA_PREFIX/share/bcftools/plugins")
  fi

  local bcftools_bin bc_root
  bcftools_bin="$(command -v bcftools 2>/dev/null || true)"
  if [ -n "$bcftools_bin" ]; then
    bc_root="$(cd -- "$(dirname -- "$bcftools_bin")/.." && pwd -P 2>/dev/null || true)"
    [ -n "$bc_root" ] && candidates+=("$bc_root/libexec/bcftools" "$bc_root/lib/bcftools/plugins" "$bc_root/share/bcftools/plugins")
  fi

  candidates+=(
    "/usr/local/libexec/bcftools" "/usr/local/lib/bcftools/plugins" "/usr/local/share/bcftools/plugins"
    "/usr/libexec/bcftools"       "/usr/lib/bcftools/plugins"       "/usr/share/bcftools/plugins"
  )

  local p
  for p in "${candidates[@]}"; do
    if [ -f "$p/gtc2vcf.so" ]; then
      export BCFTOOLS_PLUGINS="$p"
      return 0
    fi
  done

  # As a last check, see if bcftools can load it anyway
  bcftools +gtc2vcf -h >/dev/null 2>&1 && return 0

  echo "[ERROR] bcftools plugin 'gtc2vcf' not found (gtc2vcf.so)." >&2
  return 1
}

ensure_dirs() {
  mkdir -p \
    "$GTC_DIR" "$VCF_DIR" "$CNV_DIR" "$QC_DIR" "$LOG_DIR" "$TMP_RUN_DIR" \
    "$PLINK_DIR" "$PLINK_TMP_DIR" \
    "$QC_SUMMARIES_DIR" "$QC_SEXCHECK_DIR" "$QC_SEXCHECK_REPORTS_DIR" "$QC_REPORTS_DIR" \
    "$VCF_DIR/tmp" "$QC_DIR/tmp" \
    "$META_DIR" "$META_DIR/overrides" "$META_DIR/runs" \
    "$EXPR_OUT_DIR" \
    "$EXPR_DIR_MATRICES" "$EXPR_DIR_METADATA" "$EXPR_DIR_CLINICAL" "$EXPR_DIR_LISTS" "$EXPR_DIR_REPORTS" \
    "$PHENO_DIR"
}
############################################
# End of config
############################################
