#!/usr/bin/env bash
# NOTE: This is a template. Copy it to scripts/00_config.sh and edit the few lines marked EDIT ME.

############################################
# Resolve repo root (portable)
############################################
# Works when run via 'bash scripts/<...>.sh' from repo root.
# If you cd elsewhere, REPO_ROOT fallback still resolves relative to this file.
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"

############################################
# Run label and reference build
############################################
RUN="genotype_run1"        # EDIT ME if you start a new run (e.g., genotype_run2)
REF_BUILD="GRCh37"         # EDIT ME to 'GRCh38' if you switch builds later

############################################
# Input locations (relative paths)
############################################
INPUT_DIR="$REPO_ROOT/input_data"
IDAT_DIR="$INPUT_DIR/idat"
MANIFEST_DIR="$INPUT_DIR/manifest"
SAMPLE_SHEET_DIR="$INPUT_DIR/sample_sheet"
CLUSTER_DIR="$INPUT_DIR/cluster"

# EDIT ME: set your actual filenames
BPM_MANIFEST="$MANIFEST_DIR/GSAMD-24v3-0-EA_20034606_A1.bpm"     # example
CSV_MANIFEST="$MANIFEST_DIR/GSA-24v3-0_A1.csv"                   # optional but recommended
EGT_CLUSTER="$CLUSTER_DIR/GSA-24v3-0_A1_ClusterFile_custom.egt"  # example
SAMPLE_SHEET="$SAMPLE_SHEET_DIR/YourSampleSheet.csv"             # e.g., All_samples_...csv

############################################
# Reference genome
############################################
REF_DIR="$REPO_ROOT/reference/$REF_BUILD"

# Set the FASTA file per build (adjust names to what you actually place in reference/)
case "$REF_BUILD" in
  GRCh37)
    REFERENCE_FASTA="$REF_DIR/human_g1k_v37.fasta"  # e.g., from 1000G
    ;;
  GRCh38)
    # Example: NCBI no-alt analysis set with UCSC-style IDs; ensure it matches your manifest naming
    REFERENCE_FASTA="$REF_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    ;;
  *)
    echo "[ERROR] Unknown REF_BUILD: $REF_BUILD" >&2; exit 1;;
esac

############################################
# Output locations (relative paths)
############################################
OUT_DIR="$REPO_ROOT/output/$RUN"
GTC_DIR="$OUT_DIR/gtc"
VCF_DIR="$OUT_DIR/vcf"
CNV_DIR="$OUT_DIR/cnv"
QC_DIR="$OUT_DIR/qc"
LOG_DIR="$OUT_DIR/logs"
TMP_DIR="$OUT_DIR/tmp"

############################################
# Threading
############################################
THREADS="${THREADS:-16}"   # override by exporting THREADS before running scripts

############################################
# Key filenames produced by steps
############################################
# Raw multi-sample VCF from gtc2vcf
VCF_RAW="$VCF_DIR/cohort.gtc.$REF_BUILD.raw.vcf.gz"

# Sorted/left-aligned/split-normalized VCF
VCF_NORM="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"

# QC-filtered all-site VCF
VCF_QC="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.vcf.gz"

# SNP-only, biallelic, polymorphic set for PLINK import
VCF_SNP="$VCF_DIR/cohort.gtc.$REF_BUILD.qc.snps.vcf.gz"

# Tag-enriched BCF (AC/AN/AF) used during filtering
TAGS_BCF="$VCF_DIR/tmp.tags.bcf"

############################################
# Tooling helpers (no hard exits here)
############################################
activate_env() {
  # If you use conda, ensure the env is active; otherwise no-op
  # Optionally enforce a specific env name by uncommenting the next lines:
  # if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then
  #   echo "[WARN] Not in 'array-pipeline' env; current: ${CONDA_DEFAULT_ENV:-<none>}"; fi

  # Try to set bcftools plugin path if not already set
  if [ -z "${BCFTOOLS_PLUGINS:-}" ] && [ -n "${CONDA_PREFIX:-}" ]; then
    for p in "$CONDA_PREFIX/libexec/bcftools" "$CONDA_PREFIX/lib/bcftools/plugins" "$CONDA_PREFIX/share/bcftools/plugins"; do
      [ -d "$p" ] && export BCFTOOLS_PLUGINS="$p" && break
    done
  fi
}

ensure_dirs() {
  mkdir -p "$GTC_DIR" "$VCF_DIR" "$CNV_DIR" "$QC_DIR" "$LOG_DIR" "$TMP_DIR"
}

check_inputs_exist() {
  local ok=1
  for f in "$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA" "$SAMPLE_SHEET"; do
    if [ ! -s "$f" ]; then echo "[ERR] Missing file: $f" >&2; ok=0; fi
  done
  if [ $ok -eq 0 ]; then
    echo "[HINT] Place manifests under input_data/manifest/, cluster under input_data/cluster/, reference FASTA under reference/$REF_BUILD/" >&2
    return 1
  fi
  return 0
}

# This file should be sourced by numbered scripts; do not 'set -e' here.
