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
    REFERENCE_FASTA="$REF_DIR/reference.fa"
    ;;
  GRCh38)
    REFERENCE_FASTA="$REF_DIR/reference.fa"
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
# ---------- Environment & plugin setup ----------
activate_env() {
  # Optional: enforce a specific conda env name (uncomment to warn)
  # if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then
  #   echo "[WARN] Not in 'array-pipeline' env; current: ${CONDA_DEFAULT_ENV:-<none>} (this is OK if tools are available)" >&2
  # fi

  # 0) bcftools present?
  if ! command -v bcftools >/dev/null 2>&1; then
    echo "[ERROR] 'bcftools' not found on PATH. Activate your env (e.g., 'conda activate array-pipeline') or install bcftools." >&2
    return 1
  fi

  # If BCFTOOLS_PLUGINS already points to a valid gtc2vcf.so, we're done
  if [ -n "${BCFTOOLS_PLUGINS:-}" ] && [ -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then
    return 0
  fi

  # 1) Build candidate plugin paths
  local candidates=()
  # Conda-style locations
  if [ -n "${CONDA_PREFIX:-}" ]; then
    candidates+=("$CONDA_PREFIX/libexec/bcftools")
    candidates+=("$CONDA_PREFIX/lib/bcftools/plugins")
    candidates+=("$CONDA_PREFIX/share/bcftools/plugins")
  fi
  # Paths relative to the bcftools binary (covers Homebrew/system installs)
  local bcftools_bin bcftools_root
  bcftools_bin="$(command -v bcftools)"
  if [ -n "$bcftools_bin" ]; then
    # Typically .../bin/bcftools → root = one dir up
    bcftools_root="$(cd -- "$(dirname -- "$bcftools_bin")/.." && pwd -P 2>/dev/null || true)"
    if [ -n "$bcftools_root" ]; then
      candidates+=("$bcftools_root/libexec/bcftools")
      candidates+=("$bcftools_root/lib/bcftools/plugins")
      candidates+=("$bcftools_root/share/bcftools/plugins")
    fi
  fi
  # Common system-wide fallbacks
  candidates+=("/usr/local/libexec/bcftools" "/usr/local/lib/bcftools/plugins" "/usr/local/share/bcftools/plugins")
  candidates+=("/usr/libexec/bcftools" "/usr/lib/bcftools/plugins" "/usr/share/bcftools/plugins")

  # 2) Pick the first candidate that actually contains gtc2vcf.so
  local p found=""
  for p in "${candidates[@]}"; do
    if [ -f "$p/gtc2vcf.so" ]; then
      found="$p"
      break
    fi
  done

  if [ -n "$found" ]; then
    export BCFTOOLS_PLUGINS="$found"
    return 0
  fi

  # 3) As a last resort, see if bcftools can load by name (+gtc2vcf). If yes, we don't need BCFTOOLS_PLUGINS.
  if bcftools +gtc2vcf -h >/dev/null 2>&1; then
    # Works without setting BCFTOOLS_PLUGINS (some builds bake it in)
    return 0
  fi

  # 4) If we got here, plugin is missing or not discoverable
  echo "[ERROR] bcftools plugin 'gtc2vcf' not found." >&2
  echo "        Looked in:" >&2
  printf '          - %s\n' "${candidates[@]}" >&2
  echo "        Fixes:" >&2
  echo "          • Ensure you installed a bcftools build that ships plugins (Conda/bioconda recommended)." >&2
  echo "          • Verify the plugin file exists:  find \"\${CONDA_PREFIX:-/}\" -name gtc2vcf.so 2>/dev/null" >&2
  echo "          • If you built gtc2vcf yourself, set:  export BCFTOOLS_PLUGINS=/path/to/plugins" >&2
  return 1
}

# ---------- Directories ----------
ensure_dirs() {
  mkdir -p "$GTC_DIR" "$VCF_DIR" "$CNV_DIR" "$QC_DIR" "$LOG_DIR" "$TMP_DIR"
  # bcftools sort likes a temp dir; ensure one under VCF_DIR as well
  mkdir -p "$VCF_DIR/tmp" "$QC_DIR/tmp" 2>/dev/null || true
}

# ---------- Input preflight ----------
check_inputs_exist() {
  local ok=1

  # Required files
  local required=("$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA" "$SAMPLE_SHEET")
  local f
  for f in "${required[@]}"; do
    if [ ! -s "$f" ]; then
      echo "[ERR] Missing or empty file: $f" >&2
      ok=0
    fi
  done

  # Optional but recommended: CSV_MANIFEST
  if [ -n "${CSV_MANIFEST:-}" ] && [ ! -s "$CSV_MANIFEST" ]; then
    echo "[WARN] CSV_MANIFEST is set but not found or empty: $CSV_MANIFEST (continuing without it)" >&2
  fi

  if [ $ok -eq 0 ]; then
    echo "[HINT] Expected layout:" >&2
    echo "       • Manifests (BPM/CSV):    input_data/manifest/" >&2
    echo "       • Cluster (EGT):          input_data/cluster/" >&2
    echo "       • Reference FASTA:        reference/$REF_BUILD/" >&2
    echo "       • Sample sheet (CSV):     input_data/sample_sheet/" >&2
    return 1
  fi
  return 0
}

############################################
# Tell the pipeline exactly where dragena.exe lives (Windows install)
# DRAGENA_BIN_OVERRIDE="/mnt/c/Program Files/Illumina/DRAGENArrayLocal/dragena-win-x64-DAv1.3.0-rc3-sha.f3fec02ebf2c43d3f3d6327cbe3b410edbc167b4/dragena/dragena.exe"
############################################

# This file should be sourced by numbered scripts; do not 'set -e' here.