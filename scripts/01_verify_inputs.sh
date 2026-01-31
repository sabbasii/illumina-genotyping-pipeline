#!/usr/bin/env bash
# scripts/01_verify_inputs.sh
# Soft prereq checker: never exit on missing inputs; only WARN.
set -euo pipefail

############################################
# Resolve repo root from this script location
############################################
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
export REPO_ROOT
unset _SCRIPT _SCRIPT_DIR

# Sanity-check repo layout BEFORE sourcing config
if [[ ! -f "$REPO_ROOT/scripts/00_config.sh" ]]; then
  echo "[ERROR] REPO_ROOT does not look like the repo root: $REPO_ROOT" >&2
  echo "        Missing: $REPO_ROOT/scripts/00_config.sh" >&2
  exit 1
fi

# Load config (single source of truth)
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

echo "== Input sanity-check =="
echo "RUN=$RUN   REF_BUILD=$REF_BUILD"
echo "REPO_ROOT=$REPO_ROOT"
echo

############################################
# Tools (soft check)
############################################
missing_tools=()
command -v bcftools >/dev/null 2>&1 || missing_tools+=("bcftools")
command -v plink2   >/dev/null 2>&1 || missing_tools+=("plink2")
command -v samtools >/dev/null 2>&1 || missing_tools+=("samtools")
command -v python3  >/dev/null 2>&1 || missing_tools+=("python3 (optional for expr prep)")

if ((${#missing_tools[@]})); then
  echo "[WARN] Missing tools: ${missing_tools[*]}"
  echo "       (Tip) conda activate array-pipeline"
else
  echo "[OK] Tools present: bcftools, plink2, samtools, python3"
fi
echo

############################################
# Required inputs (genotyping)
############################################
echo "== Required inputs (genotyping) =="

required_missing=0

check_required() {
  local label="$1" path="$2"
  local status="MISSING"
  if [[ -s "$path" ]]; then
    status="OK"
  else
    required_missing=1
  fi
  printf "%-22s : %-7s  %s\n" "$label" "$status" "$path"
}

check_optional() {
  local label="$1" path="${2:-}"
  if [[ -n "$path" ]]; then
    local status="MISSING"
    [[ -s "$path" ]] && status="OK"
    printf "%-22s : %-7s  %s\n" "$label" "$status" "$path"
  else
    printf "%-22s : %-7s  %s\n" "$label" "SKIP" "not set"
  fi
}

check_required "BPM_MANIFEST"    "$BPM_MANIFEST"
check_optional "CSV_MANIFEST"    "${CSV_MANIFEST:-}"
check_required "EGT_CLUSTER"     "$EGT_CLUSTER"
check_required "REFERENCE_FASTA" "$REFERENCE_FASTA"
check_required "SAMPLE_SHEET"    "$SAMPLE_SHEET"
echo

if (( required_missing == 0 )); then
  echo "[OK] All required genotyping inputs are present."
else
  echo "[WARN] One or more required genotyping inputs are missing (see above)."
fi
echo

############################################
# FASTA index (create .fai when needed)
############################################
if [[ -f "$REFERENCE_FASTA" ]]; then
  if [[ -f "${REFERENCE_FASTA}.fai" ]]; then
    echo "[OK] FASTA index present: ${REFERENCE_FASTA}.fai"
  elif command -v samtools >/dev/null 2>&1; then
    echo "[INFO] Creating FASTA index: samtools faidx"
    samtools faidx "$REFERENCE_FASTA"
    echo "[OK] Created: ${REFERENCE_FASTA}.fai"
  else
    echo "[WARN] samtools not available; cannot create FASTA index automatically."
  fi
else
  echo "[WARN] REFERENCE_FASTA not found; skipping FASTA index check."
fi
echo

############################################
# IDAT presence check
############################################
echo "== IDAT scan =="
echo "[INFO] Scanning: $IDAT_DIR"

idat_n="$(find "$IDAT_DIR" -type f -iname '*.idat' 2>/dev/null | wc -l | tr -d ' ' || true)"
echo "IDAT files found: $idat_n"

if [[ "${idat_n:-0}" -eq 0 ]]; then
  echo "[WARN] No IDAT files found under: $IDAT_DIR"
else
  echo "[INFO] Example IDAT files:"
  find "$IDAT_DIR" -type f -iname '*.idat' 2>/dev/null | head -n 6 || true
fi
echo

############################################
# Expression / Microarray inputs check
############################################
echo "== Expression / microarray inputs =="

report() {
  local label="$1" path="$2"
  local status="MISSING"
  [[ -s "$path" ]] && status="OK"
  printf "%-22s : %-7s  %s\n" "$label" "$status" "$path"
}

report "EXPR_META_W"   "$EXPR_META_W"
report "EXPR_MATRIX_T" "$EXPR_MATRIX_T"

if [[ -s "$EXPR_META_W" && -s "$EXPR_MATRIX_T" ]]; then
  sz_csv="$(du -h "$EXPR_META_W" 2>/dev/null | awk '{print $1}' || true)"
  sz_trn="$(du -h "$EXPR_MATRIX_T" 2>/dev/null | awk '{print $1}' || true)"
  echo "[INFO] File sizes: EXPR_META_W=$sz_csv  EXPR_MATRIX_T=$sz_trn"
else
  echo "[INFO] Expression inputs not detected; expression-related steps will be skipped."
fi
echo

############################################
# Ensure output directories exist
############################################
ensure_dirs

echo "[OK] Output directories ensured:"
dirs=(
  "$GTC_DIR" "$VCF_DIR" "$QC_DIR" "$LOG_DIR" "$TMP_RUN_DIR"
  "$PLINK_DIR" "$PLINK_TMP_DIR"
  "$QC_SUMMARIES_DIR" "$QC_SEXCHECK_DIR" "$QC_SEXCHECK_REPORTS_DIR" "$QC_REPORTS_DIR"
  "$EXPR_OUT_DIR"
  "$EXPR_DIR_MATRICES" "$EXPR_DIR_METADATA" "$EXPR_DIR_CLINICAL" "$EXPR_DIR_LISTS" "$EXPR_DIR_REPORTS"
  "$PHENO_DIR"
)

printf '  %s\n' "${dirs[@]}"
unset dirs
echo

echo "[INFO] Run manifest target: $RUN_MANIFEST"
echo "== Done =="

# Usage:
#   bash scripts/01_verify_inputs.sh
#
# Notes:
# - This script sources scripts/00_config.sh internally for its own execution.
# - The configuration is loaded only within this script’s process.
# - When the script finishes, all variables defined by 00_config.sh
#   are discarded and do not persist in the current shell.
#
# To make configuration variables available for subsequent commands
# or other scripts, source scripts/00_config.sh explicitly:
#   source scripts/00_config.sh

# End of file