#!/usr/bin/env bash
set -euo pipefail

# Locate repo root (portable) and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"

# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

echo "== Input sanity-check =="
echo "RUN=$RUN   REF_BUILD=$REF_BUILD"
echo

# Print OK/MISSING for key files
report() {
  local label="$1" path="$2"
  local status="MISSING"
  [[ -s "$path" ]] && status="OK"
  printf "%-16s : %-7s  %s\n" "$label" "$status" "$path"
}

report "BPM_MANIFEST"   "$BPM_MANIFEST"
if [[ -n "${CSV_MANIFEST:-}" ]]; then
  report "CSV_MANIFEST" "$CSV_MANIFEST"
else
  printf "%-16s : %-7s  %s\n" "CSV_MANIFEST" "SKIP" "(not set)"
fi
report "EGT_CLUSTER"    "$EGT_CLUSTER"
report "REFERENCE_FASTA" "$REFERENCE_FASTA"
report "SAMPLE_SHEET"   "$SAMPLE_SHEET"
echo

# Fail early if required files are missing
missing=0
for f in "$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA" "$SAMPLE_SHEET"; do
  [[ -s "$f" ]] || missing=1
done
if [[ $missing -ne 0 ]]; then
  echo "[ERROR] One or more required inputs are missing. Fix paths in scripts/00_config.sh."
  exit 1
fi

# Ensure FASTA is indexed
if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
  echo "[INFO] Indexing FASTA with samtools faidx..."
  samtools faidx "$REFERENCE_FASTA"
else
  echo "[OK] FASTA index present: ${REFERENCE_FASTA}.fai"
fi
echo

# IDAT presence check
echo "[INFO] Scanning for IDAT files under: $IDAT_DIR"
idat_n="$(find "$IDAT_DIR" -type f -name '*.idat' | wc -l | tr -d ' ')"
echo "IDAT files found: $idat_n"
if [[ "$idat_n" -eq 0 ]]; then
  echo "[WARN] No IDAT files found. Place .idat files under $IDAT_DIR before running genotyping."
else
  echo "[INFO] First few IDATs:"
  find "$IDAT_DIR" -type f -name '*.idat' | head -n 6
fi
echo

# Ensure output dirs exist
ensure_dirs
echo "[OK] Output dirs ensured:"
printf "  %s\n" "$GTC_DIR" "$VCF_DIR" "$QC_DIR" "$LOG_DIR" "$TMP_DIR"
echo
echo "== Done =="
