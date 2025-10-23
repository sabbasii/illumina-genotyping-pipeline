#!/usr/bin/env bash
# Soft prereq checker: never exit on missing tools; only WARN.
set -euo pipefail

# Locate repo root and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

echo "== Input sanity-check =="
echo "RUN=$RUN   REF_BUILD=$REF_BUILD"
echo "REPO_ROOT=$REPO_ROOT"
echo

# --- Tools (soft check) ---
missing_tools=()
command -v bcftools >/dev/null 2>&1 || missing_tools+=("bcftools")
command -v plink2   >/dev/null 2>&1 || missing_tools+=("plink2")
command -v samtools >/dev/null 2>&1 || missing_tools+=("samtools")

if ((${#missing_tools[@]})); then
  echo "[WARN] Missing tools: ${missing_tools[*]}"
  echo "       (Tip) conda activate array-pipeline"
else
  echo "[OK] Tools present: bcftools, plink2, samtools"
fi
echo

# Helper: print OK/MISSING
report() {
  local label="$1" path="$2"
  local status="MISSING"
  [[ -s "$path" ]] && status="OK"
  printf "%-18s : %-7s  %s\n" "$label" "$status" "$path"
}

echo "== Required inputs =="
report "BPM_MANIFEST"     "$BPM_MANIFEST"
if [[ -n "${CSV_MANIFEST:-}" ]]; then
  report "CSV_MANIFEST"   "$CSV_MANIFEST"
else
  printf "%-18s : %-7s  %s\n" "CSV_MANIFEST" "SKIP" "(not set)"
fi
report "EGT_CLUSTER"      "$EGT_CLUSTER"
report "REFERENCE_FASTA"  "$REFERENCE_FASTA"
report "SAMPLE_SHEET"     "$SAMPLE_SHEET"
echo

# Warn (don???t exit) on missing required inputs
missing=0
for f in "$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA" "$SAMPLE_SHEET"; do
  [[ -s "$f" ]] || { echo "[WARN] Missing: $f"; missing=1; }
done
((missing==0)) && echo "[OK] All required inputs present." || echo "[WARN] Some required inputs are missing (see above)."
echo

# FASTA index (only if samtools is available)
if [[ -f "$REFERENCE_FASTA" ]]; then
  if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      echo "[INFO] Indexing FASTA with samtools faidx..."
      samtools faidx "$REFERENCE_FASTA"
      echo "[OK] Created: ${REFERENCE_FASTA}.fai"
    else
      echo "[WARN] samtools not found; cannot index FASTA automatically."
    fi
  else
    echo "[OK] FASTA index present: ${REFERENCE_FASTA}.fai"
  fi
else
  echo "[WARN] REFERENCE_FASTA not found; skipping index check."
fi
echo

# IDAT scan (case-insensitive)
echo "[INFO] Scanning for IDAT files under: $IDAT_DIR"
idat_n="$(find "$IDAT_DIR" -type f \( -iname '*.idat' \) 2>/dev/null | wc -l | tr -d ' ')"
echo "IDAT files found: $idat_n"
if [[ "$idat_n" -eq 0 ]]; then
  echo "[WARN] No IDAT files found. Place .idat files under: $IDAT_DIR"
else
  echo "[INFO] First few IDATs:"
  # Avoid SIGPIPE with head: collect then slice
  mapfile -t _first_idats < <(find "$IDAT_DIR" -type f \( -iname '*.idat' \) 2>/dev/null)
  for ((i=0; i<${#_first_idats[@]} && i<6; i++)); do
    printf '%s\n' "${_first_idats[$i]}"
  done
fi

# Ensure all output dirs exist per new layout
ensure_dirs
echo "[OK] Output dirs ensured:"
printf "  %s\n" \
  "$GTC_DIR" "$VCF_DIR" "$QC_DIR" "$LOG_DIR" \
  "$PLINK_DIR" "$TMP_DIR" \
  "$QC_SUMMARIES_DIR" "$QC_SEXCHECK_DIR" "$QC_SEXCHECK_REPORTS_DIR" "$QC_REPORTS_DIR"
echo

echo "[INFO] Run manifest target: $RUN_MANIFEST"
echo "== Done =="
