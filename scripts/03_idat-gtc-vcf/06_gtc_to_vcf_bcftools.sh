#!/usr/bin/env bash
# 06_gtc_to_vcf_cohort.sh — Build cohort-only GTC set from keep list, then convert to cohort-only VCF/BCF.
# Usage:
#   source scripts/00_config.sh
#   bash scripts/06_gtc_to_vcf_cohort.sh
# Notes:
#   - Requires: bcftools (+gtc2vcf plugin), tabix; optional samtools for FASTA index
#   - Inputs from 00_config.sh: BPM_MANIFEST, EGT_CLUSTER, REFERENCE_FASTA, CSV_MANIFEST(optional),
#                               GTC_DIR, VCF_DIR, TMP_DIR, LOG_DIR, PHENO_DIR, RUN, REF_BUILD
#   - Keep list (IIDs): $PHENO_DIR/iid_selected.keep (one IID per line, e.g., 207363850090_R12C01)

set -euo pipefail

# ---------- Resolve repo root & load config ----------
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# ---------- Prep & log ----------
activate_env || true
ensure_dirs
mkdir -p "$VCF_DIR" "$TMP_DIR" "$LOG_DIR"

LOG="$LOG_DIR/06_gtc_to_vcf_cohort.log"
exec > >(tee "$LOG") 2>&1

echo "== GTC -> VCF (BCF-first, cohort-only) =="
echo "RUN=$RUN  REF_BUILD=$REF_BUILD"
echo "GTC_DIR=$GTC_DIR"
echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"

# ---------- Tool checks ----------
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "[ERROR] tabix not found on PATH"; exit 1; }
if ! bcftools --version | grep -q 'bcftools'; then
  echo "[ERROR] bcftools not working as expected."; exit 1
fi

# ---------- Input checks ----------
miss=0
for f in "$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA"; do
  [[ -s "$f" ]] || { echo "[ERR] Missing: $f"; miss=1; }
done
[[ "$miss" -eq 0 ]] || { echo "[HINT] Fix paths in scripts/00_config.sh"; exit 1; }

# FASTA index (soft)
if [[ ! -s "${REFERENCE_FASTA}.fai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    echo "[INFO] Building FASTA index with samtools faidx"
    samtools faidx "$REFERENCE_FASTA" || true
  else
    echo "[WARN] samtools not found; ${REFERENCE_FASTA}.fai missing (bcftools may still proceed)"
  fi
fi

# CSV manifest arg (optional)
CSV_ARG=()
if [[ -n "${CSV_MANIFEST:-}" && -s "$CSV_MANIFEST" ]]; then
  CSV_ARG=(-c "$CSV_MANIFEST")
  echo "[INFO] Using CSV manifest: $CSV_MANIFEST"
else
  guess_csv="${BPM_MANIFEST%.bpm}.csv"
  if [[ -s "$guess_csv" ]]; then
    CSV_ARG=(-c "$guess_csv")
    echo "[INFO] Using CSV manifest (guessed): $guess_csv"
  else
    echo "[WARN] No CSV manifest found; proceeding without -c"
  fi
fi

# ---------- Build cohort-only GTC folder from keep list ----------
KEEP="$PHENO_DIR/iid_selected.keep"
[[ -s "$KEEP" ]] || { echo "[ERROR] Keep list missing/empty: $KEEP"; exit 1; }

SEL_DIR="$GTC_DIR/_selected_cohort"
rm -rf "$SEL_DIR"
mkdir -p "$SEL_DIR"

total_keep=$(wc -l < "$KEEP" | tr -d ' ')
echo "[INFO] iid_selected.keep: $total_keep"

missing=0
while IFS=$'\n' read -r iid; do
  [[ -n "$iid" ]] || continue
  src_exact="$GTC_DIR/${iid}.gtc"
  if [[ -s "$src_exact" ]]; then
    ln -sf "../$(basename "$src_exact")" "$SEL_DIR/${iid}.gtc"
    continue
  fi
  # fallback: fuzzy match (must be unique)
  mapfile -t cand < <(find "$GTC_DIR" -maxdepth 1 -type f -name "*${iid}*.gtc")
  if (( ${#cand[@]} == 1 )); then
    ln -sf "../$(basename "${cand[0]}")" "$SEL_DIR/${iid}.gtc"
  else
    echo "[WARN] No unique GTC match for IID: $iid"
    ((missing++))
  fi
done < "$KEEP"

sel_count=$(find -L "$SEL_DIR" -maxdepth 1 -name '*.gtc' | wc -l | tr -d ' ')
echo "[INFO] Selected GTC files detected: $sel_count"
if (( sel_count == 0 )); then
  echo "[ERROR] No selected .gtc files found in $SEL_DIR"; exit 1
fi
if (( missing > 0 )); then
  echo "[WARN] $missing IIDs had no unique GTC; proceeding with $sel_count selected."
fi

# Optional sanity: list diffs
# comm -3 <(basename -a "$SEL_DIR"/*.gtc | sed 's/\.gtc$//' | sort) <(sort "$KEEP") || true

# ---------- bcftools +gtc2vcf on cohort-only GTCs ----------
SORT_TMP="$TMP_DIR/bcfsort_${RUN}"; mkdir -p "$SORT_TMP"
gtc2vcf_plugin="+gtc2vcf"
if [[ -n "${BCFTOOLS_PLUGINS:-}" && -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]]; then
  gtc2vcf_plugin="+$BCFTOOLS_PLUGINS/gtc2vcf.so"
fi

RAW_BCF="$VCF_DIR/cohort.gtc.$REF_BUILD.raw.bcf"
NORM_BCF="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.bcf"
VCF_GZ="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"
TAGS_BCF="$VCF_DIR/cohort.gtc.$REF_BUILD.tags.bcf"

echo "[STEP] gtc2vcf -> sort (RAW BCF): $RAW_BCF"
set -x
bcftools "$gtc2vcf_plugin" \
  --no-version -Ou \
  -b "$BPM_MANIFEST" "${CSV_ARG[@]}" -e "$EGT_CLUSTER" -f "$REFERENCE_FASTA" \
  -g "$SEL_DIR" \
| bcftools sort -Ob -T "$SORT_TMP" -o "$RAW_BCF"
set +x
bcftools index -f "$RAW_BCF" || true
echo "[OK] RAW BCF created & indexed: $RAW_BCF"

echo "[STEP] Normalize (split, left-align, drop REF mismatches) -> NORM BCF"
set -x
bcftools norm --no-version \
  -m -any -c x -f "$REFERENCE_FASTA" \
  -Ob -o "$NORM_BCF" --write-index "$RAW_BCF"
set +x
echo "[OK] NORM BCF created & indexed: $NORM_BCF(.csi)"

echo "[STEP] Convert NORM BCF -> gzipped VCF & index"
set -x
bcftools view -Oz -o "$VCF_GZ" "$NORM_BCF"
tabix -f -p vcf "$VCF_GZ"
set +x
echo "[OK] Wrote: $VCF_GZ (+ .tbi)"

echo "[STEP] Fill tags (AN, AC, AF) -> $TAGS_BCF"
set -x
bcftools view -Ou "$NORM_BCF" \
| bcftools +fill-tags -Ob -o "$TAGS_BCF" -- -t AN,AC,AF
bcftools index -f "$TAGS_BCF" || true
set +x
echo "[OK] Tags BCF created & indexed: $TAGS_BCF(.csi)"

# ---------- Summaries ----------
samples_raw=$(bcftools query -l "$RAW_BCF" | wc -l | tr -d ' ')
samples_norm=$(bcftools query -l "$NORM_BCF" | wc -l | tr -d ' ')
records_norm=$(bcftools index -n "$NORM_BCF")
echo "[SUMMARY] Selected GTCs      : $sel_count (keep list: $total_keep)"
echo "[SUMMARY] RAW_BCF samples    : $samples_raw"
echo "[SUMMARY] NORM_BCF samples   : $samples_norm"
echo "[SUMMARY] NORM_BCF records   : $records_norm"
bcftools view -h "$VCF_GZ" | grep -m1 '^#CHROM' || true

echo "== Done: cohort-only GTC -> VCF =="
