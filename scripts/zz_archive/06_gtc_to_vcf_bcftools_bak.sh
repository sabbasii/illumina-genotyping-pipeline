#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root & load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

activate_env || true
ensure_dirs

LOG="$LOG_DIR/03_gtc_to_vcf_bcftools.log"
echo "== GTC -> VCF (BCF-first, bcftools +gtc2vcf with SHORT flags) ==" | tee "$LOG"
echo "RUN=$RUN  REF_BUILD=$REF_BUILD" | tee -a "$LOG"
echo "GTC_DIR=$GTC_DIR" | tee -a "$LOG"
echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}" | tee -a "$LOG"

# Tools
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH" | tee -a "$LOG"; exit 1; }

# Inputs
miss=0
for f in "$BPM_MANIFEST" "$EGT_CLUSTER" "$REFERENCE_FASTA"; do
  [[ -s "$f" ]] || { echo "[ERR] Missing: $f" | tee -a "$LOG"; miss=1; }
done
[[ "$miss" -eq 0 ]] || { echo "[HINT] Check paths in scripts/00_config.sh" | tee -a "$LOG"; exit 1; }

# FASTA index (optional but helpful)
if command -v samtools >/dev/null 2>&1 && [[ ! -s "${REFERENCE_FASTA}.fai" ]]; then
  echo "[INFO] Building FASTA index with samtools faidx" | tee -a "$LOG"
  samtools faidx "$REFERENCE_FASTA" 2>>"$LOG" || true
fi

# GTC presence
gtc_n=$(find "$GTC_DIR" -maxdepth 1 -type f -name '*.gtc' | wc -l | tr -d ' ')
echo "[INFO] GTC files detected: $gtc_n" | tee -a "$LOG"
[[ "$gtc_n" -gt 0 ]] || { echo "[ERROR] No .gtc files in $GTC_DIR" | tee -a "$LOG"; exit 1; }

# CSV (prefer provided; else guess next to BPM)
CSV_ARG=()
if [[ -n "${CSV_MANIFEST:-}" && -s "$CSV_MANIFEST" ]]; then
  CSV_ARG=(-c "$CSV_MANIFEST")
  echo "[INFO] Using CSV manifest: $CSV_MANIFEST" | tee -a "$LOG"
else
  guess_csv="${BPM_MANIFEST%.bpm}.csv"
  if [[ -s "$guess_csv" ]]; then
    CSV_ARG=(-c "$guess_csv")
    echo "[INFO] Using CSV manifest (guessed): $guess_csv" | tee -a "$LOG"
  else
    echo "[WARN] No CSV manifest found; proceeding without -c" | tee -a "$LOG"
  fi
fi

# Temp dir & plugin path
SORT_TMP="$TMP_DIR/bcfsort_${RUN}"; mkdir -p "$SORT_TMP"
gtc2vcf_plugin="+gtc2vcf"
if [[ -n "${BCFTOOLS_PLUGINS:-}" && -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]]; then
  gtc2vcf_plugin="+$BCFTOOLS_PLUGINS/gtc2vcf.so"
fi

# Outputs
RAW_BCF="$VCF_DIR/cohort.gtc.$REF_BUILD.raw.bcf"
NORM_BCF="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.bcf"

# Step 1: plugin → sort → RAW_BCF
echo "[STEP] gtc2vcf -> sort (RAW BCF): $RAW_BCF" | tee -a "$LOG"
set -x
bcftools "$gtc2vcf_plugin" \
  --no-version -Ou \
  -b "$BPM_MANIFEST" "${CSV_ARG[@]}" -e "$EGT_CLUSTER" -f "$REFERENCE_FASTA" \
  -g "$GTC_DIR" 2>>"$LOG" \
| bcftools sort -Ob -T "$SORT_TMP" -o "$RAW_BCF" 2>>"$LOG"
set +x
bcftools index -f "$RAW_BCF" 2>>"$LOG" || true
echo "[OK] RAW BCF created & indexed: $RAW_BCF" | tee -a "$LOG"

# Step 2: normalize → NORM_BCF
echo "[STEP] Normalize (split, left-align, drop REF mismatches) -> NORM BCF" | tee -a "$LOG"
set -x
bcftools norm --no-version \
  -m -any -c x -f "$REFERENCE_FASTA" \
  -Ob -o "$NORM_BCF" --write-index "$RAW_BCF" 2>>"$LOG"
set +x
echo "[OK] NORM BCF created & indexed: $NORM_BCF" | tee -a "$LOG"

# Step 3: NORM_BCF → gzipped VCF
echo "[STEP] Convert NORM BCF -> gzipped VCF & index: $VCF_NORM" | tee -a "$LOG"
set -x
bcftools view -Oz -o "$VCF_NORM" "$NORM_BCF" 2>>"$LOG"
tabix -p vcf "$VCF_NORM" 2>>"$LOG"
set +x
echo "[OK] Wrote: $VCF_NORM (+ .tbi)" | tee -a "$LOG"

# Step 4: Fill tags (AN, AC, AF)
echo "[STEP] Fill tags (AN, AC, AF) -> $TAGS_BCF" | tee -a "$LOG"
set -x
bcftools view -Ou "$NORM_BCF" 2>>"$LOG" \
| bcftools +fill-tags -Ob -o "$TAGS_BCF" -- -t AN,AC,AF 2>>"$LOG"
bcftools index -f "$TAGS_BCF" 2>>"$LOG" || true
set +x
echo "[OK] Tags BCF created & indexed: $TAGS_BCF" | tee -a "$LOG"

# Summary
samples=$(bcftools query -l "$NORM_BCF" | wc -l | tr -d ' ')
records=$(bcftools index -n "$NORM_BCF")
echo "[SUMMARY] Samples: $samples" | tee -a "$LOG"
echo "[SUMMARY] Records: $records" | tee -a "$LOG"
if [[ "$samples" -ne "$gtc_n" ]]; then
  echo "[WARN] VCF samples ($samples) != GTC files ($gtc_n). Check sample sheet & GTC names." | tee -a "$LOG"
else
  echo "[OK] VCF samples match GTC count." | tee -a "$LOG"
fi

bcftools view -h "$VCF_NORM" | grep -m1 '^#CHROM' | tee -a "$LOG" || true
echo "== Done: GTC -> VCF ==" | tee -a "$LOG"
