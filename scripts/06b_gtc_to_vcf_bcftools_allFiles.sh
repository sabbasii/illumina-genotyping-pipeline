#!/usr/bin/env bash
# 06_gtc_to_vcf_all_csv_wide.sh
# - Convert ALL .gtc → VCF/BCF under output/$RUN/vcf-all
# - Include FORMAT: GT,GQ,IGC,BAF,LRR
# - Normalize (split multiallelics, left-align, drop REF mismatches)
# - Fill AN/AC/AF
# - Export WIDE CSV matrices (variants as rows, samples as columns):
#     matrix_GT_{raw,norm}.csv          (GT protected as ="0/1" to prevent Excel date auto-convert)
#     matrix_GQ_{raw,norm}.csv
#     matrix_IGC_{raw,norm}.csv
#     matrix_BAF_{raw,norm}.csv
#     matrix_LRR_{raw,norm}.csv
#   Plus: sites_with_AF.csv
# - Self-checks read data files (avoid “.csi index alone” warning)

set -euo pipefail

# ---------- Resolve repo root & load config ----------
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

activate_env || true
ensure_dirs

# ---------- Output dirs ----------
VCF_ALL_DIR="$REPO_ROOT/output/$RUN/vcf-all"
mkdir -p "$VCF_ALL_DIR" "$TMP_DIR" "$LOG_DIR"

LOG="$LOG_DIR/06_gtc_to_vcf_all_csv_wide.log"
exec > >(tee "$LOG") 2>&1

echo "== GTC -> VCF (ALL) with GT,GQ,IGC,BAF,LRR; WIDE CSV exports =="
echo "RUN=$RUN  REF_BUILD=$REF_BUILD"
echo "GTC_DIR=$GTC_DIR"
echo "VCF_ALL_DIR=$VCF_ALL_DIR"
echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"
echo

# ---------- Tool checks ----------
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found"; exit 1; }
command -v tabix    >/dev/null 2>&1 || { echo "[ERROR] tabix not found"; exit 1; }
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
    echo "[WARN] samtools not found; ${REFERENCE_FASTA}.fai missing"
  fi
fi

# Optional CSV manifest
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

# ---------- Paths ----------
SORT_TMP="$TMP_DIR/bcfsort_${RUN}"; mkdir -p "$SORT_TMP"
gtc2vcf_plugin="+gtc2vcf"
if [[ -n "${BCFTOOLS_PLUGINS:-}" && -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]]; then
  gtc2vcf_plugin="+$BCFTOOLS_PLUGINS/gtc2vcf.so"
fi

RAW_BCF="$VCF_ALL_DIR/cohort.gtc.$REF_BUILD.raw.bcf"
NORM_BCF="$VCF_ALL_DIR/cohort.gtc.$REF_BUILD.norm.bcf"
VCF_GZ="$VCF_ALL_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"
TAGS_BCF="$VCF_ALL_DIR/cohort.gtc.$REF_BUILD.tags.bcf"

AF_CSV="$VCF_ALL_DIR/sites_with_AF.csv"

# Wide matrices:
M_RAW_GT="$VCF_ALL_DIR/matrix_GT_raw.csv"
M_NRM_GT="$VCF_ALL_DIR/matrix_GT_norm.csv"
M_RAW_GQ="$VCF_ALL_DIR/matrix_GQ_raw.csv"
M_NRM_GQ="$VCF_ALL_DIR/matrix_GQ_norm.csv"
M_RAW_IGC="$VCF_ALL_DIR/matrix_IGC_raw.csv"
M_NRM_IGC="$VCF_ALL_DIR/matrix_IGC_norm.csv"
M_RAW_BAF="$VCF_ALL_DIR/matrix_BAF_raw.csv"
M_NRM_BAF="$VCF_ALL_DIR/matrix_BAF_norm.csv"
M_RAW_LRR="$VCF_ALL_DIR/matrix_LRR_raw.csv"
M_NRM_LRR="$VCF_ALL_DIR/matrix_LRR_norm.csv"

# ---------- Convert ALL GTC -> RAW BCF (include GT,GQ,IGC,BAF,LRR) ----------
echo "[STEP] gtc2vcf -> sort (RAW BCF): $RAW_BCF"
set -x
bcftools "$gtc2vcf_plugin" \
  --no-version -Ou \
  -b "$BPM_MANIFEST" "${CSV_ARG[@]}" -e "$EGT_CLUSTER" -f "$REFERENCE_FASTA" \
  -g "$GTC_DIR" \
  -t GT,GQ,IGC,BAF,LRR \
| bcftools sort -Ob -T "$SORT_TMP" -o "$RAW_BCF"
set +x
bcftools index -f "$RAW_BCF" || true
echo "[OK] RAW BCF created & indexed: $RAW_BCF"
echo

# ---------- Normalize ----------
echo "[STEP] Normalize -> NORM BCF"
set -x
bcftools norm --no-version \
  -m -any -c x -f "$REFERENCE_FASTA" \
  -Ob -o "$NORM_BCF" --write-index "$RAW_BCF"
set +x
echo "[OK] NORM BCF created & indexed: $NORM_BCF(.csi)"
echo

# ---------- NORM BCF -> VCF.GZ ----------
echo "[STEP] NORM BCF -> VCF.GZ"
set -x
bcftools view -Oz -o "$VCF_GZ" "$NORM_BCF"
tabix -f -p vcf "$VCF_GZ"
set +x
echo "[OK] Wrote: $VCF_GZ (+ .tbi)"
echo

# ---------- Helper: write header row (variant id + samples) ----------
write_header() {
  local bcf="$1" out="$2"
  {
    # left side
    printf "CHROM,POS,ID,REF,ALT"
    # sample names
    while read -r s; do
      printf ",%s" "$s"
    done < <(bcftools query -l "$bcf")
    printf "\n"
  } > "$out"
}

# ---------- Helper: export a field into wide CSV ----------
# Args: bcf path, output.csv, FORMAT_field, protect_gt (0/1)
export_field_wide() {
  local bcf="$1" out="$2" fld="$3" protect="$4"

  write_header "$bcf" "$out"

  # Build the bcftools query format string safely: %CHROM ... [\t%GT] (or %GQ, %IGC, ...).
  # NOTE: we must not write %%$fld — that causes the error you saw.
  local fmt="%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%${fld}]\n"

  bcftools query -f "$fmt" "$bcf" \
  | awk -v PF="$protect" -v OFS=',' '
      BEGIN{FS="\t"}
      {
        # first 5 variant columns
        printf "%s,%s,%s,%s,%s", $1,$2,$3,$4,$5
        # per-sample values
        for(i=6;i<=NF;i++){
          val=$i
          if(PF=="1"){                 # Protect GT: force Excel to treat as text
            if(val=="" || val==".") printf ","
            else printf ",=\"%s\"", val
          } else {
            printf ",%s", val
          }
        }
        printf "\n"
      }' >> "$out"
  echo "[OK] Wrote ${fld} → ${out}"
}


echo "[STEP] Export WIDE matrices (RAW)"
export_field_wide "$RAW_BCF" "$M_RAW_GT"  "GT"  1
export_field_wide "$RAW_BCF" "$M_RAW_GQ"  "GQ"  0
export_field_wide "$RAW_BCF" "$M_RAW_IGC" "IGC" 0
export_field_wide "$RAW_BCF" "$M_RAW_BAF" "BAF" 0
export_field_wide "$RAW_BCF" "$M_RAW_LRR" "LRR" 0
echo

echo "[STEP] Export WIDE matrices (NORM)"
export_field_wide "$NORM_BCF" "$M_NRM_GT"  "GT"  1
export_field_wide "$NORM_BCF" "$M_NRM_GQ"  "GQ"  0
export_field_wide "$NORM_BCF" "$M_NRM_IGC" "IGC" 0
export_field_wide "$NORM_BCF" "$M_NRM_BAF" "BAF" 0
export_field_wide "$NORM_BCF" "$M_NRM_LRR" "LRR" 0
echo

# ---------- Export site-level AF CSV ----------
echo "[STEP] Export site-level frequency table (AN/AC/AF) to CSV"
echo "CHROM,POS,ID,REF,ALT,AN,AC,AF" > "$AF_CSV"
bcftools query -f '%CHROM,%POS,%ID,%REF,%ALT,%INFO/AN,%INFO/AC,%INFO/AF\n' "$TAGS_BCF" >> "$AF_CSV" || true
echo "[OK] Site AF CSV: $AF_CSV"
echo

# ---------- Self-checks ----------
echo "[CHECK] Contigs in NORM header:"
bcftools view -h "$NORM_BCF" | grep '^##contig' | head || true
echo
echo "[CHECK] FORMAT fields in RAW header (expect GT,GQ,IGC,BAF,LRR):"
bcftools view -h "$RAW_BCF" | grep '^##FORMAT' || true
echo

# ---------- Summaries ----------
samples_raw=$(bcftools query -l "$RAW_BCF" | wc -l | tr -d ' ')
samples_norm=$(bcftools query -l "$NORM_BCF" | wc -l | tr -d ' ')
records_norm=$(bcftools index -n "$NORM_BCF")
echo "[SUMMARY] RAW_BCF samples  : $samples_raw"
echo "[SUMMARY] NORM_BCF samples : $samples_norm"
echo "[SUMMARY] NORM_BCF records : $records_norm"
bcftools view -h "$VCF_GZ" | grep -m1 '^#CHROM' || true
echo
echo "== Done: ALL GTC -> vcf-all with WIDE CSV matrices (GT protected) =="
