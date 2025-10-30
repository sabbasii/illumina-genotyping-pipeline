#!/usr/bin/env bash
# 09_subset_vcf_to_cohort.sh — build cohort-only VCF using ids_selected.txt (UASGs) → map to IIDs → subset VCF
# Run order: after 03_gtc_to_vcf_bcftools.sh (VCF_NORM ready), before 10_qc_vcf.sh
# Usage:
#   source scripts/00_config.sh
#   bash scripts/09_subset_vcf_to_cohort.sh [optional:/path/to/input.vcf.gz]

set -euo pipefail

# --- Resolve repo root & load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

command -v bcftools >/dev/null 2>&1 || { echo "[ERR] bcftools not found"; exit 1; }

# --- Inputs
VCF_IN="${1:-$VCF_NORM}"                   # normalized, gzipped VCF from step 03
IDS_SELECTED="$EXPR_OUT_DIR/ids_selected.txt"   # UASG list from 02_overlap_and_select_cohorts.py
SHEET="${SHEET_PATH:-$SAMPLE_SHEET}"            # SampleSheet CSV

[[ -s "$VCF_IN" ]] || { echo "[ERR] VCF not found: $VCF_IN"; exit 1; }
[[ -s "$IDS_SELECTED" ]] || { echo "[ERR] Selected UASG list not found: $IDS_SELECTED"; exit 1; }
[[ -s "$SHEET" ]] || { echo "[ERR] SampleSheet not found: $SHEET"; exit 1; }

ensure_dirs
mkdir -p "$META_DIR" "$PHENO_DIR"

# --- Map IID <-> UASG from SampleSheet ([Data] section, case-insensitive)
MAP_TSV="$META_DIR/snp_iid_uasg.tsv"   # IID \t UASG
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

tr -d '\r' < "$SHEET" | awk -F',' -v IGNORECASE=1 '
  BEGIN{OFS="\t"; indata=0}
  $1 ~ /^\[Data\]$/ {indata=1; next}
  indata && !hdr {
    hdr=1
    for (i=1;i<=NF;i++){
      h=$i; gsub(/^"|"$/,"",h); n[i]=h
      if (h ~ /sentrix.*barcode/) b=i
      if (h ~ /sentrix.*position/) p=i
      if (h ~ /^(sample[_ ]?name|sampleid)$/) s=i
    }
    if (!b || !p || !s) { print "ERROR: Missing SentrixBarcode/Position/Sample_Name columns" > "/dev/stderr"; exit 1 }
    next
  }
  indata && hdr {
    for (i=1;i<=NF;i++) gsub(/^"|"$/,"",$i)
    iid = $b "_" $p
    uasg = $s
    print iid, uasg
  }
' | sort -u > "$MAP_TSV"

# --- Build keep list of IIDs from selected UASGs
IID_KEEP="$PHENO_DIR/iid_selected.keep"
join -1 1 -2 2 \
  <(sort "$IDS_SELECTED") \
  <(sort -k2,2 "$MAP_TSV") \
  | awk '{print $2}' > "$IID_KEEP"

sel_uasg=$(wc -l < "$IDS_SELECTED" | tr -d ' ')
keep_iid=$(wc -l < "$IID_KEEP" | tr -d ' ')
echo "[INFO] Selected UASGs: $sel_uasg"
echo "[INFO] Mapped IIDs   : $keep_iid"
if [[ "$keep_iid" -eq 0 ]]; then
  echo "[ERR] No IIDs mapped from selected UASGs. Check SampleSheet columns and ids_selected.txt"
  exit 1
fi

# Report any selected UASGs that failed to map
comm -23 <(sort "$IDS_SELECTED") <(cut -f2 "$MAP_TSV" | sort -u) | sed 's/^/[WARN] Unmapped UASG: /' || true

# --- Subset VCF to cohort IIDs
VCF_OUT="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.cohort.vcf.gz"
echo "[RUN] bcftools view -S $IID_KEEP -Oz -o $VCF_OUT $VCF_IN"
bcftools view -S "$IID_KEEP" -Oz -o "$VCF_OUT" "$VCF_IN"
tabix -p vcf "$VCF_OUT"

# --- Sanity checks
all_in=$(bcftools query -l "$VCF_IN" | wc -l | tr -d ' ')
cohort_n=$(bcftools query -l "$VCF_OUT" | wc -l | tr -d ' ')
echo "[OK] Input VCF samples : $all_in"
echo "[OK] Cohort VCF samples: $cohort_n (should match mapped IIDs: $keep_iid)"

# Save the final sample list for downstream (10_qc_vcf.sh)
bcftools query -l "$VCF_OUT" > "$META_DIR/vcf.samples.cohort"

echo
echo "[DONE] Cohort-only VCF written:"
echo "  $VCF_OUT"
echo "  $VCF_OUT.tbi"
echo "[INFO] IID keep list  : $IID_KEEP"
echo "[INFO] Cohort samples : $META_DIR/vcf.samples.cohort"
