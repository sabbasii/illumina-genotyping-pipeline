#!/bin/bash
set -euo pipefail

############################
# Colors (optional)
############################
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m'

############################
# Resolve base directory
############################
# Find the directory of this script, change into it, get its absolute path, and store it in BASE_DIR.
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"  # run pwd only if cd succeeds

############################
# Input files
############################
TARGET_SNP_LIST="${BASE_DIR}/target_snps.txt"
GENOTYPE_BCF="${BASE_DIR}/genotypes_processed.bcf"
REF_SNP_LIST="${BASE_DIR}/dbsnp_153common.bed"

############################
# checking required input files
############################
for f in "$TARGET_SNP_LIST" "$GENOTYPE_BCF" "$REF_SNP_LIST"; do
  [[ -f "$f" ]] || { echo -e "${RED}Missing file:${NC} $f"; exit 1; }
done

[[ -f "${GENOTYPE_BCF}.csi" ]] || {
  echo -e "${RED}Missing index:${NC} ${GENOTYPE_BCF}.csi"
  exit 1
}

############################
# Output directories
############################
ANALYSIS_DIR="${BASE_DIR}/analysis_target_snp"
LOG_DIR="${ANALYSIS_DIR}/log_files"

mkdir -p "$ANALYSIS_DIR" "$LOG_DIR"

LOG_FILE="${LOG_DIR}/target_snp_report.log"
OUTPUT_FILE="${ANALYSIS_DIR}/target_snp_report.tsv"

############################
# Threads
# Get the number of available CPU cores and calculate 80% of that number for threading.
############################
CORES=$(nproc --all)
THREADS=$(( CORES * 80 / 100 ))
[[ $THREADS -lt 1 ]] && THREADS=1  # Ensure at least 1 thread

############################
# Create report header
############################
bcftools view -h "$GENOTYPE_BCF" | tail -n 1 > "$OUTPUT_FILE"

############################
# Main loop
############################
not_found=0

while read -r line; do
  target="${line//[[:space:]]/}"

  echo -e "Searching for target SNP: ${YELLOW}$target${NC}" | tee -a "$LOG_FILE"

  mapfile -t -d $'\t' snp_line < <(grep -w "$target" "$REF_SNP_LIST")

  chr="${snp_line[0]}"
  start="${snp_line[1]}"
  end="${snp_line[2]}"
  ref_ID="${snp_line[3]//$'\n'/}"

  if [[ -z "$ref_ID" ]]; then
    echo -e "\t${RED}Not found in reference file${NC}\n" | tee -a "$LOG_FILE"
    ((not_found++))
    continue
  fi

  start_pos=$((start - 1))
  end_pos=$((end + 1))
  region="$chr:$start_pos-$end_pos"

  mapfile -t marker_arr < <(
    bcftools query -r "$region" -f '%CHROM %POS %ID\n' "$GENOTYPE_BCF"
  )

  rsID="$(echo "${marker_arr[0]}" | awk '{print $3}')"  

  if [[ -z "$rsID" ]]; then
    echo -e "\t${RED}No genotype found at $region${NC}\n" | tee -a "$LOG_FILE"
    ((not_found++))
    continue
  fi

  bcftools norm \
    --multiallelics - \
    --do-not-normalize \
    -r "$region" \
    "$GENOTYPE_BCF" |
    bcftools view -H |
    sed "s/$rsID/$ref_ID/g" \
    >> "$OUTPUT_FILE"

  echo -e "\t${GREEN}Added to report${NC}\n" | tee -a "$LOG_FILE"

done < "$TARGET_SNP_LIST"

############################
# Final message
############################
echo -e "${GREEN}Done.${NC} Output:"
echo -e "${YELLOW}$OUTPUT_FILE${NC}"

