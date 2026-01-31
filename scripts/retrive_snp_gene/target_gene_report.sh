#!/bin/bash
set -euo pipefail

# Generate Genotype Report for Target Genes

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
TARGET_GENE_LIST="$REPO_ROOT/input_data/target_lists/PDE_genes.txt"
GENOTYPE_BCF="$REPO_ROOT/output/genotype_run1/imputed/genotypes_processed.bcf"
REF_GENE_LIST="$REPO_ROOT/reference/ref_gene_list/sorted_gene_list_GRCh37.bed"

############################
# Default buffer values (bp)
############################
DEFAULT_BP_BEFORE=250
DEFAULT_BP_AFTER=250

############################
# Read in arguments (optional)
############################
# Usage:
#   script.sh [bp_before_gene] [bp_after_gene]
bp_before_gene="${1:-$DEFAULT_BP_BEFORE}"
bp_after_gene="${2:-$DEFAULT_BP_AFTER}"

# Validate that values are non-negative integers
if ! [[ "$bp_before_gene" =~ ^[0-9]+$ && "$bp_after_gene" =~ ^[0-9]+$ ]]; then
  echo -e "${RED}Error:${NC} Buffer values must be non-negative integers."
  exit 1
fi

############################
# Output directories
############################

ANALYSIS_DIR="${BASE_DIR}/analysis_target_gene"
LOG_DIR="${ANALYSIS_DIR}/log_files"
mkdir -p "$ANALYSIS_DIR" "$LOG_DIR"

LOG_FILE="${LOG_DIR}/target_gene_report.log"
OUTPUT_FILE="${ANALYSIS_DIR}/target_gene_report.tsv"

# Start fresh log
: > "$LOG_FILE"

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
# Instantiate array for gene locus information
############################
gene_positions=()
not_found=0

# Read the target gene list line by line.
# - Skips empty lines and comment lines starting with "#"
while IFS= read -r line || [[ -n "$line" ]]; do
  # Trim leading/trailing whitespace
  line="$(echo "$line" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"

  # Skip blanks/comments
  [[ -z "$line" ]] && continue
  [[ "$line" =~ ^# ]] && continue

  # Strip ALL internal spaces from gene name (your original behavior)
  target="${line// /}"

  # Find the first matching row in the BED reference.
  # Assumes BED columns: chr  start  end  gene
  # Uses awk for safe, exact match on 4th column (gene).
  gene_row="$(awk -v g="$target" '($4==g){print; exit}' "$REF_GENE_LIST" || true)"

  if [[ -n "$gene_row" ]]; then
    # Split into array
    read -r chr start end gene <<< "$gene_row"

    # Compute padded positions
    start_pos=$(( start - bp_before_gene ))
    end_pos=$(( end + bp_after_gene ))

    # Clamp start to >= 1 (bcftools expects 1-based regions; negative/0 breaks)
    if [[ $start_pos -lt 1 ]]; then
      start_pos=1
    fi

    position="${chr}:${start_pos}-${end_pos}"
    gene_positions+=("$position")

    echo -e "${GREEN}${gene}${NC} found in reference. Using region: ${YELLOW}${position}${NC}" \
      | tee -a "$LOG_FILE"
  else
    echo -e "${RED}${target}${NC} not found in reference. Try another gene name." \
      | tee -a "$LOG_FILE"
    ((not_found++))
  fi
done < "$TARGET_GENE_LIST"

if [[ $not_found -gt 0 ]]; then
  echo -e "${RED}${not_found}${NC} target gene(s) were not found in the reference gene list." \
    | tee -a "$LOG_FILE"

  read -rp "Do you want to continue generating the report with the found genes? (y/n): " -n 1 -r
  echo "" | tee -a "$LOG_FILE" >/dev/null

  if [[ "$REPLY" == "y" || "$REPLY" == "Y" ]]; then
    echo -e "${GREEN}Continuing with report generation...${NC}" | tee -a "$LOG_FILE"
  else
    echo -e "${RED}Report generation aborted by user.${NC}" | tee -a "$LOG_FILE"
    exit 1
  fi
fi

if [[ ${#gene_positions[@]} -eq 0 ]]; then
  echo -e "${RED}No valid gene regions were found. Nothing to report.${NC}" | tee -a "$LOG_FILE"
  exit 1
fi

# Assemble target gene positions into a comma-separated list.
positions_list="$(IFS=,; echo "${gene_positions[*]}")"

############################
# Extract variants in regions and write report body
############################
# - split multiallelics into biallelic records for easier downstream use
# - -H: no header lines (we already wrote one line)
bcftools norm \
  --multiallelics - \
  --do-not-normalize \
  --regions "$positions_list" \
  --threads "$THREADS" \
  "$GENOTYPE_BCF" \
  2>> "$LOG_FILE" \
| bcftools view \
  -H \
  2>> "$LOG_FILE" \
>> "$OUTPUT_FILE"

############################
# Report output path
############################
echo -e "${GREEN}Target gene report complete.${NC} Output file is located at: ${YELLOW}${OUTPUT_FILE}${NC}" \
  | tee -a "$LOG_FILE"

############################
# Optional: call R script to clean the data (only if present)
############################
CLEAN_R="${BASE_DIR}/clean_report.r"
if [[ -f "$CLEAN_R" ]]; then
  Rscript "$CLEAN_R" "$OUTPUT_FILE" 2>> "$LOG_FILE"
  echo -e "${GREEN}Cleaning complete.${NC}" | tee -a "$LOG_FILE"
else
  echo -e "${YELLOW}Note:${NC} clean_report.r not found at ${CLEAN_R}. Skipping cleaning step." \
    | tee -a "$LOG_FILE"
fi

############################
# Usage
############################
# Generate a genotype report for regions surrounding target genes.
#
# Syntax:
#   mdir/target_gene_report.sh <bp_before_gene> <bp_after_gene>
#
# Arguments:
#   <bp_before_gene>   Number of base pairs upstream of the gene start
#   <bp_after_gene>    Number of base pairs downstream of the gene end
#
# Example:
#   mdir/target_gene_report.sh 10000 5000
#
# This will extract all variants located:
#   - 10 kb upstream of each target gene start
#   - 5 kb downstream of each target gene end
#
# Inputs:
#   - target_genes.txt                  (one gene symbol per line)
#   - sorted_gene_list_GRCh37.bed       (chr, start, end, gene)
#   - genotypes_processed.bcf
#
# Outputs:
#   - analysis_target_gene/target_gene_report.tsv
#   - analysis_target_gene/log_files/target_gene_report.log
############################