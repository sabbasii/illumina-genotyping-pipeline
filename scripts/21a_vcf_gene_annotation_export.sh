#!/usr/bin/env bash
# vcf_gene_annotation_export.sh
#
# From an annotated VCF, this script:
#   1) Ensures a tabix index exists
#   2) Verifies that INFO/GENE is present
#   3) Exports SNP–gene–genotype tables as TSV and CSV
#
# Usage (from anywhere):
#   bash scripts/vcf_gene_annotation_export.sh
#
# Assumes:
#   - Repo layout with: output/genotype_run1/vcf/cohort.gtc.GRCh37.annotated.vcf.gz
#   - bcftools and tabix are installed and on PATH.


########################################
# Resolve repo root (one level above scripts/)
########################################
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"

VCF_DIR="$REPO_ROOT/output/genotype_run1/vcf"
ANNOTATED_VCF_BASENAME="cohort.gtc.GRCh37.annotated.vcf.gz"
ANNOTATED_VCF="$VCF_DIR/$ANNOTATED_VCF_BASENAME"
TSV_OUT="$VCF_DIR/cohort_snps_genes_genotypes.tsv"
CSV_OUT="$VCF_DIR/cohort_snps_genes_genotypes.csv"

echo "[INFO] REPO_ROOT = $REPO_ROOT"
echo "[INFO] VCF_DIR   = $VCF_DIR"
echo

########################################
# Check tools
########################################
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH."; exit 1; }
command -v tabix    >/dev/null 2>&1 || { echo "[ERROR] tabix not found on PATH."; exit 1; }

########################################
# Check annotated VCF presence
########################################
if [[ ! -s "$ANNOTATED_VCF" ]]; then
  echo "[ERROR] Annotated VCF not found:"
  echo "        $ANNOTATED_VCF"
  echo "        Make sure you ran the glist_hg19 gene annotation workflow first."
  exit 1
fi

echo "[INFO] Using annotated VCF:"
ls -lh "$ANNOTATED_VCF"
echo

cd "$VCF_DIR"

########################################
# Ensure tabix index exists
########################################
if [[ -s "${ANNOTATED_VCF_BASENAME}.tbi" ]]; then
  echo "[INFO] Found existing index: ${ANNOTATED_VCF_BASENAME}.tbi"
else
  echo "[INFO] Index not found. Building tabix index..."
  tabix -p vcf "$ANNOTATED_VCF_BASENAME"
  echo "[INFO] Created index:"
  ls -lh "${ANNOTATED_VCF_BASENAME}.tbi"
fi
echo

########################################
# Verify INFO/GENE is present in header
########################################
echo "[INFO] Checking VCF header for INFO/GENE definition..."
if bcftools view -h "$ANNOTATED_VCF_BASENAME" | grep -q 'INFO=<ID=GENE'; then
  bcftools view -h "$ANNOTATED_VCF_BASENAME" | grep 'INFO=<ID=GENE'
  echo "[INFO] INFO/GENE field detected."
else
  echo "[ERROR] No INFO/GENE definition found in header."
  echo "        The VCF does not appear to be annotated with GENE information."
  echo "        Please run the glist_hg19 gene annotation workflow first."
  exit 1
fi
echo

########################################
# Quick spot-check of variants with GENE field
########################################
echo "[INFO] Example variants (ID, CHROM, POS, REF, ALT, GENE):"
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE\n' \
  "$ANNOTATED_VCF_BASENAME" | head
echo

########################################
# Export SNP–Gene–Genotype table as TSV
########################################
echo "[INFO] Exporting SNP–Gene–Genotype table (TSV) to:"
echo "       $TSV_OUT"

bcftools query -H \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE[\t%GT]\n' \
  "$ANNOTATED_VCF_BASENAME" \
  > "$TSV_OUT"

echo "[INFO] TSV export complete. Preview (first 5 lines):"
head "$TSV_OUT" | sed -n '1,5p'
echo

########################################
# Export SNP–Gene–Genotype table as CSV
# (convert TSV → CSV for consistency with the docs)
########################################
echo "[INFO] Converting TSV → CSV:"
echo "       $CSV_OUT"

tr '\t' ',' < "$TSV_OUT" > "$CSV_OUT"

echo "[INFO] CSV export complete. Preview (first 5 lines):"
head "$CSV_OUT" | sed -n '1,5p'
echo

########################################
# Done
########################################
echo "[DONE] Generated SNP–Gene–Genotype tables:"
echo "  - $TSV_OUT"
echo "  - $CSV_OUT"
echo
echo "[NOTE] These files are ready for downstream analysis in R, Python, or Excel."
