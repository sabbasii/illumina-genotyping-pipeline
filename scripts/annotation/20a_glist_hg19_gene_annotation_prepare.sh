#!/usr/bin/env bash
# prepare_glist_hg19.sh
# Download PLINK glist-hg19, clean it, bgzip + tabix index, and create gene_header.txt
# Usage (from anywhere):
#   bash scripts/gene_annotation_prep.sh

########################################
# Resolve repo root (one level above scripts/)
########################################
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"

GENE_DIR="$REPO_ROOT/reference/gene_ranges_hg19"
mkdir -p "$GENE_DIR"
cd "$GENE_DIR"

echo "[INFO] REPO_ROOT = $REPO_ROOT"
echo "[INFO] Working in $GENE_DIR"
echo

########################################
# Check tools
########################################
have_wget=0
have_curl=0
command -v wget >/dev/null 2>&1 && have_wget=1
command -v curl >/dev/null 2>&1 && have_curl=1

command -v bgzip >/dev/null 2>&1 || { echo "[ERROR] bgzip not found (usually from htslib/tabix)."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "[ERROR] tabix not found."; exit 1; }

########################################
# Step 1 — Download glist-hg19 (if missing)
########################################
GLIST_RAW="glist-hg19"

if [[ -s "$GLIST_RAW" ]]; then
  echo "[INFO] Found existing $GLIST_RAW, skipping download."
else
  echo "[INFO] Downloading glist-hg19 from PLINK resources..."
  URL="https://www.cog-genomics.org/static/bin/plink/glist-hg19"

  if (( have_wget )); then
    wget -O "$GLIST_RAW" "$URL"
  elif (( have_curl )); then
    curl -L -o "$GLIST_RAW" "$URL"
  else
    echo "[ERROR] Neither wget nor curl is available. Please download glist-hg19 manually:"
    echo "        $URL"
    echo "        and place it as: $GENE_DIR/$GLIST_RAW"
    exit 1
  fi
fi

echo "[INFO] Preview of glist-hg19:"
head "$GLIST_RAW"
wc -l "$GLIST_RAW"
echo

########################################
# Step 2 — Ensure sorted by CHR, START
########################################
echo "[INFO] Checking sort order..."
if sort -k1,1 -k2,2n "$GLIST_RAW" | diff -q - "$GLIST_RAW" >/dev/null; then
  echo "[INFO] File is already sorted by (CHR, START)."
else
  echo "[INFO] Sorting by chromosome and start position..."
  sort -k1,1 -k2,2n "$GLIST_RAW" > "${GLIST_RAW}.sorted"
  mv "${GLIST_RAW}.sorted" "$GLIST_RAW"
fi
echo

########################################
# Step 3 — Force TAB-delimited 4-column format
########################################
echo "[INFO] Rebuilding as TAB-delimited with 4 columns (CHR, START, END, GENE)..."
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' "$GLIST_RAW" > "${GLIST_RAW}.tab"
mv "${GLIST_RAW}.tab" "$GLIST_RAW"

echo "[INFO] First few lines with visible TABs (^I):"
head "$GLIST_RAW" | cat -t
echo

########################################
# Step 4 — bgzip + tabix index
########################################
echo "[INFO] Compressing with bgzip..."
bgzip -f "$GLIST_RAW"   # produces glist-hg19.gz

echo "[INFO] Building tabix index (CHR=1, START=2, END=3)..."
tabix -s 1 -b 2 -e 3 "${GLIST_RAW}.gz"

echo "[INFO] Compressed + indexed files:"
ls -lh "${GLIST_RAW}.gz" "${GLIST_RAW}.gz.tbi"
echo

########################################
# Step 5 — Create gene_header.txt for INFO/GENE
########################################
HEADER_FILE="gene_header.txt"
echo "[INFO] Writing VCF header snippet to $HEADER_FILE"

cat << 'EOF' > "$HEADER_FILE"
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name from PLINK glist-hg19 annotation">
EOF

echo "[INFO] Header file content:"
cat "$HEADER_FILE"
echo

########################################
# Done
########################################
echo "[DONE] Prepared gene range resources:"
echo "  - $(realpath "${GLIST_RAW}.gz")"
echo "  - $(realpath "${GLIST_RAW}.gz.tbi")"
echo "  - $(realpath "$HEADER_FILE")"
echo
echo "These can now be used with bcftools annotate, e.g.:"
echo "  bcftools annotate \\"
echo "    -a reference/gene_ranges_hg19/glist-hg19.gz \\"
echo "    -h reference/gene_ranges_hg19/gene_header.txt \\"
echo "    -c CHROM,FROM,TO,GENE \\"
echo "    -Oz -o cohort.gtc.GRCh37.annotated.vcf.gz \\"
echo "    cohort.gtc.GRCh37.norm.vcf.gz"
