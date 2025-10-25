#!/usr/bin/env bash
# Build PSAM + sexmap where IID == SentrixBarcode_A "_" SentrixPosition_A (matches VCF sample IDs).
set -euo pipefail

# Try to pick up repo config for sane defaults (non-fatal if missing)
if [ -z "${REPO_ROOT:-}" ]; then
  _SCRIPT="${BASH_SOURCE[0]:-$0}"
  _SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
  REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
fi
if [ -f "${REPO_ROOT}/scripts/00_config.sh" ]; then
  # shellcheck source=/dev/null
  source "${REPO_ROOT}/scripts/00_config.sh"
fi

# Positional args
SHEET="${1:?Usage: build_psam_from_barcode.sh SampleSheet.csv [metadata/vcf.samples [PSAM_OUT [SEXMAP_OUT]]]}"
VCF_SAMPLES="${2:-${REPO_ROOT:-.}/metadata/vcf.samples}"
PSAM_OUT="${3:-${PSAM_SEX:-${REPO_ROOT:-.}/metadata/cohort.sex.psam}}"
SEXMAP_OUT="${4:-${SEXMAP_PATH:-${REPO_ROOT:-.}/metadata/sexmap.txt}}"

# Sanity checks
[ -s "$SHEET" ] || { echo "[ERROR] Sample sheet not found/empty: $SHEET" >&2; exit 1; }
[ -s "$VCF_SAMPLES" ] || { echo "[ERROR] VCF sample list not found/empty: $VCF_SAMPLES" >&2; exit 1; }

mkdir -p "$(dirname "$PSAM_OUT")"
mkdir -p "$(dirname "$SEXMAP_OUT")"

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# 1) Parse [Data] section; map IID -> SEX (1=male,2=female,0=unknown)
#    - robust to CRLF and to headers like "SentrixBarcode_A", "Sentrix Position A", "Gender"/"Sex"
tr -d '\r' < "$SHEET" | awk -F',' -v IGNORECASE=1 '
  BEGIN{OFS="\t"; indata=0}
  $1 ~ /^\[Data\]$/ {indata=1; next}
  indata && !hdr {
    hdr=1
    for (i=1;i<=NF;i++){
      h=$i; gsub(/^"|"$/,"",h); n[i]=h
      # capture typical variants
      if (h ~ /sentrix.*barcode/) b=i
      if (h ~ /sentrix.*position/) p=i
      if (h ~ /^(gender|sex)$/)   s=i
    }
    if (!b || !p) { print "ERROR: Missing SentrixBarcode/Position columns" > "/dev/stderr"; exit 1 }
    next
  }
  indata && hdr {
    for (i=1;i<=NF;i++) gsub(/^"|"$/,"",$i)
    key = $b "_" $p
    g = (s? tolower($s) : "")
    sex = (g~/^m(ale)?$/)?1:((g~/^f(emale)?$/)?2:0)
    print key, sex
  }
' | sort -u > "$tmp/map.tsv"

# 2) Emit PLINK2 PSAM with header '#FID IID PAT MAT SEX PHENOTYPE'
if [ -s "$PSAM_OUT" ]; then
  echo "[NOTE] Overwriting existing PSAM: $PSAM_OUT"
fi
{
  echo -e "#FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE"
  # Build a quick assoc array in awk for speed
  awk -F'\t' 'NR==FNR{m[$1]=$2; next} NR>FNR{
    iid=$0; sex=(iid in m? m[iid]:"0");
    printf "%s\t%s\t0\t0\t%s\t-9\n", iid, iid, sex
  }' "$tmp/map.tsv" "$VCF_SAMPLES"
} > "$PSAM_OUT"

# 3) IID→SEX map derived from the PSAM itself
if [ -s "$SEXMAP_OUT" ]; then
  echo "[NOTE] Overwriting existing sexmap: $SEXMAP_OUT"
fi
awk -F'\t' 'NR>1{print $2"\t"$5}' "$PSAM_OUT" > "$SEXMAP_OUT"

echo "Wrote:"
echo "  PSAM     : $PSAM_OUT"
echo "  sexmap   : $SEXMAP_OUT"
echo "  vcf.smpl : $VCF_SAMPLES"
