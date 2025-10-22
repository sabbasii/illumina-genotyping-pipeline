#!/usr/bin/env bash
# Build PSAM + sexmap where IID == SentrixBarcode_A "_" SentrixPosition_A (matches your VCF).
set -euo pipefail

SHEET="${1:?Usage: build_psam_from_barcode.sh SampleSheet.csv}"
VCF_SAMPLES="${2:?Usage: build_psam_from_barcode.sh SampleSheet.csv metadata/vcf.samples}"
PSAM_OUT="${3:-metadata/cohort.sex.psam}"
SEXMAP_OUT="${4:-metadata/sexmap.txt}"

mkdir -p "$(dirname "$PSAM_OUT")"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# 1) Parse [Data] (robust to [Data],,,,,,), map key -> sex
tr -d '\r' < "$SHEET" | awk -F',' -v IGNORECASE=1 '
  BEGIN{OFS="\t"; indata=0}
  $1 ~ /^\[Data\]$/ {indata=1; next}
  indata && !hdr {
    hdr=1
    for (i=1;i<=NF;i++){
      h=$i; gsub(/^"|"$/,"",h); n[i]=h
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
{
  echo -e "#FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE"
  while IFS=$'\n' read -r iid; do
    sex="0"
    hit=$(awk -v id="$iid" '$1==id{print $2; exit}' "$tmp/map.tsv")
    if [[ -n "${hit:-}" ]]; then sex="$hit"; fi
    printf "%s\t%s\t0\t0\t%s\t-9\n" "$iid" "$iid" "$sex"
  done < "$VCF_SAMPLES"
} > "$PSAM_OUT"

# 3) IID→SEX map derived from the PSAM itself
awk -F'\t' 'NR>1{print $2"\t"$5}' "$PSAM_OUT" > "$SEXMAP_OUT"

echo "Wrote: $PSAM_OUT  and  $SEXMAP_OUT"
