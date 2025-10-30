#!/usr/bin/env bash
set -euo pipefail
PSAM="${1:?Usage: 21a_psam_validate_header.sh /path/to/file.psam}"
echo "[INFO] Checking header of: $PSAM"

# show first line with visible tabs
head -1 "$PSAM" | sed 's/\t/|/g'

# fix: ensure first header token is '#FID' (not 'FID'), and strip CRs
tmp="$(mktemp)"
awk 'NR==1{ if($1=="FID") $1="#FID"; } {print}' OFS="\t" "$PSAM" | sed 's/\r$//' > "$tmp"
mv "$tmp" "$PSAM"

# verify field counts
hNF=$(awk -F"\t" 'NR==1{print NF}' "$PSAM")
bad=$(awk -F"\t" -v H="$hNF" 'NF!=H{print NR}' "$PSAM" | wc -l)
echo "[INFO] Header fields: $hNF | lines with mismatched field count: $bad"

# --------------- RUN --------------
# bash scripts/21a_psam_validate_header.sh "$PHENO_PSAM"
