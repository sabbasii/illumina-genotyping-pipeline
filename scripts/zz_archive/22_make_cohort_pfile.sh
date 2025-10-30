#!/usr/bin/env bash
# Make a cohort-restricted, phenotyped pfile from analysis.clean using your keep list + PHENO PSAM.

set -euo pipefail

# --- load repo config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- inputs (from prior steps)
PHENO_DIR="${PHENO_DIR:-$REPO_ROOT/metadata/pheno}"
PHENO_PSAM="${PHENO_PSAM:-$PHENO_DIR/cohort.pheno.psam}"
IID_KEEP="$PHENO_DIR/iid_selected.keep"

EXPR_OUT_DIR="${EXPR_OUT_DIR:-$OUT_DIR/expr/explore}"
IDS_SELECTED="$EXPR_OUT_DIR/ids_selected.txt"   # UASG list (224)

# --- output prefix
OUT_PREFIX="${PLINK_DIR}/analysis.clean.expr_overlap"   # change if you want another name

echo "=== Cohort Pfile Build ==="
echo "RUN=$RUN   REF_BUILD=$REF_BUILD"
echo "PLINK_DIR=$PLINK_DIR"
echo "PHENO_PSAM=$PHENO_PSAM"
echo "IID_KEEP=$IID_KEEP"
echo "OUT_PREFIX=$OUT_PREFIX"
echo

# --- quick checks
[[ -s "$PLINK_DIR/analysis.clean.psam" ]] || { echo "[ERR] missing: $PLINK_DIR/analysis.clean.psam"; exit 1; }
[[ -s "$PHENO_PSAM" ]] || { echo "[ERR] missing: $PHENO_PSAM"; exit 1; }
[[ -s "$IID_KEEP" ]] || { echo "[ERR] missing: $IID_KEEP"; exit 1; }
[[ -s "$IDS_SELECTED" ]] || { echo "[ERR] missing: $IDS_SELECTED"; exit 1; }

echo "[CHECK] PSAM rows (clean set):"
awk 'END{print NR-1}' "$PLINK_DIR/analysis.clean.psam"

echo "[CHECK] IIDs in keep list (phenotyped overlap):"
wc -l < "$IID_KEEP"

echo "[CHECK] Intended UASGs (ids_selected.txt):"
wc -l < "$IDS_SELECTED"

# UASGs that actually have PHENO1=1/2 in the PHENO PSAM
UASG_PHENO_TXT="$PHENO_DIR/uasg_with_pheno.txt"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="UASG") u=i; next}
            NR>1 && $NF ~ /^[12]$/{print $u}' "$PHENO_PSAM" | sort -u > "$UASG_PHENO_TXT"

echo "[CHECK] UASGs with PHENO1 (1/2) in PHENO PSAM:"
wc -l < "$UASG_PHENO_TXT"

# Show any selected UASGs missing a phenotype (should be just UASG-0337)
echo "[CHECK] Selected UASGs missing PHENO1:"
comm -23 <(sort "$IDS_SELECTED") <(cat "$UASG_PHENO_TXT") || true
echo

# --- plink2 subset: use updated PHENO PSAM and keep-list
echo "[RUN] plink2 --pfile analysis.clean --psam PHENO_PSAM --keep iid_selected.keep → ${OUT_PREFIX}.*"
plink2 \
  --pfile "$PLINK_DIR/analysis.clean" \
  --psam  "$PHENO_PSAM" \
  --keep  "$IID_KEEP" \
  --make-pgen \
  --out   "$OUT_PREFIX"

echo
echo "[DONE] Wrote cohort-restricted pfile:"
echo "  ${OUT_PREFIX}.pgen"
echo "  ${OUT_PREFIX}.pvar"
echo "  ${OUT_PREFIX}.psam"


# -----------------RUN-----------------
# source scripts/00_config.sh
# bash scripts/22_make_cohort_pfile.sh
