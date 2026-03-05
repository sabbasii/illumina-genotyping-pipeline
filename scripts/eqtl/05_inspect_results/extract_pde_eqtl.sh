#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# scripts/eqtl/05_inspect_results/06_extract_pde_eqtl.sh
#
# Purpose
#   Filter Matrix-eQTL outputs to PDE genes detected in the human
#   expression dataset, and produce quick summaries:
#     1) output/eqtl/results/inspect/pde/eqtl_cis_pde.tsv
#     2) output/eqtl/results/inspect/pde/eqtl_all_pde.tsv
#     3) counts + strongest-signal snapshots printed to stdout
#
# Inputs (default)
#   output/eqtl/results/inspect/pde/pde_genes_detected_in_expression.txt
#   output/eqtl/results/eqtl_cis.tsv
#   output/eqtl/results/eqtl_all.tsv
#
# Notes
#   - Hard-fails if the filtered outputs are empty (header-only).
#   - Assumes MatrixEQTL result tables have columns:
#       SNP (col1), gene (col2), then p-value somewhere (commonly col3).
#     We filter by gene column (col2) using exact matches to the PDE list.
#
# Run
#   source scripts/00_config.sh
#   bash scripts/eqtl/05_inspect_results/extract_pde_eqtl.sh
# ============================================================

# Resolve repo root from env (must source scripts/00_config.sh)
: "${REPO_ROOT:?REPO_ROOT is not set. Did you source scripts/00_config.sh?}"

PDE_LIST="$REPO_ROOT/output/eqtl/results/inspect/pde/pde_genes_detected_in_expression.txt"
CIS_IN="$REPO_ROOT/output/eqtl/results/eqtl_cis.tsv"
ALL_IN="$REPO_ROOT/output/eqtl/results/eqtl_all.tsv"

OUTDIR="$REPO_ROOT/output/eqtl/results/inspect/pde"
CIS_OUT="$OUTDIR/eqtl_cis_pde.tsv"
ALL_OUT="$OUTDIR/eqtl_all_pde.tsv"

mkdir -p "$OUTDIR"

# ------------------------
# Safety checks
# ------------------------
[ -s "$PDE_LIST" ] || { echo "[ERR] Missing/empty PDE list: $PDE_LIST" >&2; exit 1; }
[ -s "$CIS_IN" ]   || { echo "[ERR] Missing/empty input: $CIS_IN" >&2; exit 1; }
[ -s "$ALL_IN" ]   || { echo "[ERR] Missing/empty input: $ALL_IN" >&2; exit 1; }

# ------------------------
# Filter helper
# ------------------------
filter_eqtl_by_pde() {
  local IN="$1"
  local OUT="$2"
  local LABEL="$3"

  awk -F'\t' -v OFS='\t' '
    NR==FNR { g[$1]=1; next }   # PDE genes set
    FNR==1 { print; next }      # header
    ($2 in g)                   # keep if gene (col2) is PDE
  ' "$PDE_LIST" "$IN" > "$OUT"

  # Ensure output exists + has data rows beyond header
  [ -s "$OUT" ] || { echo "[ERR] Wrote empty file: $OUT" >&2; exit 1; }
  local nlines
  nlines="$(wc -l < "$OUT" | tr -d ' ')"
  if [ "$nlines" -le 1 ]; then
    echo "[ERR] $LABEL: No PDE rows found (header only): $OUT" >&2
    echo "      Check that column 2 is the gene column and symbols match PDE_LIST." >&2
    exit 1
  fi
}

# ------------------------
# Extract PDE rows
# ------------------------
filter_eqtl_by_pde "$CIS_IN" "$CIS_OUT" "CIS"
filter_eqtl_by_pde "$ALL_IN" "$ALL_OUT" "ALL"

cis_rows="$(( $(wc -l < "$CIS_OUT") - 1 ))"
all_rows="$(( $(wc -l < "$ALL_OUT") - 1 ))"

echo "[OK] Wrote:"
echo "  $CIS_OUT  (rows=$cis_rows)"
echo "  $ALL_OUT  (rows=$all_rows)"

# ------------------------
# Summaries
# ------------------------
echo
echo "=== Association counts ==="
echo "CIS total PDE associations: $cis_rows"
echo "ALL total PDE associations: $all_rows"

echo
echo "Top PDE genes by #cis associations:"
awk -F'\t' 'NR>1{c[$2]++} END{for(g in c) print c[g], g}' "$CIS_OUT" \
  | sort -nr | head -20

echo
echo "PDE genes WITH at least one cis-eQTL:"
awk -F'\t' 'NR>1{seen[$2]=1} END{for(g in seen) print g}' "$CIS_OUT" \
  | sort

echo
echo "=== Top 20 cis PDE hits (smallest p-values; assumes p-value is column 3) ==="
(head -1 "$CIS_OUT"; tail -n +2 "$CIS_OUT" | sort -t $'\t' -k3,3g | head -20)

echo
echo "=== Top 20 all PDE hits (smallest p-values; assumes p-value is column 3) ==="
(head -1 "$ALL_OUT"; tail -n +2 "$ALL_OUT" | sort -t $'\t' -k3,3g | head -20)