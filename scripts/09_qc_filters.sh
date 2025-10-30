#!/usr/bin/env bash
# 09_qc_filters.sh — apply SNP/sample QC, dedup, LD prune, PCA, KING relatedness
# Inputs: cohort autosomes.{pgen,pvar,psam} (from 08_pheno_and_pfiles.sh)
# Outputs: clean + pruned + unrelated datasets, PCA, keep/drop lists

set -euo pipefail

# --- Locate repo + load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

ensure_dirs
mkdir -p "$PLINK_DIR" "$QC_DIR"

# --- Parameters (can override via env)
GENO="${GENO:-0.02}"   # per-variant missingness
MIND="${MIND:-0.02}"   # per-sample missingness
MAF="${MAF:-0.01}"     # minor allele frequency
HWE="${HWE:-1e-6}"     # Hardy–Weinberg
LD_WIN="${LD_WIN:-200}" LD_STEP="${LD_STEP:-50}" LD_R2="${LD_R2:-0.2}"
NCOMP="${NCOMP:-20}"
KING_CUTOFF="${KING_CUTOFF:-0.125}"

# --- Input prefix (from 08 step)
PFILE_IN="$PLINK_DIR/plink_tmp/autosomes"
for e in pgen pvar psam; do
  [[ -s "$PFILE_IN.$e" ]] || { echo "[ERROR] Missing $PFILE_IN.$e"; exit 1; }
done

# --- Output bases
CLEAN_BASE="$PLINK_DIR/analysis.clean"
UNIQ_BASE="$CLEAN_BASE.uniq"
PRUNE_BASE="$UNIQ_BASE"
PRUNED_BASE="$UNIQ_BASE.pruned"
PCA_BASE="$UNIQ_BASE.pca"
UNREL_BASE="$UNIQ_BASE.unrel"
UNREL_PCA_BASE="$UNREL_BASE.pca"
KEEP_LIST="$UNIQ_BASE.king.keep"
DROP_LIST="$UNIQ_BASE.king.drop"

echo "== QC FILTERS on $PFILE_IN =="

# --- Step 1: filtering
plink2 --pfile "$PFILE_IN" --threads "${THREADS:-16}" \
  --geno "$GENO" --mind "$MIND" --maf "$MAF" --hwe "$HWE" \
  --make-pgen --out "$CLEAN_BASE"

# --- Step 2: unique IDs
plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --set-missing-var-ids @:#:\$r:\$a \
  --rm-dup exclude-all \
  --make-pgen --out "$UNIQ_BASE"

# --- Step 3: LD prune
plink2 --pfile "$PRUNE_BASE" --threads "${THREADS:-16}" \
  --indep-pairwise "$LD_WIN" "$LD_STEP" "$LD_R2" --out "$PRUNE_BASE"

plink2 --pfile "$PRUNE_BASE" --threads "${THREADS:-16}" \
  --extract "$PRUNE_BASE.prune.in" --make-pgen --out "$PRUNED_BASE"

# --- Step 4: PCA on pruned set
plink2 --pfile "$PRUNE_BASE" --threads "${THREADS:-16}" \
  --extract "$PRUNE_BASE.prune.in" --pca "$NCOMP" --out "$PCA_BASE"

# --- Step 5: KING relatedness
plink2 --pfile "$UNIQ_BASE" --threads "${THREADS:-16}" \
  --king-cutoff "$KING_CUTOFF" --make-pgen --out "$UNREL_BASE"

awk 'NR>1{print $2}' "$UNIQ_BASE.psam" | sort > "$UNIQ_BASE.iid"
awk 'NR>1{print $2}' "$UNREL_BASE.psam" | sort > "$UNREL_BASE.iid"
comm -12 "$UNIQ_BASE.iid" "$UNREL_BASE.iid" > "$KEEP_LIST"
comm -23 "$UNIQ_BASE.iid" "$UNREL_BASE.iid" > "$DROP_LIST"
rm -f "$UNIQ_BASE.iid" "$UNREL_BASE.iid"

# --- Step 6: PCA on unrelated set
plink2 --pfile "$UNREL_BASE" --threads "${THREADS:-16}" \
  --extract "$PRUNE_BASE.prune.in" --pca "$NCOMP" --out "$UNREL_PCA_BASE"

echo "[DONE] QC filters complete."
echo "Outputs:"
echo "  Clean:      $CLEAN_BASE.{pgen,pvar,psam}"
echo "  Uniq:       $UNIQ_BASE.{pgen,pvar,psam}"
echo "  Pruned:     $PRUNED_BASE.{pgen,pvar,psam}"
echo "  PCA (all):  $PCA_BASE.eigenvec / .eigenval"
echo "  Unrelated:  $UNREL_BASE.{pgen,pvar,psam}"
echo "  PCA (unrel):$UNREL_PCA_BASE.eigenvec / .eigenval"
echo "  Keep list:  $KEEP_LIST"
echo "  Drop list:  $DROP_LIST"

# --------------- RUN ---------------
# source scripts/00_config.sh
# bash scripts/09_qc_filters.sh