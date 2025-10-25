#!/usr/bin/env bash
# 20_qc_filters.sh ??? filters (GENO/MIND/MAF/HWE) ??? LD-prune ??? PCA ??? relatedness (KING)
# Defaults via env: GENO=0.02 MIND=0.02 MAF=0.01 HWE=1e-6 LD_WIN=200 LD_STEP=50 LD_R2=0.2 NCOMP=20 KING_CUTOFF=0.125
# Usage:
#   source scripts/00_config.sh
#   bash scripts/20_qc_filters.sh [--pfile-prefix /path/to/prefix]

# --- safe exit + strict-mode guard (prevents closing your shell if sourced)
_exit(){ local code=${1:-0}; if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then exit "$code"; else return "$code"; fi }
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then set -euo pipefail; fi
export LC_ALL=C

# -----------------------------
# Locate repo root & load config
# -----------------------------
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

ensure_dirs
mkdir -p "$PLINK_DIR" "$TMP_DIR" "$QC_SUMMARIES_DIR" "$QC_DIR"

SUMMARY_LOG="$OUT_DIR/qc/qc_filters_summary.txt"

# -----------------------------
# Parameters (env-overridable)
# -----------------------------
GENO="${GENO:-0.02}"
MIND="${MIND:-0.02}"
MAF="${MAF:-0.01}"
HWE="${HWE:-1e-6}"
LD_WIN="${LD_WIN:-200}"
LD_STEP="${LD_STEP:-50}"
LD_R2="${LD_R2:-0.2}"
NCOMP="${NCOMP:-20}"
KING_CUTOFF="${KING_CUTOFF:-0.125}"

# -----------------------------
# CLI
# -----------------------------
PFILE_IN=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pfile-prefix) PFILE_IN="$2"; shift 2;;
    -h|--help)
      echo "Usage: bash scripts/20_qc_filters.sh [--pfile-prefix /path/to/prefix]"
      echo "Env: GENO MIND MAF HWE LD_WIN LD_STEP LD_R2 NCOMP KING_CUTOFF"
      _exit 0;;
    *) echo "[ERROR] Unknown arg: $1"; _exit 2;;
  esac
done
[[ -n "$PFILE_IN" ]] || PFILE_IN="$TMP_DIR/autosomes"

# -----------------------------
# Preconditions / tools
# -----------------------------
command -v plink2 >/dev/null 2>&1 || { echo "[ERROR] plink2 not found"; _exit 1; }

need_pfile() { for e in pgen pvar psam; do [[ -s "$1.$e" ]] || { echo "[ERROR] Missing $1.$e"; return 1; }; done; }
need_pfile "$PFILE_IN" || _exit 1

# -----------------------------
# Helpers
# -----------------------------
ts() { date +"%Y%m%d_%H%M%S"; }
say() { echo "[$(date +%H:%M:%S)] $*"; }
count_samples() { awk 'NR>1{c++} END{print (c+0)}' "$1"; }
count_variants() { awk 'NR>1{c++} END{print (c+0)}' "$1"; }

# -----------------------------
# Output prefixes (no overwrite)
# -----------------------------
CLEAN_BASE="$PLINK_DIR/analysis.clean"
if [[ -e "${CLEAN_BASE}.pgen" || -e "${CLEAN_BASE}.pvar" || -e "${CLEAN_BASE}.psam" ]]; then
  CLEAN_BASE="${CLEAN_BASE}_$(ts)"
  say "[INFO] analysis.clean exists; using timestamped prefix: $(basename "$CLEAN_BASE")"
fi
PRUNE_BASE="$CLEAN_BASE"
PRUNED_BASE="${CLEAN_BASE}.pruned"
PCA_BASE="${CLEAN_BASE}.pca"
UNREL_BASE="${CLEAN_BASE}.unrel"
UNREL_PCA_BASE="${UNREL_BASE}.pca"
KEEP_LIST="${CLEAN_BASE}.king.cutoff${KING_CUTOFF}.keep"
DROP_LIST="${CLEAN_BASE}.king.cutoff${KING_CUTOFF}.removelist"

# -----------------------------
# Snapshot counts (input)
# -----------------------------
S0S=$(count_samples "${PFILE_IN}.psam")
S0V=$(count_variants "${PFILE_IN}.pvar")
say "[START] Input pfile: $PFILE_IN  samples=$S0S variants=$S0V"

# -----------------------------
# Step 1: Summaries on input (reference only)
# -----------------------------
say "[STEP] Summarizing input (missingness/freq/HWE)"
plink2 --pfile "$PFILE_IN" --threads "${THREADS:-16}" \
  --freq --missing --hardy --out "$QC_SUMMARIES_DIR/input.autosomes" >/dev/null 2>&1 || true
gzip -f "$QC_SUMMARIES_DIR"/input.autosomes.{afreq,smiss,vmiss,hardy} 2>/dev/null || true

# -----------------------------
# Step 2: Filters ??? analysis.clean
# -----------------------------
say "[STEP] Filtering: --geno $GENO --mind $MIND --maf $MAF --hwe $HWE"
plink2 --pfile "$PFILE_IN" --threads "${THREADS:-16}" \
  --geno "$GENO" --mind "$MIND" --maf "$MAF" --hwe "$HWE" \
  --make-pgen --out "$CLEAN_BASE"

plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --freq --missing --hardy --out "$CLEAN_BASE" >/dev/null 2>&1 || true
gzip -f "$CLEAN_BASE".{afreq,smiss,vmiss,hardy} 2>/dev/null || true

S1S=$(count_samples "${CLEAN_BASE}.psam")
S1V=$(count_variants "${CLEAN_BASE}.pvar")
say "[INFO] Clean set: samples=$S1S variants=$S1V -> ${CLEAN_BASE}.{pgen,pvar,psam}"

# -----------------------------
# Step 3: LD prune ??? prune.in/out + pruned pfile
# -----------------------------
say "[STEP] LD-prune: --indep-pairwise $LD_WIN $LD_STEP $LD_R2"
plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --indep-pairwise "$LD_WIN" "$LD_STEP" "$LD_R2" --out "$PRUNE_BASE"

plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --extract "${PRUNE_BASE}.prune.in" --make-pgen --out "$PRUNED_BASE"

S2V=$(count_variants "${PRUNED_BASE}.pvar" || echo 0)
say "[INFO] Pruned set: variants=$S2V -> ${PRUNED_BASE}.{pgen,pvar,psam}"

# -----------------------------
# Step 4: PCA on pruned set
# -----------------------------
say "[STEP] PCA on pruned set: NCOMP=$NCOMP"
plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --extract "${PRUNE_BASE}.prune.in" --pca "$NCOMP" approx --out "$PCA_BASE"

# -----------------------------
# Step 5: KING relatedness ??? unrelated set
# -----------------------------
say "[STEP] KING cutoff: $KING_CUTOFF"
plink2 --pfile "$CLEAN_BASE" --threads "${THREADS:-16}" \
  --king-cutoff "$KING_CUTOFF" --make-pgen --out "$UNREL_BASE"

awk 'NR>1{print $2}' "${CLEAN_BASE}.psam" | sort > "${CLEAN_BASE}.iid"
awk 'NR>1{print $2}' "${UNREL_BASE}.psam" | sort > "${UNREL_BASE}.iid"
comm -12 "${CLEAN_BASE}.iid" "${UNREL_BASE}.iid" > "$KEEP_LIST"
comm -23 "${CLEAN_BASE}.iid" "${UNREL_BASE}.iid" > "$DROP_LIST"
rm -f "${CLEAN_BASE}.iid" "${UNREL_BASE}.iid"

S3S=$(count_samples "${UNREL_BASE}.psam")
say "[INFO] Unrelated set: samples=$S3S"
say "[INFO] Keep list : $KEEP_LIST  ($(wc -l < "$KEEP_LIST" | tr -d ' ') IIDs)"
say "[INFO] Remove list: $DROP_LIST ($(wc -l < "$DROP_LIST" | tr -d ' ') IIDs)"

# -----------------------------
# Step 6: PCA on unrelateds (same pruned SNPs)
# -----------------------------
say "[STEP] PCA on unrelateds"
plink2 --pfile "$UNREL_BASE" --threads "${THREADS:-16}" \
  --extract "${PRUNE_BASE}.prune.in" --pca "$NCOMP" approx --out "$UNREL_PCA_BASE"

# -----------------------------
# Step 7: Summary log
# -----------------------------
{
  echo "========================================"
  echo "QC FILTERS ($RUN)  $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "Input pfile : $PFILE_IN"
  echo "Params      : GENO=$GENO  MIND=$MIND  MAF=$MAF  HWE=$HWE"
  echo "LD prune    : window=$LD_WIN step=$LD_STEP r2=$LD_R2"
  echo "PCA         : NCOMP=$NCOMP"
  echo "KING cutoff : $KING_CUTOFF"
  echo "Counts      :"
  printf "  before filtering:   %8d samples   %10d variants\n" "$S0S" "$S0V"
  printf "  after filtering:    %8d samples   %10d variants\n" "$S1S" "$S1V"
  printf "  after LD-prune:                  %10d variants (pruned set)\n" "$S2V"
  printf "  unrelated (KING):    %8d samples\n" "$S3S"
  echo "Outputs:"
  echo "  CLEAN       : ${CLEAN_BASE}.{pgen,pvar,psam}"
  echo "  CLEAN QC    : ${CLEAN_BASE}.{afreq.gz,smiss.gz,vmiss.gz,hardy.gz}"
  echo "  PRUNE files : ${PRUNE_BASE}.prune.in / .prune.out"
  echo "  PRUNED PFILE: ${PRUNED_BASE}.{pgen,pvar,psam}"
  echo "  PCA (clean) : ${PCA_BASE}.eigenvec / .eigenval"
  echo "  UNRELATED   : ${UNREL_BASE}.{pgen,pvar,psam}"
  echo "  PCA (unrel) : ${UNREL_PCA_BASE}.eigenvec / .eigenval"
  echo "  KEEP/DROP   : ${KEEP_LIST} / ${DROP_LIST}"
} >> "$SUMMARY_LOG"

say "[DONE] QC filters complete. Summary appended to: $SUMMARY_LOG"
