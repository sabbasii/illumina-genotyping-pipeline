#!/usr/bin/env bash
set -euo pipefail

# 07_qc_vcf_core.sh — One-time QC on cohort-only VCF
# - Validates/indexes VCF, injects contigs if needed
# - Produces summaries/stats/AF snapshot
# - **Writes canonical sample list** -> $META_DIR/vcf.samples
# - Records the exact VCF path used -> $META_DIR/current_vcf.path
#
# Usage:
#   source scripts/00_config.sh
#   bash scripts/07_qc_vcf_core.sh [optional:/path/to/cohort.gtc.$REF_BUILD.norm.vcf.gz]

# --- Locate repo root & load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- Tools
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[WARN] samtools not on PATH; FASTA .fai build skipped if missing."; }
command -v plot-vcfstats >/dev/null 2>&1 || { echo "[WARN] plot-vcfstats not on PATH; plots will be skipped."; }

# --- Inputs/outputs
VCF_INPUT="${1:-$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz}"
FASTA="$REFERENCE_FASTA"

ensure_dirs
mkdir -p "$QC_SUMMARIES_DIR" "$QC_REPORTS_DIR" "$QC_SEXCHECK_DIR" "$PLINK_DIR" "$TMP_DIR" "$META_DIR" "$PHENO_DIR"

echo "== QC on cohort VCF (core) =="
echo "[IN ] VCF         : $VCF_INPUT"
echo "[IN ] FASTA       : $FASTA"
echo "[OUT] Summaries   : $QC_SUMMARIES_DIR"
echo "[OUT] Reports     : $QC_REPORTS_DIR"
echo "[OUT] Canonical samples → $META_DIR/vcf.samples"
echo

[[ -s "$VCF_INPUT" ]] || { echo "[ERROR] Missing VCF: $VCF_INPUT"; exit 1; }

# --- Ensure FASTA index exists (soft)
if [[ ! -f "${FASTA}.fai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    echo "[QC] Building FASTA index: ${FASTA}.fai"
    samtools faidx "$FASTA" || true
  else
    echo "[WARN] samtools missing; cannot create ${FASTA}.fai automatically."
  fi
fi

# --- Ensure index exists on VCF/BCF (for -n, random access)
if [[ "$VCF_INPUT" =~ \.vcf\.gz$ ]]; then
  [[ -f "${VCF_INPUT}.tbi" ]] || bcftools index -t "$VCF_INPUT"
else
  [[ -f "${VCF_INPUT}.csi" ]] || bcftools index "$VCF_INPUT"
fi

# --- Inject contigs if header lacks them
HAS_CONTIGS=0
bcftools view -h "$VCF_INPUT" | grep -q '^##contig=<ID=' && HAS_CONTIGS=1

FIXED_VCF="$VCF_INPUT"
if [[ "$HAS_CONTIGS" -eq 0 ]]; then
  echo "[QC] No ##contig= lines found. Injecting from ${FASTA}.fai ..."
  if [[ "$VCF_INPUT" =~ \.vcf\.gz$ ]]; then
    FIXED_VCF="${VCF_INPUT%.vcf.gz}.withContigs.vcf.gz"
    bcftools reheader -f "${FASTA}.fai" -o "$FIXED_VCF" "$VCF_INPUT"
    bcftools index -t "$FIXED_VCF"
  else
    FIXED_VCF="${VCF_INPUT%.bcf}.withContigs.bcf"
    bcftools reheader -f "${FASTA}.fai" -o "$FIXED_VCF" "$VCF_INPUT"
    bcftools index "$FIXED_VCF"
  fi
  echo "[QC] Using: $FIXED_VCF"
fi

# --- 1) Header & contigs
bcftools view -h "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/header.txt" >/dev/null
bcftools view -h "$FIXED_VCF" \
  | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
  | sort -V \
  | tee "$QC_SUMMARIES_DIR/contigs.list" >/dev/null

# --- 2) Samples (write to summaries AND canonical location)
bcftools query -l "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/samples.list" > "$META_DIR/vcf.samples"
wc -l < "$META_DIR/vcf.samples" | tee "$QC_SUMMARIES_DIR/samples.count"

# --- 3) Variant count (from index)
bcftools index -n "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/variants.count"

# --- 4) Stats (+ plots if available)
bcftools stats -s - "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/bcftools.stats.txt" >/dev/null
if command -v plot-vcfstats >/dev/null 2>&1; then
  plot-vcfstats -p "$QC_REPORTS_DIR/vcfstats_plots" "$QC_SUMMARIES_DIR/bcftools.stats.txt" 2>/dev/null || true
fi

# --- 5) Variants by contig
bcftools view -H "$FIXED_VCF" | cut -f1 | sort | uniq -c | sort -k2V \
  | tee "$QC_SUMMARIES_DIR/variants_by_contig.txt"

# --- 6) REF/ALT sanity vs reference (no modifications)
bcftools norm -f "$FASTA" -c ws -Ou "$FIXED_VCF" >/dev/null 2>"$QC_SUMMARIES_DIR/refcheck.log" || true
grep -E "REF_MISMATCH|duplicate|REF_N|Bad REF" "$QC_SUMMARIES_DIR/refcheck.log" \
  | tee "$QC_SUMMARIES_DIR/ref_mismatch.summary" || true

# --- 7) AF snapshot
echo "[QC] Generating AF snapshot via +fill-tags ..."
if bcftools +fill-tags "$FIXED_VCF" -Ou -- -t AF \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' \
    | gzip -c > "$QC_SUMMARIES_DIR/af.tsv.gz"; then
  echo "[QC] AF snapshot written to $QC_SUMMARIES_DIR/af.tsv.gz"
else
  echo "[WARN] Could not generate AF snapshot." >&2
fi

# --- Record the exact VCF path used
echo "$FIXED_VCF" > "$META_DIR/current_vcf.path"
echo
echo "[DONE] Core QC complete."
echo "       Canonical samples: $META_DIR/vcf.samples"
echo "       VCF used recorded in: $META_DIR/current_vcf.path"
