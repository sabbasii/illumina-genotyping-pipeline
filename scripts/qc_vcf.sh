#!/usr/bin/env bash
set -euo pipefail

# qc_vcf.sh — reproducible QC (self-fixing header contigs + optional primary-only VCF)
#
# Usage:
#   source scripts/00_config.sh
#   bash scripts/qc_vcf.sh                   # default: $VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz
#   bash scripts/qc_vcf.sh path/to/file.vcf.gz
#   MAKE_PRIMARY=1 bash scripts/qc_vcf.sh    # also emits *.norm.primary.vcf.gz
#
# Outputs go to: $OUT_DIR/qc

# --- Locate repo root and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"

# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- Tools
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH"; exit 1; }
if ! command -v plot-vcfstats >/dev/null 2>&1; then
  echo "[QC-WARN] plot-vcfstats not found; plots will be skipped."
fi

# --- Inputs & outputs
VCF_INPUT="${1:-$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz}"
FASTA="$REFERENCE_FASTA"
mkdir -p "$OUT_DIR/qc" "$OUT_DIR/qc/plink" "$OUT_DIR/qc/vcfstats_plots"

echo "[QC] File: $VCF_INPUT"
echo "[QC] Out:  $OUT_DIR/qc"
echo "[QC] Ref:  $FASTA"

# --- Ensure FASTA index exists
if [[ ! -f "${FASTA}.fai" ]]; then
  echo "[QC] Building FASTA index: ${FASTA}.fai"
  samtools faidx "$FASTA"
fi

# --- Ensure index exists on VCF/BCF so -n and random access work
if [[ "$VCF_INPUT" =~ \.vcf\.gz$ ]]; then
  [[ -f "${VCF_INPUT}.tbi" ]] || bcftools index -t "$VCF_INPUT"
else
  [[ -f "${VCF_INPUT}.csi" ]] || bcftools index "$VCF_INPUT"
fi

# --- Check if header has contigs; if not, inject from FASTA index (header-only reheader)
HAS_CONTIGS=0
if bcftools view -h "$VCF_INPUT" | grep -q '^##contig=<ID='; then
  HAS_CONTIGS=1
fi

FIXED_VCF="$VCF_INPUT"
if [[ "$HAS_CONTIGS" -eq 0 ]]; then
  echo "[QC] No ##contig= lines found. Injecting contigs from ${FASTA}.fai ..."
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

# --- 1) Header & contigs (robust extraction)
bcftools view -h "$FIXED_VCF" | tee "$OUT_DIR/qc/header.txt" >/dev/null
bcftools view -h "$FIXED_VCF" \
  | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
  | sort -V \
  | tee "$OUT_DIR/qc/contigs.list" >/dev/null

# --- 2) Samples & count
bcftools query -l "$FIXED_VCF" | tee "$OUT_DIR/qc/samples.list" >/dev/null
bcftools query -l "$FIXED_VCF" | wc -l | tee "$OUT_DIR/qc/samples.count"

# --- 3) Variant count (from index)
bcftools index -n "$FIXED_VCF" | tee "$OUT_DIR/qc/variants.count"

# --- 4) Stats (+ plots if available)
bcftools stats -s - "$FIXED_VCF" | tee "$OUT_DIR/qc/bcftools.stats.txt" >/dev/null
if command -v plot-vcfstats >/dev/null 2>&1; then
  plot-vcfstats -p "$OUT_DIR/qc/vcfstats_plots" "$OUT_DIR/qc/bcftools.stats.txt" 2>/dev/null || true
fi

# --- 5) Variants by contig
bcftools view -H "$FIXED_VCF" | cut -f1 | sort | uniq -c | sort -k2V \
  | tee "$OUT_DIR/qc/variants_by_contig.txt"

# --- 6) REF/ALT sanity vs reference (no modifications to file)
bcftools norm -f "$FASTA" -c ws -Ou "$FIXED_VCF" >/dev/null 2>"$OUT_DIR/qc/refcheck.log" || true
grep -E "REF_MISMATCH|duplicate|REF_N|Bad REF" "$OUT_DIR/qc/refcheck.log" \
  | tee "$OUT_DIR/qc/ref_mismatch.summary" || true

# --- 7) AF snapshot (on-the-fly so it works even if AF absent in source)
echo "[QC] Generating AF snapshot via +fill-tags ..."
if bcftools +fill-tags "$FIXED_VCF" -Ou -- -t AF \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' \
    | gzip -c > "$OUT_DIR/qc/af.tsv.gz"; then
  echo "[QC] AF snapshot written to $OUT_DIR/qc/af.tsv.gz"
else
  echo "[QC-WARN] Could not generate AF snapshot." >&2
fi

# --- Prepare PSAM/sexmap (aligns IIDs to VCF sample IDs; idempotent)
mkdir -p "$REPO_ROOT/metadata"
bcftools query -l "$FIXED_VCF" > "$REPO_ROOT/metadata/vcf.samples"

bash "$REPO_ROOT/scripts/build_psam_from_barcode.sh" \
  "$SHEET_PATH" \
  "$REPO_ROOT/metadata/vcf.samples" \
  "$PSAM_SEX" \
  "$SEXMAP_PATH"

# --- 8) PLINK2 QC (2-step: make pgen -> stats)
if command -v plink2 >/dev/null 2>&1; then
  base="$OUT_DIR/qc/plink/cohort"
  tmp="$OUT_DIR/qc/plink_tmp"; mkdir -p "$tmp"

  # autosomes: use custom chr-set; keep allow-extra-chr if you need non-primary contigs
  plink2 --vcf "$FIXED_VCF" \
    --double-id --allow-extra-chr --chr-set "$CHRSET_AUTOSOMES" --chr 1-22 \
    --threads "${THREADS:-16}" \
    --make-pgen --out "$tmp/autosomes" || true

  plink2 --pfile "$tmp/autosomes" \
    --threads "${THREADS:-16}" \
    --freq --missing --hardy --out "$base" || true

  gzip -f "$base".{afreq,smiss,vmiss,hardy} 2>/dev/null || true

  # chrX: IMPORTANT — do NOT pass --chr-set/--allow-extra-chr with --split-par
  if [[ -s "$PSAM_SEX" ]]; then
    echo "[QC] chrX: using PSAM $PSAM_SEX and splitting PAR ($SPLIT_PAR)"
    plink2 --vcf "$FIXED_VCF" \
      --psam "$PSAM_SEX" \
      --split-par "$SPLIT_PAR" \
      --double-id \
      --threads "${THREADS:-16}" \
      --chr X \
      --make-pgen --out "$tmp/chrX" || true

    plink2 --pfile "$tmp/chrX" \
      --threads "${THREADS:-16}" \
      --freq --missing --hardy --out "$OUT_DIR/qc/plink/cohort.chrX" || true

    # chrX HWE lands in .hardy.x
    gzip -f "$OUT_DIR/qc/plink/cohort.chrX."{afreq,smiss,vmiss,hardy.x} 2>/dev/null || true

    # tuned sex-check (optional)
    plink2 --pfile "$tmp/chrX" \
      --threads "${THREADS:-16}" \
      --check-sex 'min-male-xf=0.8' 'max-female-yrate=0.02' \
      --out "$OUT_DIR/qc/plink/cohort.sexcheck" || true
  else
    echo "[QC-INFO] No PSAM at $PSAM_SEX; skipping chrX and sex-check."
  fi

  rm -rf "$tmp"
else
  echo "[QC-INFO] plink2 not found; skipping PLINK QC."
fi
