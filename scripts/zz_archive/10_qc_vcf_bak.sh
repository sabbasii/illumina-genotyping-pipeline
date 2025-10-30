#!/usr/bin/env bash
set -euo pipefail

# 10_qc_vcf.sh — reproducible QC (fix header contigs; summaries; PLINK2 autosomes+chrX)
#
# Usage:
#   source scripts/00_config.sh
#   bash scripts/10_qc_vcf.sh                   # default: $VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz
#   bash scripts/10_qc_vcf.sh path/to/file.vcf.gz
#   MAKE_PRIMARY=1 bash scripts/10_qc_vcf.sh    # (placeholder; not used here)

# --- Locate repo root and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- Tools
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools not found on PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[QC-WARN] samtools not on PATH; will skip FASTA indexing if needed."; }
command -v plot-vcfstats >/dev/null 2>&1 || { echo "[QC-WARN] plot-vcfstats not found; plots will be skipped."; }

# --- Inputs & outputs
VCF_INPUT="${1:-$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz}"
FASTA="$REFERENCE_FASTA"
ensure_dirs
mkdir -p "$QC_SUMMARIES_DIR" "$QC_REPORTS_DIR" "$QC_SEXCHECK_DIR" "$PLINK_DIR" "$TMP_DIR"

echo "[QC] File: $VCF_INPUT"
echo "[QC] Summaries:  $QC_SUMMARIES_DIR"
echo "[QC] Reports:    $QC_REPORTS_DIR"
echo "[QC] Sexcheck:   $QC_SEXCHECK_DIR"
echo "[QC] Ref:        $FASTA"

# --- Ensure FASTA index exists (soft)
if [[ ! -f "${FASTA}.fai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    echo "[QC] Building FASTA index: ${FASTA}.fai"
    samtools faidx "$FASTA"
  else
    echo "[QC-WARN] samtools missing; cannot create ${FASTA}.fai automatically."
  fi
fi

# --- Ensure index exists on VCF/BCF so -n and random access work
if [[ "$VCF_INPUT" =~ \.vcf\.gz$ ]]; then
  [[ -f "${VCF_INPUT}.tbi" ]] || bcftools index -t "$VCF_INPUT"
else
  [[ -f "${VCF_INPUT}.csi" ]] || bcftools index "$VCF_INPUT"
fi

# --- Check if header has contigs; if not, inject from FASTA index
HAS_CONTIGS=0
bcftools view -h "$VCF_INPUT" | grep -q '^##contig=<ID=' && HAS_CONTIGS=1

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

# --- 1) Header & contigs
bcftools view -h "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/header.txt" >/dev/null
bcftools view -h "$FIXED_VCF" \
  | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
  | sort -V \
  | tee "$QC_SUMMARIES_DIR/contigs.list" >/dev/null

# --- 2) Samples & count
bcftools query -l "$FIXED_VCF" | tee "$QC_SUMMARIES_DIR/samples.list" >/dev/null
bcftools query -l "$FIXED_VCF" | wc -l | tee "$QC_SUMMARIES_DIR/samples.count"

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

# --- 6) REF/ALT sanity vs reference (no modifications to file)
bcftools norm -f "$FASTA" -c ws -Ou "$FIXED_VCF" >/dev/null 2>"$QC_SUMMARIES_DIR/refcheck.log" || true
grep -E "REF_MISMATCH|duplicate|REF_N|Bad REF" "$QC_SUMMARIES_DIR/refcheck.log" \
  | tee "$QC_SUMMARIES_DIR/ref_mismatch.summary" || true

# --- 7) AF snapshot (on-the-fly so it works even if AF absent in source)
echo "[QC] Generating AF snapshot via +fill-tags ..."
if bcftools +fill-tags "$FIXED_VCF" -Ou -- -t AF \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AF\n' \
    | gzip -c > "$QC_SUMMARIES_DIR/af.tsv.gz"; then
  echo "[QC] AF snapshot written to $QC_SUMMARIES_DIR/af.tsv.gz"
else
  echo "[QC-WARN] Could not generate AF snapshot." >&2
fi

# --- 8) Prepare PSAM/sexmap (align IIDs to VCF sample IDs; idempotent)
bcftools query -l "$FIXED_VCF" > "$META_DIR/vcf.samples"
bash "$REPO_ROOT/scripts/11_build_psam_from_barcode.sh" \
  "$SHEET_PATH" \
  "$META_DIR/vcf.samples" \
  "$PSAM_SEX" \
  "$SEXMAP_PATH"

# --- 9) PLINK2 QC (autosomes stats to summaries; autosomes pfiles kept in TMP_DIR)
if command -v plink2 >/dev/null 2>&1; then
  base="$QC_SUMMARIES_DIR/cohort"
  mkdir -p "$TMP_DIR"

  # autosomes: keep pfiles in TMP_DIR so downstream filters can start from there
  plink2 --vcf "$FIXED_VCF" \
    --double-id --allow-extra-chr --chr-set "$CHRSET_AUTOSOMES" --chr 1-22 \
    --threads "${THREADS:-16}" \
    --make-pgen --out "$TMP_DIR/autosomes" || true

  plink2 --pfile "$TMP_DIR/autosomes" \
    --threads "${THREADS:-16}" \
    --freq --missing --hardy --out "$base" || true

  gzip -f "$base".{afreq,smiss,vmiss,hardy} 2>/dev/null || true

  # chrX: pfiles + sexcheck under QC_SEXCHECK_DIR; chrX summaries into summaries
  if [[ -s "$PSAM_SEX" ]]; then
    echo "[QC] chrX: using PSAM $PSAM_SEX and splitting PAR ($SPLIT_PAR)"
    chrxp="$QC_SEXCHECK_DIR/chrX"
    mkdir -p "$QC_SEXCHECK_DIR"

    plink2 --vcf "$FIXED_VCF" \
      --psam "$PSAM_SEX" \
      --split-par "$SPLIT_PAR" \
      --double-id \
      --threads "${THREADS:-16}" \
      --chr X \
      --make-pgen --out "$chrxp" || true

    plink2 --pfile "$chrxp" \
      --threads "${THREADS:-16}" \
      --freq --missing --hardy --out "$QC_SUMMARIES_DIR/cohort.chrX" || true
    gzip -f "$QC_SUMMARIES_DIR/cohort.chrX."{afreq,smiss,vmiss,hardy.x} 2>/dev/null || true

    plink2 --pfile "$chrxp" \
      --threads "${THREADS:-16}" \
      --check-sex 'min-male-xf=0.8' 'max-female-yrate=0.02' \
      --out "$QC_SEXCHECK_DIR/cohort.sexcheck" || true
  else
    echo "[QC-INFO] No PSAM at $PSAM_SEX; skipping chrX and sex-check."
  fi

  # NOTE: Do NOT delete $TMP_DIR; autosomes.* are needed for the next stage.
else
  echo "[QC-INFO] plink2 not found; skipping PLINK QC."
fi
