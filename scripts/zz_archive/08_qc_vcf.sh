#!/usr/bin/env bash
set -euo pipefail

# 10_qc_vcf.sh — QC + bring phenotype earlier:
# - Fix header contigs; summaries; stats
# - Build barcode PSAM (#FID IID PAT MAT SEX PHENOTYPE)
# - Run expression overlap (01/02 python) → ids_selected.txt, meta_selected.csv
# - Build IID<->UASG map + iid_selected.keep (IIDs in cohort)
# - Build strict cohort PHENO PSAM (#FID IID SEX UASG StrokeStatus Final_Diagnosis PHENO1)
# - Create cohort-only autosomes pfiles with --psam + --keep
# - chrX pfiles + sex-check on the same filtered cohort

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
command -v plink2   >/dev/null 2>&1 || { echo "[QC-WARN] plink2 not on PATH; PLINK steps will be skipped."; }

# --- Inputs & outputs
VCF_INPUT="${1:-$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz}"
FASTA="$REFERENCE_FASTA"
ensure_dirs
mkdir -p "$QC_SUMMARIES_DIR" "$QC_REPORTS_DIR" "$QC_SEXCHECK_DIR" "$PLINK_DIR" "$TMP_DIR" "$META_DIR" "$PHENO_DIR"

echo "[QC] File:       $VCF_INPUT"
echo "[QC] Summaries:  $QC_SUMMARIES_DIR"
echo "[QC] Reports:    $QC_REPORTS_DIR"
echo "[QC] Sexcheck:   $QC_SEXCHECK_DIR"
echo "[QC] Ref:        $FASTA"
echo

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

# --- 8) Prepare barcode PSAM + sexmap (align IIDs to VCF sample IDs)
bcftools query -l "$FIXED_VCF" > "$META_DIR/vcf.samples"
bash "$REPO_ROOT/scripts/11_build_psam_from_barcode.sh" \
  "$SHEET_PATH" \
  "$META_DIR/vcf.samples" \
  "$PSAM_SEX" \
  "$SEXMAP_PATH"

echo
echo "==== Cohort integration (expression + phenotype) BEFORE pgen creation ===="

# --- 8.5) Run expression overlap to get ids_selected.txt + meta_selected.csv
# (Idempotent; will overwrite outputs in $EXPR_OUT_DIR)
echo "[Cohort] Running expression overlap scripts ..."
python3 "$REPO_ROOT/scripts/01_load_and_prepare_microarray_expression_metadata.py"
python3 "$REPO_ROOT/scripts/02_overlap_and_select_cohorts.py"

# --- Build mapping: IID <-> UASG from SampleSheet; then build keep list (IIDs for selected UASGs)
echo "[Cohort] Building IID<->UASG map + cohort keep list ..."
awk -F',' 'NR>9{gsub(/\r$/,""); print $1"_"$2"\t"$3}' "$SAMPLE_SHEET" \
  | sed 's/[[:space:]]\+$//' \
  > "$META_DIR/snp_iid_uasg.tsv"

join -1 1 -2 2 \
  <(sort "$EXPR_OUT_DIR/ids_selected.txt") \
  <(sort -k2,2 "$META_DIR/snp_iid_uasg.tsv") \
  | awk '{print $2}' > "$PHENO_DIR/iid_selected.keep"

echo "[Cohort] iid_selected.keep rows: $(wc -l < "$PHENO_DIR/iid_selected.keep" | tr -d " ")"

# --- Build a strict cohort PHENO PSAM (#FID IID SEX UASG StrokeStatus Final_Diagnosis PHENO1)
# Uses: vcf.samples, snp_iid_uasg.tsv, meta_selected.csv
PHENO_STRICT="$PHENO_DIR/cohort.pheno.psam"
python3 - "$META_DIR" "$EXPR_OUT_DIR" "$PHENO_DIR" <<'PY'
import os, sys, pandas as pd
meta_dir, expr_dir, pheno_dir = sys.argv[1:]
vcf_samples = os.path.join(meta_dir, "vcf.samples")
map_path    = os.path.join(meta_dir, "snp_iid_uasg.tsv")
meta_path   = os.path.join(expr_dir, "meta_selected.csv")
out_path    = os.path.join(pheno_dir, "cohort.pheno.psam")

vcf = pd.read_csv(vcf_samples, header=None, names=["IID"], dtype=str)
vcf["#FID"] = vcf["IID"]

map_df = pd.read_csv(map_path, sep="\t", names=["IID","UASG"], dtype=str)
meta = pd.read_csv(meta_path, dtype=str).fillna("")
if "Final Diagnosis" in meta.columns:
    meta = meta.rename(columns={"Final Diagnosis":"Final_Diagnosis"})
meta["StrokeStatus"] = ""
d = meta["Final_Diagnosis"].str.strip().str.lower()
meta.loc[d.eq("control"), "StrokeStatus"] = "Control"
meta.loc[d.isin(["ischemic stroke","tia"]), "StrokeStatus"] = "Case"
meta["PHENO1"] = meta["StrokeStatus"].map({"Control":"1","Case":"2"}).fillna("")

ph = (vcf[["#FID","IID"]]
      .merge(map_df, on="IID", how="left")
      .merge(meta[["UASG","StrokeStatus","Final_Diagnosis","PHENO1"]], on="UASG", how="left"))
ph["SEX"] = ""  # unknown at this stage; PLINK treats blank as missing; fill later if desired
ph = ph[["#FID","IID","SEX","UASG","StrokeStatus","Final_Diagnosis","PHENO1"]].fillna("")

ph.to_csv(out_path, sep="\t", index=False, lineterminator="\n")
print(f"Wrote strict PHENO PSAM → {out_path}  rows={len(ph)}")
PY

# --- 9) PLINK2: Create cohort-only autosomes pfiles (use strict PSAM + keep list)
if command -v plink2 >/dev/null 2>&1; then
  echo "[PLINK] Building cohort-only autosomes pfiles with phenotype attached ..."
  mkdir -p "$TMP_DIR"

  plink2 --vcf "$FIXED_VCF" \
    --psam   "$PHENO_STRICT" \
    --keep   "$PHENO_DIR/iid_selected.keep" \
    --double-id --allow-extra-chr --chr-set "$CHRSET_AUTOSOMES" --chr 1-22 \
    --threads "${THREADS:-16}" \
    --make-pgen --out "$TMP_DIR/autosomes"

  # Summaries on the cohort autosomes set
  base="$QC_SUMMARIES_DIR/cohort"
  plink2 --pfile "$TMP_DIR/autosomes" \
    --threads "${THREADS:-16}" \
    --freq --missing --hardy --out "$base" || true
  gzip -f "$base".{afreq,smiss,vmiss,hardy} 2>/dev/null || true

  # chrX pfiles + sex-check on the SAME FILTERED COHORT
  if [[ -s "$PHENO_STRICT" ]]; then
    echo "[PLINK] chrX: cohort-restricted sex-check using strict PSAM; split PAR ($SPLIT_PAR)"
    chrxp="$QC_SEXCHECK_DIR/chrX"
    mkdir -p "$QC_SEXCHECK_DIR"

    plink2 --vcf "$FIXED_VCF" \
      --psam "$PHENO_STRICT" \
      --keep "$PHENO_DIR/iid_selected.keep" \
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
    echo "[QC-INFO] Strict PHENO PSAM not found; skipping chrX and sex-check."
  fi

  echo
  echo "[DONE] Cohort-restricted pfiles (pre-QC) ready at: $TMP_DIR/autosomes.{pgen,pvar,psam}"
else
  echo "[QC-INFO] plink2 not found; skipping PLINK steps."
fi
