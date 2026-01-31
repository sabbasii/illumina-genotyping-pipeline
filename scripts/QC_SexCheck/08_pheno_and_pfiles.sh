#!/usr/bin/env bash
set -euo pipefail

# 08_pheno_and_pfiles.sh — Build cohort pfiles, then attach phenotype/sex in a second pass.
# Reproducible + robust:
#   • Do NOT filter at import (VCF is already cohort-only)
#   • Build numeric phenotype (PHENO1 ∈ {1,2,NA}) and attach via --pheno
#   • Apply sex via --update-sex
#   • Defensive checks for outputs at each step

# --- Locate repo root & load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- Tools
command -v plink2 >/dev/null 2>&1 || { echo "[ERROR] plink2 not found on PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found on PATH"; exit 1; }

# --- Dirs
ensure_dirs
mkdir -p "$TMP_DIR" "$PHENO_DIR" "$META_DIR" "$QC_SEXCHECK_DIR" "$QC_SUMMARIES_DIR"

# --- Resolve the exact VCF path (written by 07_qc_vcf_core.sh)
if [[ -s "$META_DIR/current_vcf.path" ]]; then
  FIXED_VCF="$(cat "$META_DIR/current_vcf.path")"
else
  FIXED_VCF="$VCF_DIR/cohort.gtc.$REF_BUILD.norm.vcf.gz"
fi
[[ -s "$FIXED_VCF" ]] || { echo "[ERROR] VCF not found: $FIXED_VCF"; exit 1; }
[[ -s "$META_DIR/vcf.samples" ]] || { echo "[ERROR] Canonical samples missing: $META_DIR/vcf.samples (run 07_qc_vcf_core.sh)."; exit 1; }

echo "== Phenotype + pfiles (decoupled) =="
echo "[VCF] $FIXED_VCF"
echo "[CANONICAL SAMPLES] $META_DIR/vcf.samples"
echo

# --- 1) Refresh expression overlap + keep list (idempotent; fast)
if [[ -f "$REPO_ROOT/scripts/03_expr_prepare_metadata.py" ]]; then python3 "$REPO_ROOT/scripts/03_expr_prepare_metadata.py"; fi
if [[ -f "$REPO_ROOT/scripts/04_expr_overlap_and_select.py" ]]; then python3 "$REPO_ROOT/scripts/04_expr_overlap_and_select.py"; fi
if [[ -f "$REPO_ROOT/scripts/05_make_cohort_id_lists.py" ]]; then python3 "$REPO_ROOT/scripts/05_make_cohort_id_lists.py"; fi

# Keep list is still useful for diagnostics; we do not filter at import anymore.
[[ -s "$PHENO_DIR/iid_selected.keep" ]] || { echo "[ERROR] Missing keep list: $PHENO_DIR/iid_selected.keep"; exit 1; }

PHENO_PSAM="$PHENO_DIR/cohort.pheno.psam"       # reference copy
PHENO_UPDATE="$PHENO_DIR/pheno_update.txt"      # FID IID PHENO1 (headered)
SEX_UPDATE="$PHENO_DIR/sex_update.txt"          # IID SEX (no header)

# --- 2) Build strict PSAM + emit update files (numeric PHENO1; sanitize whitespace/CRs)
python3 - "$META_DIR" "$EXPR_OUT_DIR" "$PHENO_DIR" "$SAMPLE_SHEET" <<'PY'
import os, sys, pandas as pd
meta_dir, expr_dir, pheno_dir, sample_sheet = sys.argv[1:]
vcf_samples = os.path.join(meta_dir, "vcf.samples")
u2i_path    = os.path.join(pheno_dir, "uasg_to_iid.tsv")    # from 05_make_cohort_id_lists.py
keep_path   = os.path.join(pheno_dir, "iid_selected.keep")
meta_path   = os.path.join(expr_dir, "meta_selected.csv")
psam_out    = os.path.join(pheno_dir, "cohort.pheno.psam")
pheno_upd   = os.path.join(pheno_dir, "pheno_update.txt")   # FID IID PHENO1 (headered)
sex_upd     = os.path.join(pheno_dir, "sex_update.txt")     # IID SEX (no header)

# Load IIDs from VCF (canonical)
vcf = pd.read_csv(vcf_samples, header=None, names=["IID"], dtype=str)
vcf["#FID"] = vcf["IID"]

# IID <-> UASG map
if os.path.isfile(u2i_path):
    u2i = pd.read_csv(u2i_path, sep="\t", names=["UASG","IID"], dtype=str)
    i2u = u2i.drop_duplicates("IID")[["IID","UASG"]]
else:
    ss0 = pd.read_csv(sample_sheet, sep=None, engine="python", skiprows=8, dtype=str)
    ss0.columns = ss0.columns.str.strip()
    i2u = pd.DataFrame({
        "IID":  ss0["SentrixBarcode_A"].str.strip() + "_" + ss0["SentrixPosition_A"].str.strip(),
        "UASG": ss0["Sample_Name"].str.strip()
    })

# Metadata and phenotype mapping
meta = pd.read_csv(meta_path, dtype=str).fillna("")
meta.columns = meta.columns.str.strip()
if "Final Diagnosis" in meta.columns:
    meta = meta.rename(columns={"Final Diagnosis":"Final_Diagnosis"})
d = meta["Final_Diagnosis"].astype(str).str.strip().str.lower()
meta["StrokeStatus"] = ""
meta.loc[d.eq("control"), "StrokeStatus"] = "Control"
meta.loc[d.isin(["ischemic stroke","tia"]), "StrokeStatus"] = "Case"

# PHENO1 numeric-only with NA for missing
meta["PHENO1"] = meta["StrokeStatus"].map({"Control": "1", "Case": "2"}).fillna("NA")

# SEX from SampleSheet (1=male, 2=female)
ss = pd.read_csv(sample_sheet, sep=None, engine="python", skiprows=8, dtype=str)
ss.columns = ss.columns.str.strip()
ss["IID"] = ss["SentrixBarcode_A"].astype(str).str.strip() + "_" + ss["SentrixPosition_A"].astype(str).str.strip()
gcol = next((c for c in ["Gender","Sex","gender","sex"] if c in ss.columns), None)
if gcol:
    g = ss[gcol].astype(str).str.strip().str.lower()
    ss["SEX"] = g.map({"m":"1","male":"1","f":"2","female":"2"}).fillna("")
else:
    ss["SEX"] = ""
sexmap = ss[["IID","SEX"]].drop_duplicates()

# Assemble PSAM on kept cohort (intersection with canonical VCF samples)
keep = pd.read_csv(keep_path, header=None, names=["IID"], dtype=str)
ph = (vcf.merge(keep, on="IID", how="inner")[["#FID","IID"]]
        .merge(sexmap, on="IID", how="left")
        .merge(i2u,    on="IID", how="left")
        .merge(meta[["UASG","StrokeStatus","Final_Diagnosis","PHENO1"]], on="UASG", how="left")).fillna("")

# Sanitize fields (strip CRs/quotes/space)
for c in ["#FID","IID","SEX","UASG","StrokeStatus","Final_Diagnosis","PHENO1"]:
    ph[c] = (ph[c].astype(str)
                  .str.replace('\r','', regex=False)
                  .str.replace('"','',  regex=False)
                  .str.strip())

# Enforce PHENO1 tokens
ph["PHENO1"] = ph["PHENO1"].where(ph["PHENO1"].isin(["1","2","NA"]), "NA")

# Write PSAM (reference)
psam_cols = ["#FID","IID","SEX","UASG","StrokeStatus","Final_Diagnosis","PHENO1"]
ph[psam_cols].to_csv(psam_out, sep="\t", index=False, lineterminator="\n")

# Emit update files for PLINK second-pass
# pheno_update: headered FID IID PHENO1
upd = ph[["#FID","IID","PHENO1"]].rename(columns={"#FID":"FID"})
upd.to_csv(pheno_upd, sep="\t", index=False, lineterminator="\n")

# sex_update: IID SEX (no header)
sex = ph[["IID","SEX"]]
sex.to_csv(sex_upd, sep="\t", index=False, header=False, lineterminator="\n")

print(f"[BUILD] PSAM={psam_out} rows={len(ph)} | pheno_update={pheno_upd} | sex_update={sex_upd}")
PY

# Validate PHENO1 values in update file (must be only 1/2/NA)
awk -F'\t' 'NR>1{ if ($3!="1" && $3!="2" && $3!="NA") { bad=1; print "[ERR] Bad PHENO1 in pheno_update at line " NR ": \"" $3 "\"" } } END{ if (bad) exit 1 }' "$PHENO_UPDATE"

# --- 3) Import autosomes from VCF — NO PSAM and NO --keep (VCF is already cohort-only)
plink2 --vcf "$FIXED_VCF" \
  --double-id --allow-extra-chr --chr-set "${CHRSET_AUTOSOMES:-37}" --chr 1-22 \
  --threads "${THREADS:-16}" \
  --make-pgen --out "$TMP_DIR/autosomes_raw"

# Guard: ensure import completed
[[ -s "$TMP_DIR/autosomes_raw.psam" ]] || { echo "[ERROR] Import did not produce $TMP_DIR/autosomes_raw.psam"; exit 1; }

# --- 4) Attach phenotype/sex and rewrite final autosomes PGEN
plink2 --pfile "$TMP_DIR/autosomes_raw" \
  --update-sex "$SEX_UPDATE" \
  --pheno "$PHENO_UPDATE" --pheno-name PHENO1 \
  --no-categorical \
  --threads "${THREADS:-16}" \
  --make-pgen --out "$TMP_DIR/autosomes"

# --- 5) Summaries on final autosomes
base="$QC_SUMMARIES_DIR/cohort"
plink2 --pfile "$TMP_DIR/autosomes" \
  --threads "${THREADS:-16}" \
  --freq --missing --hardy --out "$base" || true
gzip -f "$base".{afreq,smiss,vmiss,hardy} 2>/dev/null || true

# --- 6) chrX pfiles (import WITH sex) + sex-check ---------------------
# Use PSAM (has SEX) at import time; otherwise PLINK refuses chrX.
PHENO_PSAM="$PHENO_DIR/cohort.pheno.psam"

# Guard: PSAM must exist and have a SEX column
[[ -s "$PHENO_PSAM" ]] || { echo "[ERROR] Missing PSAM: $PHENO_PSAM"; exit 1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="SEX") found=1} END{if(!found) exit 1}' "$PHENO_PSAM" \
  || { echo "[ERROR] PSAM lacks SEX column: $PHENO_PSAM"; exit 1; }

# Import chrX with sex provided via --psam
plink2 --vcf "$FIXED_VCF" \
  --psam "$PHENO_PSAM" \
  --double-id \
  --allow-extra-chr --chr-set "${CHRSET_AUTOSOMES:-37}" \
  --split-par "${SPLIT_PAR:-b37}" \
  --threads "${THREADS:-16}" \
  --chr X \
  --make-pgen --out "$QC_SEXCHECK_DIR/chrX"

# Guard: ensure chrX was produced
[[ -s "$QC_SEXCHECK_DIR/chrX.psam" ]] || { 
  echo "[ERROR] chrX import failed; $QC_SEXCHECK_DIR/chrX.psam not found"; 
  echo "Check: $QC_SEXCHECK_DIR/chrX.log"; exit 1; 
}

# Optional summaries
plink2 --pfile "$QC_SEXCHECK_DIR/chrX" \
  --threads "${THREADS:-16}" \
  --freq --missing --hardy --out "$QC_SUMMARIES_DIR/cohort.chrX" || true
gzip -f "$QC_SUMMARIES_DIR/cohort.chrX."{afreq,smiss,vmiss,hardy.x} 2>/dev/null || true

# Sex-check (sex already present via PSAM)
plink2 --pfile "$QC_SEXCHECK_DIR/chrX" \
  --threads "${THREADS:-16}" \
  --check-sex 'min-male-xf=0.8' 'max-female-yrate=0.02' \
  --out "$QC_SEXCHECK_DIR/cohort.sexcheck"
