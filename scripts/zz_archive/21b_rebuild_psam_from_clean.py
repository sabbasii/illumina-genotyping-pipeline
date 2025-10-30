#!/usr/bin/env python3
"""
21b_rebuild_psam_from_clean.py
Rebuild a strict PLINK2 PSAM with exactly 7 columns:
#FID  IID  SEX  UASG  StrokeStatus  Final_Diagnosis  PHENO1

- Reads base PSAM:   $PLINK_DIR/analysis.clean.psam
- Reads current PHENO PSAM (to grab the phenotype columns already computed)
- Left-joins by IID onto the base PSAM (so row order and tabs are clean)
- Writes to: $PHENO_PSAM  (overwrites)
"""

import os
import pandas as pd

REPO_ROOT    = os.environ.get("REPO_ROOT", os.getcwd())
PLINK_DIR    = os.environ.get("PLINK_DIR", os.path.join(REPO_ROOT, "output", "genotype_run1", "qc", "plink"))
PHENO_DIR    = os.environ.get("PHENO_DIR", os.path.join(REPO_ROOT, "metadata", "pheno"))
PHENO_PSAM   = os.environ.get("PHENO_PSAM", os.path.join(PHENO_DIR, "cohort.pheno.psam"))

base_psam_path  = os.path.join(PLINK_DIR, "analysis.clean.psam")
pheno_psam_path = PHENO_PSAM  # use the current one only as a source of columns

# Load base PSAM (guaranteed PLINK2-valid)
base = pd.read_csv(base_psam_path, sep="\t", dtype=str)
base.columns = base.columns.str.strip()

# Ensure base has required IDs
if "#FID" not in base.columns:
    if "FID" in base.columns:
        base = base.rename(columns={"FID": "#FID"})
    else:
        raise KeyError("Base PSAM missing #FID/FID")
if "IID" not in base.columns:
    raise KeyError("Base PSAM missing IID")
if "SEX" not in base.columns:
    base["SEX"] = ""

# Load current phenotype PSAM (source of phenotype cols)
ph = pd.read_csv(pheno_psam_path, sep="\t", dtype=str, keep_default_na=False)
ph.columns = ph.columns.str.strip()

# Normalize column names (spaces → underscores for safety)
ph = ph.rename(columns={"Final Diagnosis": "Final_Diagnosis"})

# Keep just what we need, keyed by IID
need = ["IID", "UASG", "StrokeStatus", "Final_Diagnosis", "PHENO1"]
for c in need:
    if c not in ph.columns:
        ph[c] = ""
ph = ph[need].copy()

# Left join onto base by IID to preserve row order and ensure all tabs
out = base[["#FID", "IID", "SEX"]].merge(ph, on="IID", how="left")

# Coerce missing to empty strings, not 'nan'
out = out.fillna("")

# Enforce final column order and write with strict tabs + LF
cols = ["#FID", "IID", "SEX", "UASG", "StrokeStatus", "Final_Diagnosis", "PHENO1"]
out = out[cols]
os.makedirs(PHENO_DIR, exist_ok=True)
out.to_csv(PHENO_PSAM, sep="\t", index=False, lineterminator="\n")

print("Rebuilt PHENO PSAM →", PHENO_PSAM)
print("Rows:", len(out), "| Cols:", len(out.columns))
