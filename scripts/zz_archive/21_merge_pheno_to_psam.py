#!/usr/bin/env python3
"""
21_merge_pheno_to_psam.py  —  deterministic, PLINK2-safe PSAM writer

- Reads: analysis.clean(.unrel).psam, SampleSheet, meta_selected.csv, ids_selected.txt
- Merges phenotype columns
- Writes a PSAM whose FIRST header column is '#FID' (PLINK2 requirement)
- Guarantees column order + identical field count on every line
- Also writes iid_selected.keep (IIDs with PHENO1 in {1,2})
"""

import os, argparse
import pandas as pd

# --- env (set by 00_config.sh) ---
REPO_ROOT   = os.environ.get("REPO_ROOT", os.getcwd())
PLINK_DIR   = os.environ.get("PLINK_DIR", os.path.join(REPO_ROOT, "output", "genotype_run1", "qc", "plink"))
SAMPLE_SHEET= os.environ.get("SAMPLE_SHEET")
EXPR_OUT_DIR= os.environ.get("EXPR_OUT_DIR", os.path.join(REPO_ROOT, "output", "genotype_run1", "expr", "explore"))
PHENO_DIR   = os.environ.get("PHENO_DIR", os.path.join(REPO_ROOT, "metadata", "pheno"))
PHENO_PSAM  = os.environ.get("PHENO_PSAM", os.path.join(PHENO_DIR, "cohort.pheno.psam"))

META_SELECTED = os.path.join(EXPR_OUT_DIR, "meta_selected.csv")
IDS_SELECTED  = os.path.join(EXPR_OUT_DIR, "ids_selected.txt")

os.makedirs(PHENO_DIR, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument("--psam", default=os.path.join(PLINK_DIR, "analysis.clean.psam"),
                    help="Input .psam (e.g., analysis.clean.psam or analysis.clean.unrel.psam)")
args = parser.parse_args()

# --- load inputs ---
psam = pd.read_csv(args.psam, sep="\t")
ss   = pd.read_csv(SAMPLE_SHEET, sep=None, engine="python", skiprows=8, encoding="utf-8")
meta = pd.read_csv(META_SELECTED)
with open(IDS_SELECTED, "r", encoding="utf-8") as f:
    ids_selected = [x.strip() for x in f if x.strip()]

# normalize headers
psam.columns = psam.columns.str.strip()
ss.columns   = ss.columns.str.strip()
meta.columns = meta.columns.str.strip()

# ensure IID present; synthesize FID if missing
if "IID" not in psam.columns:
    raise KeyError(".psam missing column: IID")
if "FID" not in psam.columns:
    psam.insert(0, "FID", psam["IID"].astype(str))

# SampleSheet essentials
for col in ("Sample_Name","SentrixBarcode_A","SentrixPosition_A"):
    if col not in ss.columns:
        raise KeyError(f"SampleSheet missing column: {col}")

# Build IID <-> UASG map
ss["IID"]  = (ss["SentrixBarcode_A"].astype(str).str.strip() + "_" +
              ss["SentrixPosition_A"].astype(str).str.strip())
ss["UASG"] = ss["Sample_Name"].astype(str).str.strip()
map_iid_uasg = ss[["IID","UASG"]].drop_duplicates()

# Restrict to selected UASGs
sel = set(ids_selected)
map_iid_uasg = map_iid_uasg[map_iid_uasg["UASG"].isin(sel)].copy()

# Derive StrokeStatus + PHENO1
if "UASG" not in meta.columns:              raise KeyError("meta_selected.csv missing 'UASG'")
if "Final Diagnosis" not in meta.columns:   raise KeyError("meta_selected.csv missing 'Final Diagnosis'")
diag = meta["Final Diagnosis"].astype("string").str.strip().str.lower()
meta = meta.assign(StrokeStatus=pd.Series(pd.NA, index=meta.index, dtype="string"))
meta.loc[diag.eq("control"), "StrokeStatus"] = "Control"
meta.loc[diag.isin(["ischemic stroke","tia"]), "StrokeStatus"] = "Case"
meta_keep = meta[meta["StrokeStatus"].isin(["Control","Case"])].copy()

pheno = (map_iid_uasg
         .merge(meta_keep[["UASG","Final Diagnosis","StrokeStatus"]], on="UASG", how="inner"))
pheno["PHENO1"] = pheno["StrokeStatus"].map({"Control":1, "Case":2}).astype("Int64")

# Merge into PSAM by IID
psam_merged = psam.merge(pheno[["IID","UASG","StrokeStatus","Final Diagnosis","PHENO1"]],
                         on="IID", how="left")

# Reorder columns; enforce PLINK2 header begins with '#FID'
cols_base = list(psam.columns)  # includes FID, IID, SEX, ...
# Replace first column name with '#FID' while keeping data in the same column
cols_final = [*cols_base, "UASG", "StrokeStatus", "Final Diagnosis", "PHENO1"]
psam_merged.columns = [("#FID" if c=="FID" and i==0 else c)
                       for i,c in enumerate(psam_merged.columns)]

# Ensure all final columns exist (create if missing)
for c in ("UASG","StrokeStatus","Final Diagnosis","PHENO1"):
    if c not in psam_merged.columns:
        psam_merged[c] = pd.NA

# Keep only the intended order
psam_out = psam_merged[[("#FID" if i==0 else c) for i,c in enumerate(cols_base)] + 
                       ["UASG","StrokeStatus","Final Diagnosis","PHENO1"]]

# Tab-delimited with consistent newlines; no stray CR
psam_out.to_csv(PHENO_PSAM, sep="\t", index=False, lineterminator="\n")

# Keep-list of phenotyped IIDs present in PSAM
iid_keep = psam_out.loc[psam_out["PHENO1"].isin([1,2]), "IID"].dropna().astype(str)
iid_keep_path = os.path.join(PHENO_DIR, "iid_selected.keep")
iid_keep.to_csv(iid_keep_path, index=False, header=False)

print("=== Merge PHENO → PSAM ===")
print(f"Input PSAM     : {args.psam}")
print(f"PSAM rows      : {len(psam_out)}")
print(f"Selected UASG  : {len(sel)}")
print(f"Phenotyped IIDs: {iid_keep.size}")
print(f"Wrote PHENO PSAM: {PHENO_PSAM}")
print(f"Wrote IID keep  : {iid_keep_path}")


# ------------ Run (after QC) ------------
    # source scripts/00_config.sh

# Use clean set
    # python3 scripts/21_merge_pheno_to_psam.py

# Or target unrelateds
    # python3 scripts/21_merge_pheno_to_psam.py --psam "$PLINK_DIR/analysis.clean.unrel.psam"
