#!/usr/bin/env python3
# 05_build_cohort_pheno_psam.py
# Build a strict PSAM with 7 tab-delimited columns:
# #FID  IID  SEX  UASG  StrokeStatus  Final_Diagnosis  PHENO1
#
# Inputs (from env/repo layout):
#   - $REPO_ROOT/scripts/00_config.sh provides paths
#   - $PHENO_DIR/iid_selected.keep (IIDs for cohort; one per line)
#   - SampleSheet.csv (skiprows=8; needs SentrixBarcode_A, SentrixPosition_A, Sample_Name, Gender/Sex)
#   - $EXPR_OUT_DIR/meta_selected.csv (must include UASG, Final Diagnosis, StrokeStatus [preferred])
# Optional:
#   - $PLINK_DIR/analysis.clean*.psam (used as a base if present; else we build base from keep list)
#
# Output:
#   - $PHENO_DIR/cohort.pheno.psam (7 columns, tabs, LF line endings)

import os, glob, sys
import pandas as pd

def valid_repo(root: str) -> bool:
    return bool(
        root
        and os.path.isdir(root)
        and os.path.isfile(os.path.join(root, "scripts", "00_config.sh"))
        and os.path.isdir(os.path.join(root, "output"))
    )

# --- Resolve repo root robustly (ENV if valid; else fallback to script dir/..)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
env_repo = os.environ.get("REPO_ROOT", "").strip()
fallback_repo = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
if valid_repo(env_repo):
    REPO = env_repo
    repo_note = f"Using REPO_ROOT from environment: {REPO}"
else:
    REPO = fallback_repo
    repo_note = (
        f"Ignoring invalid REPO_ROOT='{env_repo}' and using script-based repo: {REPO}"
        if env_repo else f"Using script-based repo: {REPO}"
    )

RUN   = os.environ.get("RUN", "genotype_run1").strip() or "genotype_run1"
OUT   = os.path.join(REPO, "output", RUN)
PLINK = os.path.join(OUT, "qc", "plink")
EXPR  = os.path.join(OUT, "expr", "explore")
META  = os.path.join(REPO, "metadata")
PHENO = os.path.join(META, "pheno")
os.makedirs(PHENO, exist_ok=True)

# SampleSheet: trust env if file exists; else fallback path
SSHEET_ENV = os.environ.get("SAMPLE_SHEET", "").strip()
SSHEET_FALLBACK = os.path.join(
    REPO, "input_data", "sample_sheet",
    "All_samples_Examine_SNPs_GWAS studies_GJ-P01_infiniumSampleSheet.csv"
)
SSHEET = SSHEET_ENV if (SSHEET_ENV and os.path.isfile(SSHEET_ENV)) else SSHEET_FALLBACK
if not os.path.isfile(SSHEET):
    raise FileNotFoundError(f"SampleSheet not found: {SSHEET}")

# Cohort keep list (IIDs)
iid_keep_path = os.path.join(PHENO, "iid_selected.keep")
if not os.path.isfile(iid_keep_path):
    raise FileNotFoundError(f"Missing keep list: {iid_keep_path} (run 02_overlap_and_select_cohorts.py first)")
keep = pd.read_csv(iid_keep_path, header=None, names=["IID"], dtype=str)
keep["IID"] = keep["IID"].str.strip()
if keep["IID"].duplicated().any():
    raise ValueError("Duplicate IIDs in iid_selected.keep")

# Selected metadata (phenotype)
meta_path = os.path.join(EXPR, "meta_selected.csv")
if not os.path.isfile(meta_path):
    raise FileNotFoundError(f"Missing metadata: {meta_path} (run 02_overlap_and_select_cohorts.py first)")
meta = pd.read_csv(meta_path, dtype=str).fillna("")
meta.columns = meta.columns.str.strip()

# Normalize/derive StrokeStatus if needed
if "UASG" not in meta.columns:
    raise KeyError("meta_selected.csv missing 'UASG'")
if "Final Diagnosis" in meta.columns and "Final_Diagnosis" not in meta.columns:
    meta = meta.rename(columns={"Final Diagnosis": "Final_Diagnosis"})
if "StrokeStatus" not in meta.columns:
    diag_norm = meta["Final_Diagnosis"].astype(str).str.strip().str.lower()
    stroke = pd.Series(pd.NA, index=meta.index, dtype="string")
    stroke = stroke.mask(diag_norm.eq("control"), "Control")
    stroke = stroke.mask(diag_norm.isin(["ischemic stroke", "tia"]), "Case")
    meta["StrokeStatus"] = stroke.fillna("")
# Keep only needed phenotype cols
meta = meta[["UASG", "StrokeStatus", "Final_Diagnosis"]].copy()

# SampleSheet → map IID↔UASG and SEX (1=male, 2=female, 0=unknown)
ss = pd.read_csv(SSHEET, sep=None, engine="python", skiprows=8, dtype=str)
ss.columns = ss.columns.str.strip()
for col in ("SentrixBarcode_A", "SentrixPosition_A", "Sample_Name"):
    if col not in ss.columns:
        raise KeyError(f"SampleSheet missing required column: {col}")

# Determine the gender/sex column name
sex_col = None
for candidate in ("Gender", "gender", "Sex", "sex"):
    if candidate in ss.columns:
        sex_col = candidate
        break
if sex_col is None:
    # If truly absent, treat all as unknown
    ss["__SEX_RAW__"] = ""
    sex_col = "__SEX_RAW__"

ss["IID"]  = ss["SentrixBarcode_A"].astype(str).str.strip() + "_" + ss["SentrixPosition_A"].astype(str).str.strip()
ss["UASG"] = ss["Sample_Name"].astype(str).str.strip()
sex_raw = ss[sex_col].astype(str).str.strip().str.lower()
sex_num = (sex_raw.map({"m":1, "male":1, "f":2, "female":2}).fillna(0)).astype(int)
map_iid = ss[["IID", "UASG"]].drop_duplicates()
map_sex = pd.DataFrame({"IID": ss["IID"], "SEX": sex_num}).drop_duplicates()

# Base PSAM: prefer newest analysis.clean*.psam; else build minimal from keep list
psam_candidates = []
base_psam = os.path.join(PLINK, "analysis.clean.psam")
if os.path.isfile(base_psam):
    psam_candidates.append(base_psam)
psam_candidates.extend(sorted(glob.glob(os.path.join(PLINK, "analysis.clean_*.psam"))))

if psam_candidates:
    psam_path = max(psam_candidates, key=os.path.getmtime)
    base = pd.read_csv(psam_path, sep="\t", dtype=str)
    base.columns = base.columns.str.strip()
    if "#FID" not in base.columns and "FID" in base.columns:
        base = base.rename(columns={"FID": "#FID"})
    if "#FID" not in base.columns:
        base.insert(0, "#FID", base["IID"].astype(str))
    if "IID" not in base.columns:
        raise KeyError("Base PSAM missing IID")
    base = base[["#FID", "IID"]].copy()
else:
    psam_path = "<built-from-keep>"
    base = keep.copy()
    base.insert(0, "#FID", base["IID"])

# Restrict to cohort IIDs, preserve keep order
base = base.merge(keep, on="IID", how="inner").drop_duplicates(subset=["IID"])
# Merge sex from SampleSheet (override any prior sex), default to 0
base = base.merge(map_sex, on="IID", how="left")
base["SEX"] = base["SEX"].fillna(0).astype(int)

# Merge IID↔UASG, then phenotype by UASG
ps = base.merge(map_iid, on="IID", how="left")
ps = ps.merge(meta, on="UASG", how="left")

# PHENO1: Control→1, Case→2, else 0
ps["PHENO1"] = ps["StrokeStatus"].map({"Control": 1, "Case": 2}).fillna(0).astype(int)

# Final 7 columns
final = ps[["#FID", "IID", "SEX", "UASG", "StrokeStatus", "Final_Diagnosis", "PHENO1"]].copy()

# Fill missing text with empty strings (SEX/PHENO1 are already numeric)
for col in ["UASG", "StrokeStatus", "Final_Diagnosis"]:
    final[col] = final[col].fillna("")

# Assertions / sanity checks
assert final.shape[1] == 7, "Output PSAM must have exactly 7 columns"
assert final["IID"].nunique() == len(final), "Duplicate IIDs in output PSAM"
assert len(final) == len(keep), f"Row count mismatch: PSAM={len(final)} vs keep={len(keep)}"

# Write with strict tabs + LF
out_path = os.path.join(PHENO, "cohort.pheno.psam")
final.to_csv(out_path, sep="\t", index=False, lineterminator="\n")

# Quick summary
n_missing_uasg = (final["UASG"] == "").sum()
n_pheno = final["PHENO1"].isin([1, 2]).sum()
print(repo_note)
print("RUN label     :", RUN)
print("Base PSAM     :", psam_path)
print("Wrote         :", out_path)
print("Rows          :", len(final))
print("Phenotyped    :", int(n_pheno))
print("Missing UASG  :", int(n_missing_uasg))
