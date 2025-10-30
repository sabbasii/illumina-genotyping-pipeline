#!/usr/bin/env python3
"""
02_overlap_and_select_cohorts.py

- Intersect dfT2 (transpose_numbers + added samples) with df_snp by sample IDs
- Keep only samples whose Final Diagnosis is in {Control, Ischemic Stroke, TIA}
- Save:
  * final expression matrix (dfT2 subset with selected samples)
  * full metadata (df_csv columns 0..179) for selected samples
"""

import os
import pandas as pd

# ---- env (from 00_config.sh) ----
REPO_ROOT = os.environ.get("REPO_ROOT", os.getcwd())
RUN       = os.environ.get("RUN", "genotype_run1")
EXPR_OUT_DIR   = os.environ.get("EXPR_OUT_DIR", os.path.join(REPO_ROOT, "output", RUN, "expr", "explore"))
EXP_TGA_CSV    = os.environ.get("EXP_TGA_CSV")
EXP_TRANSPOSE  = os.environ.get("EXP_TRANSPOSE")  # original transpose (not needed directly here)
SAMPLE_SHEET   = os.environ.get("SAMPLE_SHEET")
DFT2_PATH      = os.path.join(EXPR_OUT_DIR, "dfT2_transpose_with_added_samples.csv")

# ---- outputs ----
SELECTED_EXPR_PATH   = os.path.join(EXPR_OUT_DIR, "expr_selected.csv") # final expression matrix
SELECTED_META_PATH   = os.path.join(EXPR_OUT_DIR, "meta_selected.csv") # full metadata (cols 0–179) for those samples
SELECTED_IDS_LIST    = os.path.join(EXPR_OUT_DIR, "ids_selected.txt") # sample IDs (for plink --keep)

# ---- constants ----
FIRST_SAMPLE_COL_IDX = 7
CSV_PROBE_START_NAME = "TC0100006437.hg.1"
KEEP_DIAG = {"control", "ischemic stroke", "tia"}  # case-insensitive compare

def read_auto(path, *, skiprows=None):
    return pd.read_csv(path, sep=None, engine="python", encoding="utf-8",
                       on_bad_lines="warn", skiprows=skiprows)

def norm_ids(x):
    return pd.Series(x, dtype="string").str.strip()

# ---- load inputs ----
if not os.path.isfile(DFT2_PATH):
    raise FileNotFoundError(f"dfT2 not found: {DFT2_PATH} (run 01_* prepare script first)")
dfT2 = read_auto(DFT2_PATH)

df_csv = read_auto(EXP_TGA_CSV)
if "UASG" not in df_csv.columns:
    raise KeyError("df_csv missing 'UASG'")
if CSV_PROBE_START_NAME not in df_csv.columns:
    raise KeyError(f"df_csv missing probe start '{CSV_PROBE_START_NAME}'")
csv_meta_end_idx = df_csv.columns.get_loc(CSV_PROBE_START_NAME)  # should be 180

df_snp = read_auto(SAMPLE_SHEET, skiprows=8)
df_snp.columns = df_snp.columns.str.strip()
if "Sample_Name" not in df_snp.columns:
    raise KeyError("df_snp missing 'Sample_Name'")

# ---- build ID sets ----
t2_ids = norm_ids(pd.Index(dfT2.columns[FIRST_SAMPLE_COL_IDX:]))   # samples in dfT2
snp_ids = norm_ids(df_snp["Sample_Name"])                           # samples in sample sheet
csv_ids = norm_ids(df_csv["UASG"])                                  # samples in df_csv

overlap_t2_snp = set(t2_ids) & set(snp_ids)

# ---- filter df_csv to allowed diagnoses ----
if "Final Diagnosis" not in df_csv.columns:
    raise KeyError("df_csv missing 'Final Diagnosis'")

df_csv["_UASG_"] = csv_ids
diag_norm = df_csv["Final Diagnosis"].astype("string").str.strip().str.lower()
is_keep_diag = diag_norm.isin(KEEP_DIAG)
df_csv_keep = df_csv.loc[is_keep_diag].copy()

# ---- final selected IDs = overlap ∩ keep_diag ----
selected_ids = sorted(overlap_t2_snp & set(df_csv_keep["_UASG_"]))
if not selected_ids:
    raise ValueError("No samples after overlap + diagnosis filter. Check labels in 'Final Diagnosis'.")

# ---- build final expression matrix (dfT2 subset) ----
keep_cols = list(dfT2.columns[:FIRST_SAMPLE_COL_IDX]) + selected_ids
expr_selected = dfT2.loc[:, keep_cols]
expr_selected.to_csv(SELECTED_EXPR_PATH, index=False)

# ---- build full metadata (cols 0..179) for selected samples ----
meta_cols = list(df_csv.columns[:csv_meta_end_idx])  # indices 0..179 inclusive of 0-based slice end-exclusive
meta_selected = (df_csv_keep[meta_cols + ["_UASG_"]]
                 .set_index("_UASG_")
                 .loc[selected_ids]
                 .reset_index()
                 .rename(columns={"_UASG_": "UASG"}))
meta_selected.to_csv(SELECTED_META_PATH, index=False)

# ---- save ID list (handy for plink --keep) ----
with open(SELECTED_IDS_LIST, "w", encoding="utf-8") as f:
    f.write("\n".join(selected_ids) + "\n")

print("=== Overlap + Selection Summary ===")
print(f"dfT2 samples         : {len(t2_ids)}")
print(f"SNP sheet samples    : {len(snp_ids)}")
print(f"Overlap (dfT2 ∩ SNP) : {len(overlap_t2_snp)}")
print(f"Selected diagnoses   : {sorted(KEEP_DIAG)}")
print(f"Selected samples     : {len(selected_ids)}")
print("")
print(f"Wrote expression  → {SELECTED_EXPR_PATH}")
print(f"Wrote metadata    → {SELECTED_META_PATH}")
print(f"Wrote ID list     → {SELECTED_IDS_LIST}")

# How to run (after sourcing your config):
# from repo root
    # source scripts/00_config.sh
    # python3 scripts/02_overlap_and_select_cohorts.py

# Three output files inside $EXPR_OUT_DIR:
    # expr_selected.csv — final expression matrix
    # meta_selected.csv — full metadata
    # ids_selected.txt — sample IDs (for plink --keep)