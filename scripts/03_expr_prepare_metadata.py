#!/usr/bin/env python3
"""
01_load_and_prepare_microarray_expression_metadata.py

Step 1+2: Load inputs (df_csv, df_transpose, df_snp), validate,
extract mRS (row0) from df_transpose, and build dfT2 by adding samples
present in df_csv but missing from df_transpose.
Outputs go to $EXPR_OUT_DIR.
"""

import os
import pandas as pd
from datetime import datetime

# ---- env (set by 00_config.sh) ----
REPO_ROOT = os.environ.get("REPO_ROOT", os.getcwd())
RUN       = os.environ.get("RUN", "genotype_run1")
EXP_TGA_CSV   = os.environ.get("EXP_TGA_CSV")
EXP_TRANSPOSE = os.environ.get("EXP_TRANSPOSE")
SAMPLE_SHEET  = os.environ.get("SAMPLE_SHEET")
EXPR_OUT_DIR  = os.environ.get("EXPR_OUT_DIR", os.path.join(REPO_ROOT, "output", RUN, "expr", "explore"))
os.makedirs(EXPR_OUT_DIR, exist_ok=True)

REPORT_PATH   = os.path.join(EXPR_OUT_DIR, "01_prepare_report.txt")
MRS_PATH      = os.path.join(EXPR_OUT_DIR, "mrs_outcome_by_sample.csv")
DFT2_PATH     = os.path.join(EXPR_OUT_DIR, "dfT2_transpose_with_added_samples.csv")

CSV_PROBE_START_NAME = "TC0100006437.hg.1"
FIRST_SAMPLE_COL_IDX = 7

def read_auto(path, *, skiprows=None):
    if not path or not os.path.isfile(path):
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep=None, engine="python", encoding="utf-8", on_bad_lines="warn", skiprows=skiprows)

def norm_ids(x):
    return pd.Series(x, dtype="string").str.strip()

def write_report(lines):
    with open(REPORT_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

# ---- load ----
df_csv       = read_auto(EXP_TGA_CSV)
df_transpose = read_auto(EXP_TRANSPOSE)
df_snp       = read_auto(SAMPLE_SHEET, skiprows=8)
df_snp.columns = df_snp.columns.str.strip()

# ---- checks ----
if "UASG" not in df_csv.columns:
    raise KeyError("df_csv missing 'UASG'")
if CSV_PROBE_START_NAME not in df_csv.columns:
    raise KeyError(f"df_csv missing first probe column '{CSV_PROBE_START_NAME}'")
csv_expr_start_idx = df_csv.columns.get_loc(CSV_PROBE_START_NAME)

if "ID" not in df_transpose.columns:
    raise KeyError("df_transpose missing 'ID'")
if "Sample_Name" not in df_snp.columns:
    raise KeyError("df_snp missing 'Sample_Name' (ensure skiprows=8)")

# ---- extract row0 clinical (e.g., mRS good/bad) ----
if df_transpose.shape[0] < 1 or df_transpose.shape[1] <= FIRST_SAMPLE_COL_IDX:
    raise ValueError("df_transpose shape unexpected; cannot extract row0")
mrs = df_transpose.iloc[0, FIRST_SAMPLE_COL_IDX:]
mrs.index = norm_ids(mrs.index)
mrs.name = "mRS_outcome"
mrs.to_csv(MRS_PATH, header=True)

# ---- build dfT2: add samples present in df_csv but missing in df_transpose ----
csv_ids = set(norm_ids(df_csv["UASG"]))
t_ids   = set(norm_ids(pd.Index(df_transpose.columns[FIRST_SAMPLE_COL_IDX:])))

only_in_csv = sorted(csv_ids - t_ids)  # typically ['UASG-1436','UASG-0256'] in your notes
dfT2 = df_transpose.copy(deep=True)

# probe order from df_transpose (rows 1+: expression)
probe_order = dfT2.loc[1:, "ID"].astype(str).str.strip().tolist()

# expression columns in df_csv (drop helpers like UASG_clean if present)
csv_expr_cols = df_csv.columns[csv_expr_start_idx:]
csv_expr_cols = [c for c in csv_expr_cols if c in set(probe_order)]

df_csv["_UASG_idx_"] = df_csv["UASG"].astype(str).str.strip()
csv_by_uasg = df_csv.set_index("_UASG_idx_", drop=False)

added, skipped = [], []
for sid in only_in_csv:
    if sid not in csv_by_uasg.index:
        skipped.append((sid, "not found in df_csv"))
        continue
    # vector of expression values aligned to probe_order
    row = csv_by_uasg.loc[sid, csv_expr_cols]
    expr_vec = row.reindex(probe_order)
    # create the new column in dfT2
    dfT2[sid] = pd.NA
    dfT2.loc[1:, sid] = expr_vec.values
    added.append(sid)

# ---- overlaps (for info only at this step) ----
snp_ids = set(norm_ids(df_snp["Sample_Name"]))
overlap_csv_t   = csv_ids & t_ids
overlap_all     = overlap_csv_t & snp_ids
t2_ids_final    = set(norm_ids(pd.Index(dfT2.columns[FIRST_SAMPLE_COL_IDX:])))

# ---- save outputs ----
dfT2.to_csv(DFT2_PATH, index=False)

# ---- report ----
lines = []
lines.append(f"Prepare Report — {datetime.now().isoformat(timespec='seconds')}")
lines.append(f"RUN={RUN}")
lines.append("")
lines.append("Files:")
lines.append(f"  df_csv        : {EXP_TGA_CSV}")
lines.append(f"  df_transpose  : {EXP_TRANSPOSE}")
lines.append(f"  df_snp        : {SAMPLE_SHEET}")
lines.append("")
lines.append("Shapes:")
lines.append(f"  df_csv        : {df_csv.shape}")
lines.append(f"  df_transpose  : {df_transpose.shape}")
lines.append(f"  df_snp        : {df_snp.shape}")
lines.append("")
lines.append("Key checks:")
lines.append(f"  df_csv probe start name/index : '{CSV_PROBE_START_NAME}' / {csv_expr_start_idx}")
lines.append(f"  first sample col in transpose : {FIRST_SAMPLE_COL_IDX}")
lines.append("")

lines.append("ID counts (before add):")
lines.append(f"  csv_ids={len(csv_ids)}  t_ids={len(t_ids)}  snp_ids={len(snp_ids)}")
lines.append(f"  overlap csv∩t={len(overlap_csv_t)}  overlap all={len(overlap_all)}")
lines.append("")
lines.append(f"Added samples from df_csv not in df_transpose: {len(added)} -> {added}")
if skipped:
    lines.append(f"Skipped (not found): {skipped}")
lines.append(f"T2 final sample count (cols ≥ {FIRST_SAMPLE_COL_IDX}): {len(t2_ids_final)}")
lines.append("")
lines.append(f"Wrote mRS outcomes: {MRS_PATH}")
lines.append(f"Wrote dfT2:         {DFT2_PATH}")

write_report(lines)
print("\n".join(lines))


# How to run (after sourcing your config):
# from repo root
    # source scripts/00_config.sh
    # python3 scripts/expr_prep/01_load_and_check_microarray_expression_metadata.py

# Two output files:
    # 01_load_and_check_report.txt (in $EXPR_OUT_DIR) — a quick report: file paths, shapes, detected probe-start index, key checks, and overlap counts (df_csv, df_transpose, df_snp).
    # mrs_outcome_by_sample.csv

