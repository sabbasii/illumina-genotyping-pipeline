#!/usr/bin/env python3
"""
05_make_cohort_id_lists.py
Map selected UASGs (ids_selected.txt) -> IIDs via SampleSheet, and write:
  - metadata/pheno/iid_selected.keep      (IIDs for bcftools/plink --keep/-S)
  - metadata/pheno/uasg_to_iid.tsv        (diagnostic mapping)
  - metadata/vcf.samples.cohort           (same IIDs list, useful for logs)

Env (set by scripts/00_config.sh):
  REPO_ROOT, RUN, SAMPLE_SHEET, EXPR_OUT_DIR, PHENO_DIR

Run:
  source scripts/00_config.sh
  python3 scripts/05_make_cohort_id_lists.py
"""

import os
import sys
import pandas as pd

def norm(s: pd.Series) -> pd.Series:
    return pd.Series(s, dtype="string").str.strip()

# --- Resolve paths from environment (with safe fallbacks) ---
REPO_ROOT   = os.environ.get("REPO_ROOT", os.getcwd())
RUN         = os.environ.get("RUN", "genotype_run1")
SAMPLE_SHEET= os.environ.get("SAMPLE_SHEET", "")
EXPR_OUT_DIR= os.environ.get("EXPR_OUT_DIR", os.path.join(REPO_ROOT, "output", RUN, "expr", "explore"))
PHENO_DIR   = os.environ.get("PHENO_DIR", os.path.join(REPO_ROOT, "metadata", "pheno"))

IDS_SELECTED = os.path.join(EXPR_OUT_DIR, "ids_selected.txt")  # UASG list from 04_expr_overlap_and_select.py
IID_KEEP_OUT = os.path.join(PHENO_DIR, "iid_selected.keep")
U2I_TSV_OUT  = os.path.join(PHENO_DIR, "uasg_to_iid.tsv")
COHORT_SAMPLES = os.path.join(REPO_ROOT, "metadata", "vcf.samples.cohort")

os.makedirs(PHENO_DIR, exist_ok=True)
os.makedirs(os.path.dirname(COHORT_SAMPLES), exist_ok=True)

# --- Guardrails ---
if not os.path.isfile(IDS_SELECTED):
    sys.exit(f"[ERR] ids_selected.txt not found: {IDS_SELECTED}\n       Run 04_expr_overlap_and_select.py first.")
if not os.path.isfile(SAMPLE_SHEET):
    sys.exit(f"[ERR] SAMPLE_SHEET not found: {SAMPLE_SHEET}\n       Check scripts/00_config.sh or env.")

# --- Load inputs ---
with open(IDS_SELECTED, "r", encoding="utf-8") as f:
    sel_uasg = [x.strip() for x in f if x.strip()]
sel_uasg_s = set(sel_uasg)

# SampleSheet header starts after 8 lines
ss = pd.read_csv(SAMPLE_SHEET, sep=None, engine="python", skiprows=8, dtype=str)
ss.columns = ss.columns.str.strip()

needed = ["Sample_Name", "SentrixBarcode_A", "SentrixPosition_A"]
missing = [c for c in needed if c not in ss.columns]
if missing:
    sys.exit(f"[ERR] SampleSheet missing columns: {missing}")

ss["Sample_Name"]      = norm(ss["Sample_Name"])
ss["SentrixBarcode_A"] = norm(ss["SentrixBarcode_A"])
ss["SentrixPosition_A"]= norm(ss["SentrixPosition_A"])
ss["IID"]              = ss["SentrixBarcode_A"] + "_" + ss["SentrixPosition_A"]
ss["UASG"]             = ss["Sample_Name"]

# Deduplicate any repeated UASG rows favoring first occurrence, then unique (UASG -> IID) map
ss_u = ss.dropna(subset=["UASG", "IID"]).drop_duplicates(subset=["UASG"], keep="first")
u2i = ss_u.set_index("UASG")["IID"].to_dict()

# --- Map and validate ---
mapped = []
missing_uasg = []
dup_iid_counts = {}

for u in sel_uasg:
    iid = u2i.get(u)
    if not iid or pd.isna(iid) or iid == "":
        missing_uasg.append(u)
        continue
    mapped.append((u, iid))
    dup_iid_counts[iid] = dup_iid_counts.get(iid, 0) + 1

mapped_uasg = [u for u, _ in mapped]
mapped_iids = [iid for _, iid in mapped]

# Detect duplicate IIDs (should not happen; warn if it does)
dup_iids = sorted([iid for iid, c in dup_iid_counts.items() if c > 1])

# --- Write outputs ---
# 1) keep list (IIDs only, one per line)
with open(IID_KEEP_OUT, "w", encoding="utf-8") as f:
    for iid in mapped_iids:
        f.write(f"{iid}\n")

# 2) mapping table (UASG\tIID)
pd.DataFrame(mapped, columns=["UASG", "IID"]).to_csv(U2I_TSV_OUT, sep="\t", index=False)

# 3) cohort samples convenience copy (same content as keep list)
with open(COHORT_SAMPLES, "w", encoding="utf-8") as f:
    for iid in mapped_iids:
        f.write(f"{iid}\n")

# --- Summary ---
print("=== 05_make_cohort_id_lists.py ===")
print(f"RUN                : {RUN}")
print(f"Selected UASGs in  : {IDS_SELECTED}")
print(f"SampleSheet        : {SAMPLE_SHEET}")
print(f"Mapped IIDs (keep) : {len(mapped_iids)}")
print(f"Wrote keep list    : {IID_KEEP_OUT}")
print(f"Wrote map (UASG→IID): {U2I_TSV_OUT}")
print(f"Wrote cohort list  : {COHORT_SAMPLES}")

if missing_uasg:
    print(f"[WARN] UASGs with no IID in SampleSheet ({len(missing_uasg)}):")
    for u in missing_uasg[:10]:
        print("  -", u)
    if len(missing_uasg) > 10:
        print("  ...", len(missing_uasg) - 10, "more")

if dup_iids:
    print(f"[WARN] Duplicate IIDs mapped from multiple UASGs ({len(dup_iids)}). Check {U2I_TSV_OUT}")
    for iid in dup_iids[:10]:
        print("  -", iid, f"(count={dup_iid_counts[iid]})")
    if len(dup_iids) > 10:
        print("  ...", len(dup_iids) - 10, "more")

# Helpful sanity assertions (soft)
if len(mapped_iids) != len(set(mapped_iids)):
    print("[WARN] keep list contains duplicate IIDs (non-fatal but unexpected).")
if len(mapped_iids) == 0:
    sys.exit("[ERR] No IIDs mapped. Aborting.")

# --------------- RUN ---------------
# python3 scripts/05_make_cohort_id_lists.py