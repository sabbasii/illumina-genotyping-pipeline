#!/usr/bin/env python3
"""
build_GE_from_expr.py

Build Matrix-eQTL expression matrix (GE.txt) from expr_table.tsv.

Input (tab-delimited):
  input_data/expr_array/expr_table.tsv
    - Must include column: SYMBOL
    - Must include sample columns like: UASG_####

Output (tab-delimited):
  output/genotype_run1/eqtl/GE.txt
    - First column: geneid (SYMBOL)
    - Remaining columns: UASG_#### samples
    - Duplicate SYMBOL rows are collapsed by MAX per sample (not mean).
"""

from __future__ import annotations

import os
import sys
import pandas as pd


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def main() -> int:
    repo_root = os.environ.get("REPO_ROOT", "")
    if not repo_root:
        die("REPO_ROOT is not set. Run: source scripts/00_config.sh")

    in_tsv = os.path.join(repo_root, "input_data", "expr_array", "expr_table.tsv")
    out_txt = os.path.join(repo_root, "output", "genotype_run1", "eqtl", "GE.txt")
    os.makedirs(os.path.dirname(out_txt), exist_ok=True)

    if not os.path.isfile(in_tsv):
        die(f"Missing input: {in_tsv}")

    df = pd.read_csv(in_tsv, sep="\t", dtype=str)  # read as str first; convert samples after
    if "SYMBOL" not in df.columns:
        die("expr_table.tsv must have a 'SYMBOL' column.")

    uasg_cols = [c for c in df.columns if c.startswith("UASG_")]
    if not uasg_cols:
        die("No sample columns found starting with 'UASG_'.")

    # Keep only needed columns
    ge = df[["SYMBOL"] + uasg_cols].copy()
    ge = ge.rename(columns={"SYMBOL": "geneid"})

    # Drop empty geneid
    ge["geneid"] = ge["geneid"].astype(str).str.strip()
    ge = ge.dropna(subset=["geneid"])
    ge = ge[ge["geneid"] != ""]

    # Convert sample columns to numeric (non-numeric -> NA)
    for c in uasg_cols:
        ge[c] = pd.to_numeric(ge[c], errors="coerce")

    # Collapse duplicate SYMBOL rows by MAX per sample
    # (If all values are NA for a gene/sample, result stays NA)
    ge = ge.groupby("geneid", as_index=False)[uasg_cols].max()

    # Write
    ge.to_csv(out_txt, sep="\t", index=False, na_rep="NA")
    print(f"Wrote {out_txt} with shape {ge.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Run
# source scripts/00_config.sh
# python scripts/eqtl/build_GE_from_expr.py