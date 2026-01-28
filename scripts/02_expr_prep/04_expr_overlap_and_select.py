#!/usr/bin/env python3
"""
04_expr_overlap_and_select.py

Select samples for downstream analysis.

What this step does:
- Intersect the prepared expression matrix (from step 03) with the SNP SampleSheet
- Keep only samples with allowed clinical diagnoses
- Write filtered expression matrix, metadata table, and ID list

Inputs are taken ONLY from environment variables (set by scripts/00_config.sh):
  EXPR_META_W   : path to wide expression CSV (metadata + probes; must include UASG + Final Diagnosis)
  SAMPLE_SHEET  : path to Illumina SampleSheet CSV

This step defines the final cohort used in later analyses.
"""

from __future__ import annotations

import os
import pandas as pd

from .utils.config import load_config
from .utils.ids import norm_ids
from .utils.io import ensure_file, read_auto, write_lines
from .utils.log import banner, kv
from .utils.paths import get_paths
from .utils.samplesheet import read_samplesheet


def _require_env_file(name: str) -> str:
    """Get env var value and ensure it points to an existing file."""
    val = os.environ.get(name, "").strip()
    if not val:
        raise EnvironmentError(
            f"Missing environment variable '{name}'. "
            f"Did you run: source scripts/00_config.sh ?"
        )
    return ensure_file(val)


def main() -> None:
    banner("04 — Overlap (expr ∩ SNP) and select by diagnosis")

    # ---- setup ----
    P = get_paths()
    cfg = load_config()

    FIRST_SAMPLE_COL_IDX = int(cfg["FIRST_SAMPLE_COL_IDX"])
    CSV_PROBE_START_NAME = str(cfg["CSV_PROBE_START_NAME"])
    KEEP_DIAG = {str(d).strip().lower() for d in cfg["KEEP_DIAG"]}
    SAMPLESHEET_SKIPROWS = int(cfg["SAMPLESHEET_SKIPROWS"])

    wide_path = _require_env_file("EXPR_META_W")
    samplesheet_path = _require_env_file("SAMPLE_SHEET")

    if not P.final_expr.is_file():
        raise FileNotFoundError(
            f"Missing final expression matrix: {P.final_expr}\n"
            f"Run step 03 first:\n"
            f"  python3 -m scripts.expr_prep.03_expr_metadata_prep"
        )

    kv("Wide expression", wide_path)
    kv("SampleSheet", samplesheet_path)
    kv("Final expr (03)", P.final_expr)

    # ---- load inputs ----
    # Step 03 output
    df_final = read_auto(P.final_expr, dtype=str)

    # Wide file (we use metadata columns and diagnosis)
    df_wide = read_auto(wide_path, dtype=str)

    # SampleSheet
    df_snp = read_samplesheet(samplesheet_path, skiprows=SAMPLESHEET_SKIPROWS)

    # ---- checks ----
    if "UASG" not in df_wide.columns:
        raise KeyError("Wide file missing required column: 'UASG'")
    if "Final Diagnosis" not in df_wide.columns:
        raise KeyError("Wide file missing required column: 'Final Diagnosis'")
    if CSV_PROBE_START_NAME not in df_wide.columns:
        raise KeyError(
            f"Probe start column '{CSV_PROBE_START_NAME}' not found in wide file"
        )

    if df_final.shape[1] <= FIRST_SAMPLE_COL_IDX:
        raise ValueError(
            "Final expression matrix does not contain sample columns at the expected index. "
            f"FIRST_SAMPLE_COL_IDX={FIRST_SAMPLE_COL_IDX}, n_cols={df_final.shape[1]}"
        )

    if "Sample_Name" not in df_snp.columns:
        raise KeyError("SampleSheet missing required column: 'Sample_Name'")

    # ---- identify metadata columns in wide file ----
    # Metadata columns are everything before the first probe column
    meta_end_idx = df_wide.columns.get_loc(CSV_PROBE_START_NAME)
    meta_cols = [str(c) for c in df_wide.columns[:meta_end_idx]]

    # IMPORTANT: avoid duplicate 'UASG' after rebuilding it from normalized IDs
    meta_cols = [c for c in meta_cols if c != "UASG"]

    # ---- normalize IDs across sources ----
    final_ids = set(norm_ids(df_final.columns[FIRST_SAMPLE_COL_IDX:]))
    snp_ids = set(norm_ids(df_snp["Sample_Name"]))
    wide_ids_norm = norm_ids(df_wide["UASG"])

    # overlap between expression (03) and SNP samplesheet
    overlap_final_snp = final_ids & snp_ids

    # ---- filter wide metadata by diagnosis ----
    df_wide["_UASG_NORM_"] = wide_ids_norm
    diag_norm = (
        df_wide["Final Diagnosis"]
        .astype("string")
        .str.strip()
        .str.lower()
    )
    df_wide_keep = df_wide.loc[diag_norm.isin(KEEP_DIAG)].copy()

    # final selected IDs must be in:
    # (expr ∩ SNP) AND diagnosis-kept metadata
    selected_ids = sorted(overlap_final_snp & set(df_wide_keep["_UASG_NORM_"]))

    if not selected_ids:
        raise ValueError(
            "No samples remain after overlap and diagnosis filtering.\n"
            "Things to check:\n"
            "  1) Do df_final sample column names match SampleSheet Sample_Name values?\n"
            "  2) Are 'Final Diagnosis' values spelled as expected in KEEP_DIAG?\n"
        )

    # ---- filter expression matrix ----
    keep_cols = list(df_final.columns[:FIRST_SAMPLE_COL_IDX]) + selected_ids
    expr_filtered = df_final.loc[:, keep_cols]
    expr_filtered.to_csv(P.filtered_expr, index=False)

    # ---- filter metadata (safe selection) ----
    meta_indexed = (
        df_wide_keep[meta_cols + ["_UASG_NORM_"]]
        .set_index("_UASG_NORM_")
    )

    # reindex is safer than .loc[] because it won't crash if something is missing
    meta_filtered = (
        meta_indexed.reindex(selected_ids)
        .reset_index()
        .rename(columns={"_UASG_NORM_": "UASG"})
    )

    # Count missing metadata rows robustly (works even if something weird happens)
    n_missing_meta = int(meta_filtered["UASG"].isna().to_numpy().sum()) if "UASG" in meta_filtered.columns else 0
    if n_missing_meta:
        banner("Warning")
        print(f"{n_missing_meta} selected IDs had no matching metadata rows (ID mismatch).")

    meta_filtered.to_csv(P.filtered_meta, index=False)

    # ---- write ID list ----
    write_lines(P.filtered_ids, selected_ids)

    # ---- summary ----
    banner("Overlap and selection summary")
    print(f"FINAL expr samples   : {len(final_ids)}")
    print(f"SNP sheet samples    : {len(snp_ids)}")
    print(f"Overlap (expr ∩ SNP) : {len(overlap_final_snp)}")
    print(f"Selected samples     : {len(selected_ids)}")
    print("")
    print(f"Wrote expression → {P.filtered_expr}")
    print(f"Wrote metadata   → {P.filtered_meta}")
    print(f"Wrote ID list    → {P.filtered_ids}")


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------
# How to run
#
# From the repository root:
#
#   source scripts/00_config.sh
#   python3 -m scripts.expr_prep.04_expr_overlap_and_select
#
# Outputs (under $EXPR_OUT_DIR):
#   matrices/expression_matrix_filtered.csv
#   metadata/metadata_filtered.csv
#   lists/ids_filtered.txt
# ---------------------------------------------------------------------