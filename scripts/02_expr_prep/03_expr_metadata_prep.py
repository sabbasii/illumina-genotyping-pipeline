#!/usr/bin/env python3
"""
03_expr_metadata_prep.py

Prepare expression metadata and matrices (alignment + missing-sample fill).

What this step does:
- Loads:
    1) wide expression CSV (must include a 'UASG' column)
    2) transposed expression matrix CSV (must include an 'ID' column)
    3) Illumina SNP SampleSheet (for reporting/consistency checks)
- Extracts row-0 clinical values (e.g., mRS) from the transposed matrix
- Adds samples that exist in the wide file but are missing from the transposed matrix
- Writes:
    - clinical/mrs_outcome_by_sample.csv
    - matrices/expression_matrix_transposed_FINAL.csv
    - reports/expr_metadata_check_report.txt

Inputs are taken ONLY from environment variables (set by scripts/00_config.sh):
  EXPR_META_W     : path to wide expression CSV
  EXPR_MATRIX_T   : path to transposed expression CSV
  SAMPLE_SHEET    : path to Illumina SampleSheet CSV
"""

from __future__ import annotations

import os
from datetime import datetime

import pandas as pd

from .utils.config import load_config
from .utils.ids import norm_ids
from .utils.io import ensure_file, read_auto
from .utils.log import banner, kv, write_report
from .utils.paths import get_paths
from .utils.samplesheet import read_samplesheet


def _require_env(name: str) -> str:
    """Return env var value or raise a clear error."""
    val = os.environ.get(name, "").strip()
    if not val:
        raise EnvironmentError(
            f"Missing environment variable '{name}'. "
            f"Did you run: source scripts/00_config.sh ?"
        )
    return ensure_file(val)


def main() -> None:
    banner("03 — Expression metadata prep (align + fill missing samples)")

    # ---- setup ----
    P = get_paths()
    cfg = load_config()

    FIRST_SAMPLE_COL_IDX = int(cfg["FIRST_SAMPLE_COL_IDX"])
    CSV_PROBE_START_NAME = str(cfg["CSV_PROBE_START_NAME"])
    SAMPLESHEET_SKIPROWS = int(cfg["SAMPLESHEET_SKIPROWS"])

    wide_path = _require_env("EXPR_META_W")
    transpose_path = _require_env("EXPR_MATRIX_T")
    samplesheet_path = _require_env("SAMPLE_SHEET")

    kv("Wide expression", wide_path)
    kv("Transposed matrix", transpose_path)
    kv("SampleSheet", samplesheet_path)
    kv("Output dir", P.expr_out_dir)

    # ---- load inputs ----
    df_wide = read_auto(wide_path, dtype=str)
    df_transpose = read_auto(transpose_path, dtype=str)
    df_snp = read_samplesheet(samplesheet_path, skiprows=SAMPLESHEET_SKIPROWS)

    # ---- basic checks ----
    if "UASG" not in df_wide.columns:
        raise KeyError("Wide expression file is missing required column: 'UASG'")

    if CSV_PROBE_START_NAME not in df_wide.columns:
        raise KeyError(
            f"Probe start column '{CSV_PROBE_START_NAME}' not found in wide file."
        )

    if "ID" not in df_transpose.columns:
        raise KeyError("Transposed matrix is missing required column: 'ID'")

    if df_transpose.shape[1] <= FIRST_SAMPLE_COL_IDX:
        raise ValueError(
            "Transposed matrix does not contain sample columns at the expected index. "
            f"FIRST_SAMPLE_COL_IDX={FIRST_SAMPLE_COL_IDX}, n_cols={df_transpose.shape[1]}"
        )

    # ---- extract row-0 clinical values (e.g. mRS) ----
    # Row 0 holds clinical values for sample columns (as per your dataset convention)
    mrs = df_transpose.iloc[0, FIRST_SAMPLE_COL_IDX:]
    mrs.index = norm_ids(mrs.index)
    mrs.name = "mRS_outcome"
    mrs.to_csv(P.mrs_csv, header=True)

    # ---- determine probe order ----
    csv_expr_start_idx = df_wide.columns.get_loc(CSV_PROBE_START_NAME)
    probe_order = (
        df_transpose.loc[1:, "ID"]
        .astype(str)
        .str.strip()
        .tolist()
    )

    # ---- determine sample IDs present in each file ----
    wide_ids = set(norm_ids(df_wide["UASG"]))
    transpose_ids = set(norm_ids(df_transpose.columns[FIRST_SAMPLE_COL_IDX:]))
    only_in_wide = sorted(wide_ids - transpose_ids)

    # ---- prepare wide expression subset ----
    df_wide["_UASG_NORM_"] = norm_ids(df_wide["UASG"])
    df_wide = df_wide.set_index("_UASG_NORM_", drop=False)

    # Expression columns in the wide file start at CSV_PROBE_START_NAME
    wide_expr_cols_all = [str(c).strip() for c in df_wide.columns[csv_expr_start_idx:]]
    wide_expr_cols = [c for c in wide_expr_cols_all if c in set(probe_order)]

    if not wide_expr_cols:
        raise ValueError(
            "No overlapping probe columns found between wide file and transposed matrix. "
            "Check CSV_PROBE_START_NAME and the probe IDs in both files."
        )

    # ---- add missing samples (only those present in wide but missing in transpose) ----
    df_final = df_transpose.copy(deep=True)

    added: list[str] = []
    skipped: list[tuple[str, str]] = []

    for sid in only_in_wide:
        if sid not in df_wide.index:
            skipped.append((sid, "not found in wide file index"))
            continue

        # Reindex to probe_order so probe alignment matches transposed matrix
        expr_vec = (
            df_wide.loc[sid, wide_expr_cols]
            .reindex(probe_order)
            .to_numpy()
        )

        # Add new sample column and fill rows 1: with expression values
        df_final[sid] = pd.NA
        df_final.loc[1:, sid] = expr_vec
        added.append(sid)

    # ---- save final matrix ----
    df_final.to_csv(P.final_expr, index=False)

    # ---- build report ----
    report_lines = [
        f"Expression prep report — {datetime.now().isoformat(timespec='seconds')}",
        "",
        "Inputs:",
        f"  wide expression : {wide_path}",
        f"  transposed      : {transpose_path}",
        f"  SNP samplesheet : {samplesheet_path}",
        "",
        "Input shapes:",
        f"  wide file       : {df_wide.shape}",
        f"  transposed      : {df_transpose.shape}",
        f"  SNP sheet       : {df_snp.shape}",
        "",
        "Sample counts:",
        f"  wide samples       : {len(wide_ids)}",
        f"  transposed samples : {len(transpose_ids)}",
        f"  only in wide       : {len(only_in_wide)}",
        f"  added samples      : {len(added)}",
        "",
        f"Added sample IDs ({len(added)}): {added}",
    ]

    if skipped:
        report_lines += ["", f"Skipped ({len(skipped)}): {skipped}"]

    report_lines += [
        "",
        "Outputs:",
        f"  mRS file        : {P.mrs_csv}",
        f"  final matrix    : {P.final_expr}",
        f"  report          : {P.report_txt}",
    ]

    write_report(P.report_txt, report_lines)

    # Also print the report to console
    print("\n" + "\n".join(report_lines))


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------
# How to run
#
# From the repository root:
#
#   source scripts/00_config.sh
#   python3 -m scripts.expr_prep.03_expr_metadata_prep
#
# Outputs (under $EXPR_OUT_DIR):
#   clinical/mrs_outcome_by_sample.csv
#   matrices/expression_matrix_transposed_FINAL.csv
#   reports/expr_metadata_check_report.txt
# ---------------------------------------------------------------------