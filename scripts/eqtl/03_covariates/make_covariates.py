#!/usr/bin/env python3
"""
make_covariates.py

Build Matrix eQTL-style Covariates.txt from clinical metadata for the FINAL overlap samples.

Inputs (defaults)
-----------------
Overlap list:
  $REPO_ROOT/output/genotype_run1/eqtl/uasg_ge_snp_meta_overlap.txt
Metadata CSV:
  $REPO_ROOT/metadata/clinical_data.csv   (must contain 'sample_name')

Output
------
  $REPO_ROOT/output/genotype_run1/eqtl/Covariates.txt

Format (Matrix eQTL covariates)
-------------------------------
Rows = covariates, Cols = samples, first column = 'id'
"""

from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path
import pandas as pd


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--overlap", default=None, help="Path to overlap list (one sample per line)")
    p.add_argument("--meta", default=None, help="Path to clinical metadata CSV")
    p.add_argument("--sample-col", default="sample_name", help="Sample ID column in metadata")
    p.add_argument("--out", default=None, help="Output Covariates.txt path")
    return p.parse_args()


def read_overlap(path: Path) -> list[str]:
    keep: list[str] = []
    with path.open() as f:
        for line in f:
            s = line.strip().strip('"').rstrip("\r")
            if s:
                keep.append(s)
    if not keep:
        die(f"Overlap list is empty: {path}")
    return keep


def main() -> int:
    args = parse_args()

    repo_root = os.environ.get("REPO_ROOT", "")
    if not repo_root:
        die("REPO_ROOT not set. Run: source scripts/00_config.sh")

    base_eqtl = Path(repo_root) / "output/genotype_run1/eqtl"

    overlap_path = Path(args.overlap) if args.overlap else (base_eqtl / "uasg_ge_snp_meta_overlap.txt")
    meta_path = Path(args.meta) if args.meta else (Path(repo_root) / "metadata/clinical_data.csv")
    out_path = Path(args.out) if args.out else (base_eqtl / "Covariates.txt")

    if not overlap_path.exists():
        die(f"Missing overlap list: {overlap_path}")
    if not meta_path.exists():
        die(f"Missing metadata CSV: {meta_path}")

    samples = read_overlap(overlap_path)

    # Load metadata (keep types as strings; we’ll encode later)
    meta = pd.read_csv(meta_path, dtype=str)
    if args.sample_col not in meta.columns:
        die(f"Metadata must contain column '{args.sample_col}'. Found: {list(meta.columns)}")

    # Subset to overlap samples + preserve order
    meta[args.sample_col] = meta[args.sample_col].astype(str)
    meta_sub = meta[meta[args.sample_col].isin(samples)].copy()

    # Reindex to exact overlap order
    meta_sub = meta_sub.set_index(args.sample_col).reindex(samples)

    # Safety: if any overlap sample is missing from metadata, reindex creates NA rows
    missing_in_meta = meta_sub.index[meta_sub.isna().all(axis=1)].tolist()
    if missing_in_meta:
        die(
            "Some overlap samples are missing from metadata (would create empty covariate columns):\n"
            + "\n".join(missing_in_meta[:20])
            + ("\n..." if len(missing_in_meta) > 20 else "")
        )

    # Optional cleanup: drop columns that are entirely missing
    meta_sub = meta_sub.dropna(axis=1, how="all")

    # Transpose to Matrix eQTL covariate format
    cov = meta_sub.T
    cov.insert(0, "id", cov.index.astype(str))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    cov.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    print(f"Wrote: {out_path}")
    print(f"Samples (columns excluding id): {cov.shape[1] - 1}")
    print(f"Covariates (rows excluding header): {cov.shape[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Run:
source scripts/00_config.sh
python3 scripts/eqtl/03_covariates/make_covariates.py
"""