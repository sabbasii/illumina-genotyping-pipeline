#!/usr/bin/env python3
"""
preprocess_covariates.py

Convert Covariates.txt (mixed types) into a strictly numeric covariate matrix
that MatrixEQTL can load.

Input (default)
---------------
  $REPO_ROOT/output/genotype_run1/eqtl/Covariates.txt

Outputs
-------
  $REPO_ROOT/output/genotype_run1/eqtl/Covariates_numeric.txt
  $REPO_ROOT/output/genotype_run1/eqtl/covariates_preprocess_summary.txt

Rules (matches your prior R logic)
----------------------------------
- Encode common values:
    sex: female/f -> 0, male/m -> 1
    yes/no: no/n/false/f -> 0, yes/y/true/t -> 1
- One-hot encode:
    diagnosis -> dx_*  (reference dropped: 'control' if present else most frequent)
    ancestry  -> anc_* (reference dropped: most frequent)
- Drop covariate rows that become all-NA after conversion
- Preserve sample column order exactly
"""

from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path
import pandas as pd

# ---- make scripts/ importable so "import eqtl..." works ----
HERE = Path(__file__).resolve()
SCRIPTS_DIR = HERE.parents[2]  # .../scripts
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from eqtl.utils.covariate_helpers import sanitize_covariates_numeric


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="in_path", default=None, help="Input Covariates.txt path")
    p.add_argument("--out", dest="out_path", default=None, help="Output Covariates_numeric.txt path")
    p.add_argument("--summary", dest="sum_path", default=None, help="Output summary txt path")
    p.add_argument("--no-ancestry", action="store_true", help="Disable ancestry one-hot encoding")
    p.add_argument("--no-diagnosis", action="store_true", help="Disable diagnosis one-hot encoding")
    p.add_argument("--dx-ref", default="control", help="Reference diagnosis level to drop (default: control)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    repo_root = os.environ.get("REPO_ROOT", "")
    if not repo_root:
        die("REPO_ROOT not set. Run: source scripts/00_config.sh")

    base_eqtl = Path(repo_root) / "output/genotype_run1/eqtl"

    in_path = Path(args.in_path) if args.in_path else (base_eqtl / "Covariates.txt")
    out_path = Path(args.out_path) if args.out_path else (base_eqtl / "Covariates_numeric.txt")
    sum_path = Path(args.sum_path) if args.sum_path else (base_eqtl / "covariates_preprocess_summary.txt")

    if not in_path.exists():
        die(f"Missing input: {in_path}")

    # Covariates.txt: col1=id; remaining columns are samples; rows are covariates
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    if "id" not in df.columns:
        die("Covariates.txt must have a first column named 'id'.")

    df = df.set_index("id")
    if df.shape[1] == 0:
        die("Covariates.txt has no sample columns.")

    # Core conversion
    cov_num, summary_lines = sanitize_covariates_numeric(
        cov=df,
        include_diagnosis=(not args.no_diagnosis),
        include_ancestry=(not args.no_ancestry),
        diagnosis_ref=args.dx_ref,
    )

    # Write Covariates_numeric.txt in Matrix eQTL covariate format (id + samples)
    cov_out = cov_num.copy()
    cov_out.insert(0, "id", cov_out.index.astype(str))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    cov_out.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    # Summary
    lines = []
    lines.append(f"Input:  {in_path}")
    lines.append(f"Output: {out_path}")
    lines.append(f"Samples: {cov_num.shape[1]}")
    lines.append(f"Covariate rows (numeric): {cov_num.shape[0]}")
    lines.append(f"Includes diagnosis one-hot: {not args.no_diagnosis}")
    lines.append(f"Includes ancestry one-hot:  {not args.no_ancestry}")
    lines.append("")
    lines.extend(summary_lines)

    sum_path.write_text("\n".join(lines).rstrip() + "\n")

    print(f"Wrote:\n  {out_path}\n  {sum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Run:
source scripts/00_config.sh
python3 scripts/eqtl/03_covariates/preprocess_covariates.py
"""