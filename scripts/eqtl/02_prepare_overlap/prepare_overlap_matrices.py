#!/usr/bin/env python3
"""
prepare_overlap_matrices.py

Purpose
-------
Prepare the final, aligned Matrix-eQTL matrices by:
1) Finding the sample overlap across:
     - GE.txt header samples
     - SNP.txt header samples
     - clinical metadata sample_name
2) Writing the final overlap list (preserving sample order)
3) Writing GE_overlap.txt and SNP_overlap.txt (subset to overlap samples)
4) Writing overlap_summary.txt with counts

Default inputs
--------------
GE   : $REPO_ROOT/output/genotype_run1/eqtl/GE.txt
SNP  : $REPO_ROOT/output/genotype_run1/eqtl/SNP.txt
META : $REPO_ROOT/metadata/clinical_data.csv  (must contain 'sample_name')

Default outputs (in eqtl dir)
-----------------------------
uasg_ge_snp_meta_overlap.txt
GE_overlap.txt
SNP_overlap.txt
overlap_summary.txt

Run
---
source scripts/00_config.sh
python3 scripts/eqtl/02_prepare_overlap/prepare_overlap_matrices.py
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Sequence


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def read_tsv_header_cols(path: Path) -> List[str]:
    with path.open() as f:
        header = f.readline()
    if not header:
        die(f"Empty file or missing header: {path}")
    return header.rstrip("\n").split("\t")


def header_sample_ids(path: Path) -> List[str]:
    cols = read_tsv_header_cols(path)
    if len(cols) < 2:
        die(f"Header has <2 columns: {path}")
    # cols[0] = geneid/snpid; remaining = sample IDs
    return [c.rstrip("\r") for c in cols[1:]]


def read_meta_sample_order(meta_csv: Path, sample_col: str = "sample_name") -> List[str]:
    """
    Read sample IDs from metadata CSV in file order.
    Requires a column named `sample_col`.
    """
    out: List[str] = []
    with meta_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            die(f"Metadata CSV has no header: {meta_csv}")
        if sample_col not in reader.fieldnames:
            die(
                f"Metadata CSV must contain column '{sample_col}'. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            s = (row.get(sample_col) or "").strip().strip('"').rstrip("\r")
            if s:
                out.append(s)

    if not out:
        die(f"No sample IDs found in metadata column '{sample_col}': {meta_csv}")
    return out


def build_overlap_order(
    ge_samples: Sequence[str],
    snp_samples: Sequence[str],
    meta_samples: Sequence[str],
    order_source: str = "meta",
) -> List[str]:
    """
    Compute (GE ∩ SNP ∩ META) but preserve a chosen order.
    order_source:
      - 'meta' : preserve metadata file order (recommended)
      - 'ge'   : preserve GE header order
    """
    ge_set = set(ge_samples)
    snp_set = set(snp_samples)

    ge_snp = ge_set.intersection(snp_set)
    if not ge_snp:
        die("GE ∩ SNP overlap is empty (no shared sample IDs). Check headers.")

    if order_source == "meta":
        # keep meta order, but only those present in GE∩SNP
        final = [s for s in meta_samples if s in ge_snp]
    elif order_source == "ge":
        meta_set = set(meta_samples)
        final = [s for s in ge_samples if (s in ge_snp and s in meta_set)]
    else:
        die(f"Unknown order_source: {order_source}")

    if not final:
        die("Final overlap (GE ∩ SNP ∩ META) is empty. Check sample IDs & metadata column.")
    return final


def subset_matrix(
    in_path: Path,
    out_path: Path,
    keep_samples: List[str],
    label: str,
) -> None:
    """
    Stream-subset a Matrix-eQTL matrix:
      col1 = id
      col2.. = samples
    Preserves keep_samples order.
    Safety:
      - refuses if 0 kept samples
      - refuses if output would have no data rows
    """
    with in_path.open() as f:
        header = f.readline()
        if not header:
            die(f"{label} input is empty: {in_path}")

        cols = header.rstrip("\n").split("\t")
        if len(cols) < 2:
            die(f"{label} header has <2 columns: {in_path}")

        col_index: Dict[str, int] = {c.rstrip("\r"): i for i, c in enumerate(cols)}

        kept_idx = [0]  # keep row id col
        present = 0
        missing = 0
        missing_ids: List[str] = []

        for s in keep_samples:
            if s in col_index:
                kept_idx.append(col_index[s])
                present += 1
            else:
                missing += 1
                if len(missing_ids) < 10:
                    missing_ids.append(s)

        if present == 0:
            die(
                f"{label}: none of the overlap sample IDs were found in header.\n"
                f"  in_file: {in_path}\n"
                f"  first_missing_examples: {missing_ids}\n"
            )

        out_path.parent.mkdir(parents=True, exist_ok=True)

        n_rows = 0
        with out_path.open("w") as out:
            out.write("\t".join(cols[i] for i in kept_idx) + "\n")
            for line in f:
                if not line:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < len(cols):
                    parts += [""] * (len(cols) - len(parts))
                out.write("\t".join(parts[i] for i in kept_idx) + "\n")
                n_rows += 1

        if n_rows == 0:
            die(f"{label}: wrote header but no data rows to {out_path}")

        print(
            f"{label}: wrote {out_path} "
            f"(kept_samples={present}, missing_in_header={missing}, rows={n_rows})"
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--ge", default=None, help="Path to GE.txt")
    p.add_argument("--snp", default=None, help="Path to SNP.txt")
    p.add_argument("--meta", default=None, help="Path to clinical_data.csv")
    p.add_argument("--meta-sample-col", default="sample_name", help="Metadata sample column name")
    p.add_argument(
        "--order",
        choices=["meta", "ge"],
        default="meta",
        help="Which file's order to preserve for the final overlap list",
    )
    p.add_argument("--outdir", default=None, help="Output directory (eqtl dir)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    repo_root = os.environ.get("REPO_ROOT", "")
    if not repo_root:
        die("REPO_ROOT not set. Run: source scripts/00_config.sh")

    base_eqtl = Path(repo_root) / "output/genotype_run1/eqtl"

    ge_path = Path(args.ge) if args.ge else (base_eqtl / "GE.txt")
    snp_path = Path(args.snp) if args.snp else (base_eqtl / "SNP.txt")
    meta_path = Path(args.meta) if args.meta else (Path(repo_root) / "metadata/clinical_data.csv")
    outdir = Path(args.outdir) if args.outdir else base_eqtl

    for pth in (ge_path, snp_path, meta_path):
        if not pth.exists():
            die(f"Missing file: {pth}")

    outdir.mkdir(parents=True, exist_ok=True)

    # --- read sample IDs ---
    ge_samples = header_sample_ids(ge_path)
    snp_samples = header_sample_ids(snp_path)
    meta_samples = read_meta_sample_order(meta_path, sample_col=args.meta_sample_col)

    ge_set = set(ge_samples)
    snp_set = set(snp_samples)
    meta_set = set(meta_samples)

    ge_n = len(ge_set)
    snp_n = len(snp_set)
    ge_snp_n = len(ge_set.intersection(snp_set))
    triple = build_overlap_order(ge_samples, snp_samples, meta_samples, order_source=args.order)
    triple_n = len(triple)

    # --- write overlap list (preserving chosen order) ---
    overlap_path = outdir / "uasg_ge_snp_meta_overlap.txt"
    overlap_path.write_text("\n".join(triple) + "\n")

    # --- subset matrices ---
    ge_out = outdir / "GE_overlap.txt"
    snp_out = outdir / "SNP_overlap.txt"
    subset_matrix(ge_path, ge_out, triple, "GE")
    subset_matrix(snp_path, snp_out, triple, "SNP")

    # --- write summary ---
    summary_path = outdir / "overlap_summary.txt"
    summary = (
        f"GE_UASG_cols:          {ge_n}\n"
        f"SNP_UASG_cols:         {snp_n}\n"
        f"GE_SNP_overlap:        {ge_snp_n}\n"
        f"GE_SNP_META_overlap:   {triple_n}\n\n"
        f"Order preserved from:   {args.order}\n"
        f"Wrote overlap list:\n  {overlap_path}\n"
        f"Wrote matrices:\n  {ge_out}\n  {snp_out}\n"
    )
    summary_path.write_text(summary)

    print(summary.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())