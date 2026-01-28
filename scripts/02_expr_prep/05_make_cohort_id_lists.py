#!/usr/bin/env python3
"""
05_make_cohort_id_lists.py

Map selected UASG IDs (from step 04) to SNP IIDs using the Illumina SampleSheet.

Writes:
- metadata/pheno/iid_filtered.keep          (IIDs only; one per line; deduplicated)
- metadata/pheno/mapping_uasg_to_iid.tsv    (UASG ↔ IID mapping; for tracing)
- metadata/vcf_samples_filtered.txt         (same IID list; convenient copy)

Inputs are taken ONLY from environment variables (set by scripts/00_config.sh):
  SAMPLE_SHEET : path to Illumina SampleSheet CSV
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

from .utils.config import load_config
from .utils.ids import norm_ids
from .utils.io import ensure_file, write_lines
from .utils.log import banner, kv
from .utils.paths import get_paths
from .utils.samplesheet import read_samplesheet, uasg_to_iid_map


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
    banner("05 — Make cohort ID lists (UASG → IID)")

    # ---- setup ----
    P = get_paths()
    cfg = load_config()
    SAMPLESHEET_SKIPROWS = int(cfg["SAMPLESHEET_SKIPROWS"])

    samplesheet_path = _require_env_file("SAMPLE_SHEET")
    kv("SampleSheet", samplesheet_path)

    # Input from step 04
    ids_filtered = P.filtered_ids
    if not ids_filtered.is_file():
        sys.exit(
            f"[ERR] ids_filtered.txt not found: {ids_filtered}\n"
            f"      Run 04 first:\n"
            f"        python3 -m scripts.expr_prep.04_expr_overlap_and_select"
        )

    # Output dirs (tracked)
    repo_root = P.repo_root
    pheno_dir = Path(os.environ.get("PHENO_DIR", str(repo_root / "metadata" / "pheno"))).resolve()
    meta_dir = Path(os.environ.get("META_DIR", str(repo_root / "metadata"))).resolve()

    pheno_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)

    iid_keep_out = pheno_dir / "iid_filtered.keep"
    u2i_tsv_out = pheno_dir / "mapping_uasg_to_iid.tsv"
    cohort_samples_out = meta_dir / "vcf_samples_filtered.txt"

    # ---- load selected UASGs ----
    raw_uasg: list[str] = ids_filtered.read_text(encoding="utf-8").splitlines()
    sel_uasg = [u for u in norm_ids(raw_uasg) if u]  # normalize to match mapping keys

    if not sel_uasg:
        sys.exit(f"[ERR] {ids_filtered} is empty after cleaning. Aborting.")

    # ---- load SampleSheet and build mapping ----
    ss = read_samplesheet(samplesheet_path, skiprows=SAMPLESHEET_SKIPROWS)
    u2i = uasg_to_iid_map(ss)  # keys are normalized UASG

    # ---- map and validate ----
    mapped: list[tuple[str, str]] = []
    missing_uasg: list[str] = []

    for u in sel_uasg:
        iid = u2i.get(u)
        if not iid:
            missing_uasg.append(u)
            continue
        mapped.append((u, iid))

    if not mapped:
        sys.exit("[ERR] No IIDs mapped. Aborting.")

    # Deduplicate IIDs for keep-list tools (plink/bcftools)
    mapped_iids_all = [iid for _, iid in mapped]
    mapped_iids_unique = sorted(set(mapped_iids_all))

    # Track duplicates (informational)
    iid_counts: dict[str, int] = {}
    for iid in mapped_iids_all:
        iid_counts[iid] = iid_counts.get(iid, 0) + 1
    dup_iids = sorted([iid for iid, c in iid_counts.items() if c > 1])

    # ---- write outputs ----
    write_lines(iid_keep_out, mapped_iids_unique)

    pd.DataFrame(mapped, columns=["UASG", "IID"]).to_csv(
        u2i_tsv_out, sep="\t", index=False
    )

    write_lines(cohort_samples_out, mapped_iids_unique)

    # ---- summary ----
    banner("05_make_cohort_id_lists summary")
    print(f"RUN                 : {P.run}")
    print(f"Selected UASGs in    : {ids_filtered}")
    print(f"SampleSheet          : {samplesheet_path}")
    print(f"Mapped rows (U→IID)  : {len(mapped)}")
    print(f"Unique IIDs written  : {len(mapped_iids_unique)}")
    print("")
    print(f"Wrote keep list      : {iid_keep_out}")
    print(f"Wrote mapping TSV    : {u2i_tsv_out}")
    print(f"Wrote cohort list    : {cohort_samples_out}")

    if missing_uasg:
        print(f"\n[WARN] UASGs with no IID in SampleSheet ({len(missing_uasg)}):")
        for u in missing_uasg[:10]:
            print("  -", u)
        if len(missing_uasg) > 10:
            print("  ...", len(missing_uasg) - 10, "more")

    if dup_iids:
        print(f"\n[WARN] Duplicate IIDs mapped from multiple UASGs ({len(dup_iids)}):")
        for iid in dup_iids[:10]:
            print("  -", iid, f"(count={iid_counts[iid]})")
        if len(dup_iids) > 10:
            print("  ...", len(dup_iids) - 10, "more")


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------
# How to run
#
# From the repository root:
#
#   source scripts/00_config.sh
#   python3 -m scripts.expr_prep.05_make_cohort_id_lists
#
# Inputs:
#   $EXPR_OUT_DIR/lists/ids_filtered.txt
#
# Outputs:
#   metadata/pheno/iid_filtered.keep
#   metadata/pheno/mapping_uasg_to_iid.tsv
#   metadata/vcf_samples_filtered.txt
# ---------------------------------------------------------------------