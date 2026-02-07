#!/usr/bin/env python3
"""
matrix_helpers.py

Helpers for Matrix-eQTL matrix files (GE.txt, SNP.txt, Covariates.txt).

- read_tsv_header: read header columns from a TSV
- header_sample_ids: get sample IDs from matrix header (skipping first column)
- read_overlap_list: read 1-column overlap list (preserve order)
- intersect_preserve_order: keep items that exist in a set, preserving original order
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence


def read_tsv_header(path: Path) -> List[str]:
    with path.open() as f:
        header = f.readline()
    if not header:
        raise SystemExit(f"ERROR: empty file or missing header: {path}")
    return header.rstrip("\n").split("\t")


def header_sample_ids(path: Path) -> List[str]:
    """
    For Matrix-eQTL style matrices:
      col1 = id (geneid/snpid/id)
      col2.. = sample IDs
    Returns sample IDs in file header order.
    """
    cols = read_tsv_header(path)
    if len(cols) < 2:
        raise SystemExit(f"ERROR: header has <2 columns: {path}")
    # strip CR if any
    return [c.rstrip("\r") for c in cols[1:]]


def read_overlap_list(path: Path) -> List[str]:
    keep: List[str] = []
    with path.open() as f:
        for line in f:
            s = line.strip().strip('"').rstrip("\r")
            if s:
                keep.append(s)
    if not keep:
        raise SystemExit(f"ERROR: overlap list is empty: {path}")
    return keep


def intersect_preserve_order(items: Sequence[str], allowed: Iterable[str]) -> List[str]:
    allowed_set = set(allowed)
    return [x for x in items if x in allowed_set]
