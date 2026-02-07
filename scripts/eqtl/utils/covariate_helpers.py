#!/usr/bin/env python3
"""
covariate_helpers.py

Helpers for converting Matrix-eQTL Covariates.txt (mixed types) into numeric form.

Assumes pandas is available (you already use it in the pipeline).

Functions:
- encode_common_strings: sex + yes/no conversions
- one_hot_encode_row: one-hot encode a categorical covariate row
- sanitize_covariates_numeric: core transformation used by preprocess_covariates.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class OneHotResult:
    dummies: pd.DataFrame          # rows=dummy covariates, cols=samples
    ref_level: Optional[str]
    counts: Dict[str, int]


def _lower_trim(x: str) -> str:
    return str(x).strip().lower()


def encode_common_strings(vals: pd.Series) -> pd.Series:
    """
    Convert common clinical encodings to numeric-like strings:
      - yes/no -> 1/0
      - sex female/male -> 0/1
      - blanks -> NA (kept as NA/NaN)
    Returns a Series of strings / NA.
    """
    v = vals.astype("string").copy()
    v = v.str.strip()

    # blanks -> <NA>
    v = v.replace("", pd.NA)

    low = v.str.lower()

    # yes/no
    v = v.mask(low.isin(["yes", "y", "true", "t"]), "1")
    v = v.mask(low.isin(["no", "n", "false", "f"]), "0")

    # sex
    v = v.mask(low.isin(["female", "f"]), "0")
    v = v.mask(low.isin(["male", "m"]), "1")

    return v


def one_hot_encode_row(
    vals: pd.Series,
    prefix: str,
    ref_level: Optional[str] = None,
) -> OneHotResult:
    """
    One-hot encode a categorical covariate row (vals across samples).
    - reference level is dropped (to avoid collinearity)
    - if ref_level not provided or missing, uses most frequent level
    """
    v = vals.astype("string").map(_lower_trim)
    v = v.replace("", pd.NA)

    counts = v.value_counts(dropna=True).to_dict()
    if not counts:
        return OneHotResult(dummies=pd.DataFrame(), ref_level=None, counts={})

    if ref_level is None or _lower_trim(ref_level) not in counts:
        ref = next(iter(counts.keys()))  # most frequent
    else:
        ref = _lower_trim(ref_level)

    levels = [lv for lv in counts.keys() if lv != ref]
    d = {}
    for lv in levels:
        safe = "".join(ch if ch.isalnum() else "_" for ch in lv)
        name = f"{prefix}_{safe}"
        d[name] = (v == lv).astype("float")  # 1.0/0.0; NA->False->0.0
        d[name] = d[name].where(~v.isna(), pd.NA)

    dummies = pd.DataFrame(d).T  # rows=dummies, cols=samples
    dummies.columns = vals.index

    return OneHotResult(dummies=dummies, ref_level=ref, counts={k: int(v) for k, v in counts.items()})


def sanitize_covariates_numeric(
    cov: pd.DataFrame,
    include_diagnosis: bool = True,
    include_ancestry: bool = True,
    diagnosis_ref: str = "control",
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Core conversion used by preprocess_covariates.py.

    Input cov:
      - index = covariate names (from 'id' column)
      - columns = sample IDs
      - values = strings/numbers

    Output:
      - numeric covariate matrix (index=covariates, columns=samples)
      - summary lines (list of strings)
    """
    summary: List[str] = []
    samples = list(cov.columns)

    def pop_row(name: str) -> Optional[pd.Series]:
        if name in cov.index:
            s = cov.loc[name].copy()
            cov.drop(index=name, inplace=True)
            return s
        return None

    diagnosis = pop_row("diagnosis") if include_diagnosis else None
    ancestry = pop_row("ancestry") if include_ancestry else None
    _ = pop_row("subtype")  # ignored by default (matches your R script)

    # encode remaining rows and coerce numeric
    cov2 = cov.copy()
    for rn in cov2.index:
        cov2.loc[rn] = encode_common_strings(cov2.loc[rn])

    cov_num = cov2.apply(pd.to_numeric, errors="coerce")

    # drop all-NA rows
    all_na = cov_num.isna().all(axis=1)
    dropped = cov_num.index[all_na].tolist()
    if dropped:
        summary.append("Dropped non-numeric covariate rows (all NA after conversion):")
        summary.extend([f"  - {x}" for x in dropped])
        summary.append("")
        cov_num = cov_num.loc[~all_na]

    # add diagnosis dummies
    if diagnosis is not None:
        res = one_hot_encode_row(diagnosis, prefix="dx", ref_level=diagnosis_ref)
        if res.counts:
            summary.append("Diagnosis one-hot encoding:")
            summary.append("  Levels (counts): " + ", ".join([f"{k}={v}" for k, v in res.counts.items()]))
            summary.append(f"  Reference dropped: {res.ref_level}")
            summary.append("")
        if not res.dummies.empty:
            cov_num = pd.concat([cov_num, res.dummies.apply(pd.to_numeric, errors="coerce")], axis=0)
        else:
            summary.append("  (No diagnosis dummies created; diagnosis missing or single-level.)")
            summary.append("")

    # add ancestry dummies
    if ancestry is not None:
        res = one_hot_encode_row(ancestry, prefix="anc", ref_level=None)
        if res.counts:
            summary.append("Ancestry one-hot encoding:")
            summary.append("  Levels (counts): " + ", ".join([f"{k}={v}" for k, v in res.counts.items()]))
            summary.append(f"  Reference dropped: {res.ref_level}")
            summary.append("")
        if not res.dummies.empty:
            cov_num = pd.concat([cov_num, res.dummies.apply(pd.to_numeric, errors="coerce")], axis=0)
        else:
            summary.append("  (No ancestry dummies created; ancestry missing or single-level.)")
            summary.append("")

    # preserve column order
    cov_num = cov_num[samples]
    return cov_num, summary
