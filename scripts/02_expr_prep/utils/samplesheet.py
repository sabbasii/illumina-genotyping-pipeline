# scripts/expr_prep/utils/samplesheet.py
"""
Helpers for working with Illumina SampleSheet files.

These functions standardize how SampleSheets are read and how
sample identifiers are derived, so all steps behave the same way.
"""

from __future__ import annotations

import pandas as pd

from .io import read_auto
from .ids import norm_ids


def read_samplesheet(path: str, *, skiprows: int = 8) -> pd.DataFrame:
    """
    Read a SampleSheet CSV and clean column names.

    Most Illumina SampleSheets have a short header section;
    skiprows should usually be 8.
    """
    df = read_auto(path, skiprows=skiprows, dtype=str)
    df.columns = df.columns.astype(str).str.strip()
    return df


def add_iid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add IID and UASG columns derived from SampleSheet fields.

    IID format:
        <SentrixBarcode_A>_<SentrixPosition_A>
    """
    required = ["Sample_Name", "SentrixBarcode_A", "SentrixPosition_A"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"SampleSheet missing columns: {missing}")

    out = df.copy()
    out["Sample_Name"] = norm_ids(out["Sample_Name"])
    out["SentrixBarcode_A"] = norm_ids(out["SentrixBarcode_A"])
    out["SentrixPosition_A"] = norm_ids(out["SentrixPosition_A"])

    out["IID"] = out["SentrixBarcode_A"] + "_" + out["SentrixPosition_A"]
    out["UASG"] = out["Sample_Name"]

    return out


def uasg_to_iid_map(df: pd.DataFrame) -> dict:
    """
    Build a one-to-one mapping from UASG to IID.

    If a UASG appears more than once, the first occurrence is kept.
    """
    df2 = add_iid(df)
    df2 = (
        df2.dropna(subset=["UASG", "IID"])
           .drop_duplicates(subset=["UASG"], keep="first")
    )
    return df2.set_index("UASG")["IID"].to_dict()