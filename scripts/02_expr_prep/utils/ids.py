# scripts/expr_prep/utils/ids.py
"""
Small helpers for handling sample identifiers.

These functions keep ID handling consistent across steps.
"""

from __future__ import annotations

import pandas as pd


def norm_ids(x) -> pd.Series:
    """
    Normalize IDs by converting to string and stripping whitespace.
    """
    return pd.Series(x, dtype="string").str.strip()