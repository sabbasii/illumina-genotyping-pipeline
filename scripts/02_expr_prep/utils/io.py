# scripts/expr_prep/utils/io.py
"""
I/O helpers used across the expression prep steps.

Goals:
- Read CSV/TSV files without hardcoding the delimiter.
- Fail early with clear errors when files are missing.
- Keep behavior consistent across all step scripts.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import pandas as pd


def ensure_file(path: str | os.PathLike) -> str:
    """
    Return an absolute path and raise a clear error if the file is missing.
    """
    p = Path(path).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"Missing file: {p}")
    return str(p.resolve())


def read_auto(
    path: str | os.PathLike,
    *,
    skiprows: Optional[int] = None,
    dtype: Optional[object] = None,
    encoding: str = "utf-8",
) -> pd.DataFrame:
    """
    Read a delimited text file (CSV/TSV) and auto-detect the delimiter.

    Notes:
    - Uses pandas' Python parser so the delimiter can be inferred.
    - If your input is a SampleSheet with header lines, pass skiprows=8.
    - For identifiers (barcodes, positions, sample IDs), dtype=str is often safer.
    """
    path = ensure_file(path)
    return pd.read_csv(
        path,
        sep=None,            # infer delimiter
        engine="python",     # required for sep=None inference
        encoding=encoding,
        on_bad_lines="warn",
        skiprows=skiprows,
        dtype=dtype,
    )


def write_text(path: str | os.PathLike, text: str, *, encoding: str = "utf-8") -> None:
    """
    Write a text file and create parent folders if needed.
    """
    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding=encoding)


def write_lines(
    path: str | os.PathLike, lines: list[str], *, encoding: str = "utf-8"
) -> None:
    """
    Write one line per item (adds a final newline).
    """
    write_text(path, "\n".join(lines) + "\n", encoding=encoding)
