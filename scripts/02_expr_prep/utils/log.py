# scripts/expr_prep/utils/log.py
"""
Small logging helpers for pipeline steps.

These functions keep console output and text reports
simple, consistent, and easy to scan.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterable


def timestamp() -> str:
    """
    Return the current time in a compact, readable format.
    """
    return datetime.now().isoformat(timespec="seconds")


def banner(title: str) -> None:
    """
    Print a clear section header to the console.
    """
    print(f"\n=== {title} ===")


def kv(label: str, value) -> None:
    """
    Print a simple key–value line.
    """
    print(f"{label:<24}: {value}")


def write_report(path: str | Path, lines: Iterable[str]) -> None:
    """
    Write a plain text report file.

    Each item in `lines` is written on its own line.
    Parent folders are created if needed.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
