#!/usr/bin/env python3
"""
genotype_helpers.py

Reusable helpers for building Matrix-eQTL SNP inputs.

- normalize_chr: ensure 'chr' prefix (1 -> chr1)
- alt_count_from_gt: convert VCF GT to ALT-allele count (0/1/2)
- stable_snpid: robust snpid (use ID if present else chr:pos:ref:alt)
"""

from __future__ import annotations

from typing import Optional

MISSING = "NA"


def normalize_chr(chrom: str) -> str:
    """Ensure chromosome has 'chr' prefix (1 -> chr1, chrX -> chrX)."""
    chrom = (chrom or "").strip()
    if not chrom:
        return chrom
    if chrom.lower().startswith("chr"):
        # keep original casing except normalize to 'chr' prefix style
        # (e.g., "CHR1" -> "chr1")
        return "chr" + chrom[3:]
    return f"chr{chrom}"


def alt_count_from_gt(gt: str, missing: str = MISSING) -> str:
    """
    Convert VCF GT to ALT-allele count:
      0/0 -> 0
      0/1 -> 1
      1/1 -> 2
    Handles phased (|) and missing (./. or .|.)
    Multi-allelic: counts any non-zero allele as ALT (0/2->1; 2/2->2).
    """
    gt = (gt or "").strip()
    if gt in {"", ".", "./.", ".|."}:
        return missing

    sep = "|" if "|" in gt else ("/" if "/" in gt else None)
    if sep is None:
        return missing

    a, b = gt.split(sep, 1)
    if a == "." or b == ".":
        return missing

    try:
        ia, ib = int(a), int(b)
    except ValueError:
        return missing

    return str((1 if ia != 0 else 0) + (1 if ib != 0 else 0))


def stable_snpid(
    vid: str,
    chrom: str,
    pos: str,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
) -> str:
    """
    Prefer rsID/ID if present, else fall back to CHR:POS:REF:ALT.
    chrom should already be normalized if you want chr-prefix.
    """
    vid = (vid or "").strip()
    if vid and vid != ".":
        return vid

    ref = (ref or "").strip()
    alt = (alt or "").strip()
    chrom = (chrom or "").strip()
    pos = (pos or "").strip()

    # If ref/alt unavailable, still produce a stable fallback
    if ref and alt:
        return f"{chrom}:{pos}:{ref}:{alt}"
    return f"{chrom}:{pos}"
