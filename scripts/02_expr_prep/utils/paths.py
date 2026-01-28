# scripts/expr_prep/utils/paths.py
"""
Central place for resolving output paths used by the expression pipeline.

All directories are created here so individual step scripts do not
need to worry about filesystem setup.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from .log import kv



@dataclass(frozen=True)
class Paths:
    """
    Container for commonly used paths in expr_prep.
    """
    repo_root: Path
    run: str

    expr_out_dir: Path
    matrices: Path
    metadata: Path
    lists: Path
    reports: Path
    clinical: Path

    # Standard files written by the pipeline
    final_expr: Path
    filtered_expr: Path
    filtered_meta: Path
    filtered_ids: Path
    mrs_csv: Path
    report_txt: Path


def get_paths() -> Paths:
    """
    Resolve all paths from environment variables with safe defaults.
    """
    if "REPO_ROOT" not in os.environ:
        kv("Warning", "REPO_ROOT not set; using current working directory")
    repo_root = Path(os.environ.get("REPO_ROOT", os.getcwd())).resolve()
    run = os.environ.get("RUN", "genotype_run1")

    expr_out_dir = Path(
        os.environ.get(
            "EXPR_OUT_DIR",
            repo_root / "output" / run / "expr" / "prep",
        )
    ).resolve()

    matrices = Path(os.environ.get("EXPR_DIR_MATRICES", expr_out_dir / "matrices"))
    metadata = Path(os.environ.get("EXPR_DIR_METADATA", expr_out_dir / "metadata"))
    lists = Path(os.environ.get("EXPR_DIR_LISTS", expr_out_dir / "lists"))
    reports = Path(os.environ.get("EXPR_DIR_REPORTS", expr_out_dir / "reports"))
    clinical = Path(os.environ.get("EXPR_DIR_CLINICAL", expr_out_dir / "clinical"))

    # Ensure directories exist
    for d in (expr_out_dir, matrices, metadata, lists, reports, clinical):
        d.mkdir(parents=True, exist_ok=True)

    return Paths(
        repo_root=repo_root,
        run=run,
        expr_out_dir=expr_out_dir,
        matrices=matrices,
        metadata=metadata,
        lists=lists,
        reports=reports,
        clinical=clinical,
        final_expr=matrices / "expression_matrix_transposed_FINAL.csv",
        filtered_expr=matrices / "expression_matrix_filtered.csv",
        filtered_meta=metadata / "metadata_filtered.csv",
        filtered_ids=lists / "ids_filtered.txt",
        mrs_csv=clinical / "mrs_outcome_by_sample.csv",
        report_txt=reports / "expr_metadata_check_report.txt",
    )