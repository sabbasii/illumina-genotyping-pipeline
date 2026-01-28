# scripts/expr_prep/utils/config.py
"""
Lightweight configuration loader for the expression pipeline.

Order of precedence:
1) Built-in defaults
2) Optional config.json (next to expr_prep/)
3) Environment variables

This lets the pipeline run out of the box, while still allowing
small adjustments without editing the code.
"""

from __future__ import annotations

import json
import os
from pathlib import Path


# Reasonable defaults used by all steps unless overridden
DEFAULTS = {
    "FIRST_SAMPLE_COL_IDX": 7,
    "CSV_PROBE_START_NAME": "TC0100006437.hg.1",
    "KEEP_DIAG": ["Control", "Ischemic Stroke", "TIA"],
    "SAMPLESHEET_SKIPROWS": 8,
}


def load_config() -> dict:
    """
    Load configuration values and return them as a dictionary.

    Values are taken from defaults, then optionally overridden by
    config.json and environment variables.
    """
    cfg = dict(DEFAULTS)

    # Optional JSON config (kept local to expr_prep)
    config_path = Path(__file__).resolve().parents[1] / "config.json"
    if config_path.is_file():
        with open(config_path, "r", encoding="utf-8") as f:
            user_cfg = json.load(f)
        for k, v in user_cfg.items():
            if v is not None:
                cfg[k] = v

    # Environment variable overrides (useful on HPC / CI)
    if "FIRST_SAMPLE_COL_IDX" in os.environ:
        cfg["FIRST_SAMPLE_COL_IDX"] = int(os.environ["FIRST_SAMPLE_COL_IDX"])

    if "CSV_PROBE_START_NAME" in os.environ:
        cfg["CSV_PROBE_START_NAME"] = os.environ["CSV_PROBE_START_NAME"]

    if "SAMPLESHEET_SKIPROWS" in os.environ:
        cfg["SAMPLESHEET_SKIPROWS"] = int(os.environ["SAMPLESHEET_SKIPROWS"])

    if "KEEP_DIAG" in os.environ:
        cfg["KEEP_DIAG"] = [
            x.strip()
            for x in os.environ["KEEP_DIAG"].split(",")
            if x.strip()
        ]

    return cfg