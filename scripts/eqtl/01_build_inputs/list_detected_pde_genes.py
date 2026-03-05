#!/usr/bin/env python3
"""
List PDE genes detected in human expression dataset (GE.txt).

Uses authoritative GE.txt generated earlier in pipeline.

Output:
  output/eqtl/results/inspect/pde/pde_genes_detected_in_expression.txt
"""

import os
from pathlib import Path
import pandas as pd

REPO_ROOT = os.environ.get("REPO_ROOT", "")
if not REPO_ROOT:
    raise SystemExit("REPO_ROOT not set. Run: source scripts/00_config.sh")

ge_path = Path(REPO_ROOT) / "output/eqtl/GE.txt"

if not ge_path.exists():
    raise SystemExit(f"Missing GE.txt: {ge_path}")

# human PDE 
pde_genes = {
"PDE10A","PDE11A","PDE12","PDE1A","PDE1B","PDE1C",
"PDE2A","PDE3A","PDE3B","PDE4A","PDE4B","PDE4C",
"PDE4D","PDE4DIP","PDE5A","PDE6A","PDE6B","PDE6C",
"PDE6D","PDE6G","PDE6H","PDE7A","PDE7B","PDE8A",
"PDE8B","PDE9A"
}

ge = pd.read_csv(ge_path, sep="\t")

if "geneid" not in ge.columns:
    raise SystemExit("GE.txt must contain column 'geneid'")

detected = sorted(set(ge["geneid"]).intersection(pde_genes))

out_dir = Path(REPO_ROOT) / "output/eqtl/results/inspect/pde"
out_dir.mkdir(parents=True, exist_ok=True)

out_path = out_dir / "pde_genes_detected_in_expression.txt"

with open(out_path, "w") as f:
    for g in detected:
        f.write(g + "\n")

print("Detected PDE genes:", len(detected))
print("\n".join(detected))
print(f"\nSaved: {out_path}")

# run
# source scripts/00_config.sh
# python3 scripts/eqtl/01_build_inputs/list_detected_pde_genes.py