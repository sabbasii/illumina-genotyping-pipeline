"""
make_snp_and_snpsloc_from_bcf.py

Convert genotype data from a BCF file into Matrix eQTL input files.

What this script does
---------------------
• Streams variants from BCF using bcftools
• Converts GT → ALT allele counts (0 / 1 / 2)
• Normalizes chromosome names (1 → chr1)
• Creates stable SNP IDs (uses ID or chr:pos:REF:ALT)

Inputs
------
BCF:
    output/imputed/genotypes_processed.bcf

Outputs
-------
output/eqtl_allSNPs/

• SNP.txt
    Matrix eQTL SNP matrix (ALT allele counts)

• snpsloc.txt
    SNP genomic locations (snpid, chr, pos)

Run
---
source scripts/00_config.sh
python3 scripts/eqtl/01_build_inputs/make_snp_and_snpsloc_from_bcf.py
"""

#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path
from typing import List

# ---- make scripts/ importable so "import eqtl..." works ----
HERE = Path(__file__).resolve()
SCRIPTS_DIR = HERE.parents[2]  # .../scripts
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from eqtl.utils.genotype_helpers import alt_count_from_gt, normalize_chr, stable_snpid

# ---------- CONFIG (repo-relative) ----------
BCF_REL = "output/imputed/genotypes_processed.bcf"
OUTDIR_REL = "output/eqtl_allSNPs"
SNP_OUT = "SNP.txt"
SNPLOC_OUT = "snpsloc.txt"


def run(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def repo_root_from_this_script() -> Path:
    # scripts/eqtl/01_build_inputs/<thisfile> -> up 3 levels to repo root
    return HERE.parents[3]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def main() -> None:
    repo_root = repo_root_from_this_script()
    bcf = repo_root / BCF_REL
    outdir = repo_root / OUTDIR_REL
    ensure_dir(outdir)

    if not bcf.exists():
        die(f"BCF not found: {bcf}")

    # index check (optional)
    if not ((bcf.with_suffix(bcf.suffix + ".csi")).exists() or (bcf.with_suffix(bcf.suffix + ".tbi")).exists()):
        print("WARN: No .csi/.tbi index found next to BCF. Continuing anyway.", file=sys.stderr)

    # bcftools presence
    cp = run(["bcftools", "--version"])
    if cp.returncode != 0:
        die("bcftools not found in PATH (activate env, or install bcftools).")

    # sample IDs
    cp = run(["bcftools", "query", "-l", str(bcf)])
    if cp.returncode != 0:
        die(f"bcftools query -l failed:\n{cp.stderr}")
    samples = [s.strip() for s in cp.stdout.splitlines() if s.strip()]
    if not samples:
        die("No samples found in BCF (bcftools query -l returned empty).")

    snp_path = outdir / SNP_OUT
    snploc_path = outdir / SNPLOC_OUT

    # Query fields:
    # CHROM POS ID REF ALT [GT per sample]
    fmt = r"%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", fmt, str(bcf)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    written = 0
    with snp_path.open("w", encoding="utf-8") as f_snp, snploc_path.open("w", encoding="utf-8") as f_loc:
        # headers
        f_snp.write("snpid\t" + "\t".join(samples) + "\n")
        f_loc.write("snpid\tchr\tpos\n")

        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            chrom_raw, pos, vid, ref, alt = parts[:5]
            gts = parts[5:]

            chrom = normalize_chr(chrom_raw)
            snpid = stable_snpid(vid, chrom, pos, ref, alt)

            # snpsloc row
            f_loc.write(f"{snpid}\t{chrom}\t{pos}\n")

            # SNP row (ALT counts)
            vals = [alt_count_from_gt(gt) for gt in gts]
            f_snp.write(snpid + "\t" + "\t".join(vals) + "\n")

            written += 1
            if written % 100000 == 0:
                print(f"...processed {written:,} variants", file=sys.stderr)

    stderr = proc.stderr.read() if proc.stderr is not None else ""
    ret = proc.wait()
    if ret != 0:
        die(f"bcftools query failed (exit {ret}).\n{stderr}")

    print("Done.")
    print(f"BCF:   {bcf}")
    print(f"Wrote: {snp_path}")
    print(f"Wrote: {snploc_path}")
    print("SNP.txt contains ALT-allele counts from GT (0/1/2); missing = NA.")


if __name__ == "__main__":
    main()


# Note
# I will later add --outdir and --bcf CLI args (so it can reuse for other runs without editing constants)