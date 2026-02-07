#!/usr/bin/env python3
"""
make_snp_and_snpsloc_from_tsv.py

Convert a VCF-like TSV (e.g. target_gene_report.tsv) into Matrix-eQTL inputs:
  - SNP.txt     (ALT allele counts 0/1/2)
  - snpsloc.txt (snpid, chr, pos) with chr prefix (chr1, chr2, ...)

Defaults (REPO_ROOT required):
  Input:
    output/genotype_run1/target_lists/target_gene_report.tsv
  Outputs:
    output/genotype_run1/eqtl/SNP.txt
    output/genotype_run1/eqtl/snpsloc.txt

Optional:
  --subset N
    Writes only first N variant rows to:
      SNP.subsetN.txt
      snpsloc.subsetN.txt
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# NOTE: requires scripts/eqtl to be importable
HERE = Path(__file__).resolve()
SCRIPTS_DIR = HERE.parents[2]  # .../scripts
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
    
# Run from repo root with: source scripts/00_config.sh
from eqtl.utils.genotype_helpers import alt_count_from_gt, normalize_chr, stable_snpid


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def subset_out_path(full_out: Path, n: int) -> Path:
    """SNP.txt -> SNP.subset100.txt (or snpsloc.txt -> snpsloc.subset100.txt)."""
    root = full_out.stem
    ext = full_out.suffix or ".txt"
    return full_out.with_name(f"{root}.subset{n}{ext}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("input_tsv", nargs="?", default=None, help="Input TSV (VCF-like).")
    p.add_argument("--subset", type=int, default=None, help="Write only first N variant rows.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    repo_root = os.environ.get("REPO_ROOT", "")
    if not repo_root:
        die("REPO_ROOT not set. Run: source scripts/00_config.sh")

    default_in = Path(repo_root) / "output/genotype_run1/target_lists/target_gene_report.tsv"
    out_dir = Path(repo_root) / "output/genotype_run1/eqtl"
    out_dir.mkdir(parents=True, exist_ok=True)

    snp_out_full = out_dir / "SNP.txt"
    loc_out_full = out_dir / "snpsloc.txt"

    in_path = Path(args.input_tsv) if args.input_tsv else default_in
    if not in_path.exists():
        die(f"Input TSV not found: {in_path}")

    if args.subset is not None and args.subset <= 0:
        die("--subset must be > 0")

    if args.subset is None:
        snp_out = snp_out_full
        loc_out = loc_out_full
        max_rows = None
        mode = "full"
    else:
        snp_out = subset_out_path(snp_out_full, args.subset)
        loc_out = subset_out_path(loc_out_full, args.subset)
        max_rows = args.subset
        mode = f"subset ({args.subset} rows)"

    print(f"[INFO] Input   : {in_path}")
    print(f"[INFO] SNP out : {snp_out}")
    print(f"[INFO] LOC out : {loc_out}")
    print(f"[INFO] Mode    : {mode}")

    written = 0

    with in_path.open("r", encoding="utf-8") as fin, \
         snp_out.open("w", encoding="utf-8") as fout_snp, \
         loc_out.open("w", encoding="utf-8") as fout_loc:

        header_line = fin.readline()
        if not header_line:
            die(f"Empty TSV: {in_path}")

        header = header_line.rstrip("\n").split("\t")

        # Required columns
        required = ["#CHROM", "POS", "ID", "FORMAT"]
        missing = [c for c in required if c not in header]
        if missing:
            die(f"Missing required columns in TSV header: {missing}")

        chrom_idx = header.index("#CHROM")
        pos_idx = header.index("POS")
        id_idx = header.index("ID")
        fmt_idx = header.index("FORMAT")

        # Optional but recommended for stable snpid fallback
        ref_idx = header.index("REF") if "REF" in header else None
        alt_idx = header.index("ALT") if "ALT" in header else None

        sample_ids = header[fmt_idx + 1:]
        if not sample_ids:
            die("No sample columns found after FORMAT.")

        # Headers
        fout_snp.write("snpid\t" + "\t".join(sample_ids) + "\n")
        fout_loc.write("snpid\tchr\tpos\n")

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= fmt_idx:
                continue

            chrom_raw = fields[chrom_idx]
            pos = fields[pos_idx]
            vid = fields[id_idx]

            ref = fields[ref_idx] if ref_idx is not None and ref_idx < len(fields) else ""
            alt = fields[alt_idx] if alt_idx is not None and alt_idx < len(fields) else ""

            chrom = normalize_chr(chrom_raw)
            snpid = stable_snpid(vid, chrom, pos, ref, alt)

            fmt_keys = fields[fmt_idx].split(":")
            gt_i = fmt_keys.index("GT") if "GT" in fmt_keys else None

            # snpsloc row
            fout_loc.write(f"{snpid}\t{chrom}\t{pos}\n")

            # SNP matrix row
            row = [snpid]
            sample_cells = fields[fmt_idx + 1: fmt_idx + 1 + len(sample_ids)]
            if len(sample_cells) < len(sample_ids):
                # pad if line is short
                sample_cells += ["."] * (len(sample_ids) - len(sample_cells))

            for cell in sample_cells:
                parts = cell.split(":")
                gt = parts[gt_i] if gt_i is not None and gt_i < len(parts) else "."
                row.append(alt_count_from_gt(gt))

            fout_snp.write("\t".join(row) + "\n")
            written += 1

            if max_rows is not None and written >= max_rows:
                break

    print(f"[INFO] Done: wrote {written} variant rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


"""
========================
How to run
========================

From the repo root (illumina-genotyping-pipeline/):

1) Load env vars (sets REPO_ROOT)
   source scripts/00_config.sh

2) Full outputs (ALL rows)
   python3 scripts/eqtl/01_build_inputs/make_snp_and_snpsloc_from_tsv.py

   Outputs:
     output/genotype_run1/eqtl/SNP.txt
     output/genotype_run1/eqtl/snpsloc.txt

3) Subset-only outputs (FIRST N rows, e.g. 100)
   python3 scripts/eqtl/01_build_inputs/make_snp_and_snpsloc_from_tsv.py --subset 100

   Outputs:
     output/genotype_run1/eqtl/SNP.subset100.txt
     output/genotype_run1/eqtl/snpsloc.subset100.txt
"""
