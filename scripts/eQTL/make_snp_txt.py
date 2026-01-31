#!/usr/bin/env python3
"""
make_snp_txt.py

Convert a VCF-like TSV (target_gene_report.tsv) into a Matrix-eQTL style SNP matrix.

Default:
  - Write ALL SNP rows to: output/eqtl/SNP.txt

Optional:
  - --subset N
    Write ONLY first N SNP rows to: output/eqtl/SNP.subsetN.txt
    (does NOT create SNP.txt)

Output format:
  - First column: snpid (from TSV "ID")
  - Header: sample IDs (UASG_#### ...)
  - Values: ALT allele count (0/1/2) from GT
"""

import os
import sys


# ---------------- helpers ---------------- #

def alt_count_from_gt(gt: str) -> str:
    """Convert GT (0|0, 0/1, 1|1, etc.) to ALT allele count (0/1/2)."""
    gt = (gt or "").strip()
    if gt in {".", "./.", ".|.", ""}:
        return "NA"

    sep = "|" if "|" in gt else ("/" if "/" in gt else None)
    if sep is None:
        return "NA"

    a, b = gt.split(sep, 1)
    if a == "." or b == ".":
        return "NA"

    try:
        ia, ib = int(a), int(b)
    except ValueError:
        return "NA"

    return str((ia != 0) + (ib != 0))


def parse_args(argv):
    """
    Usage:
      make_snp_txt.py [input.tsv] [--subset N]
    """
    in_path = None
    subset_n = None

    i = 0
    while i < len(argv):
        if argv[i] == "--subset":
            if i + 1 >= len(argv):
                sys.exit("ERROR: --subset requires an integer")
            subset_n = int(argv[i + 1])
            if subset_n <= 0:
                sys.exit("ERROR: --subset must be > 0")
            i += 2
        else:
            if in_path is not None:
                sys.exit("ERROR: Too many positional arguments")
            in_path = argv[i]
            i += 1

    return in_path, subset_n


def subset_out_path(full_out: str, n: int) -> str:
    """SNP.txt → SNP.subset100.txt"""
    d = os.path.dirname(full_out)
    base = os.path.basename(full_out)
    root, ext = os.path.splitext(base)
    return os.path.join(d, f"{root}.subset{n}{ext or '.txt'}")


def write_matrix(in_path: str, out_path: str, max_rows=None) -> int:
    """Write SNP matrix; optionally stop after max_rows."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    written = 0

    with open(in_path) as fin, open(out_path, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")

        id_idx = header.index("ID")
        fmt_idx = header.index("FORMAT")
        sample_ids = header[fmt_idx + 1:]

        fout.write("snpid\t" + "\t".join(sample_ids) + "\n")

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            snpid = fields[id_idx]

            fmt_keys = fields[fmt_idx].split(":")
            gt_i = fmt_keys.index("GT") if "GT" in fmt_keys else None

            row = [snpid]
            for cell in fields[fmt_idx + 1:fmt_idx + 1 + len(sample_ids)]:
                parts = cell.split(":")
                gt = parts[gt_i] if gt_i is not None and gt_i < len(parts) else "."
                row.append(alt_count_from_gt(gt))

            fout.write("\t".join(row) + "\n")
            written += 1

            if max_rows is not None and written >= max_rows:
                break

    return written


# ---------------- main ---------------- #

def main():
    REPO_ROOT = os.environ.get("REPO_ROOT")
    if not REPO_ROOT:
        sys.exit("ERROR: REPO_ROOT not set. Run: source scripts/00_config.sh")

    DEFAULT_IN = os.path.join(
        REPO_ROOT,
        "output/genotype_run1/target_lists/target_gene_report.tsv",
    )
    DEFAULT_OUT = os.path.join(
        REPO_ROOT,
        "output/eqtl/SNP.txt",
    )

    in_arg, subset_n = parse_args(sys.argv[1:])
    in_path = in_arg or DEFAULT_IN

    if subset_n is None:
        out_path = DEFAULT_OUT
        max_rows = None
        mode = "full"
    else:
        out_path = subset_out_path(DEFAULT_OUT, subset_n)
        max_rows = subset_n
        mode = f"subset ({subset_n} rows)"

    print(f"[INFO] Input : {in_path}")
    print(f"[INFO] Output: {out_path}")
    print(f"[INFO] Mode  : {mode}")

    n = write_matrix(in_path, out_path, max_rows)
    print(f"[INFO] Done  : wrote {n} SNP rows")


if __name__ == "__main__":
    main()

"""
========================
How to run
========================

From the repo root (illumina-genotyping-pipeline/):

1) Load environment variables (sets REPO_ROOT)
   source scripts/00_config.sh

2) Create full SNP matrix (ALL rows)
   python3 scripts/eQTL/make_snp_txt.py

   Output:
     output/eqtl/SNP.txt

3) Create subset only (FIRST N rows, e.g. 100)
   python3 scripts/eQTL/make_snp_txt.py --subset 100

   Output:
     output/eqtl/SNP.subset100.txt
"""