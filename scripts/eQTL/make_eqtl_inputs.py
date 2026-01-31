#!/usr/bin/env python3
"""
make_eqtl_inputs.py

Convert a VCF-like TSV (target_gene_report.tsv) into Matrix-eQTL input files.

Default:
  - Write ALL SNP rows to:
      output/genotype_run1/eqtl/SNP.txt
      output/genotype_run1/eqtl/snpsloc.txt

Optional:
  - --subset N
    Write ONLY first N SNP rows to:
      output/genotype_run1/eqtl/SNP.subsetN.txt
      output/genotype_run1/eqtl/snpsloc.subsetN.txt

SNP.txt format:
  - First column: snpid (from TSV "ID")
  - Header: sample IDs (UASG_#### ...)
  - Values: ALT allele count (0/1/2) from GT

snpsloc.txt format:
  - snpid, chr, pos  (from TSV ID, #CHROM, POS)
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
      make_eqtl_inputs.py [input.tsv] [--subset N]
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
    """SNP.txt -> SNP.subset100.txt (or snpsloc.txt -> snpsloc.subset100.txt)."""
    d = os.path.dirname(full_out)
    base = os.path.basename(full_out)
    root, ext = os.path.splitext(base)
    return os.path.join(d, f"{root}.subset{n}{ext or '.txt'}")


def write_matrices_and_snpsloc(in_path: str, snp_out: str, snpsloc_out: str, max_rows=None) -> int:
    """
    Write:
      - SNP matrix (snp_out)
      - snpsloc (snpsloc_out)
    If max_rows is None: write all SNP rows.
    If max_rows is an int: write only first max_rows SNP rows.
    Returns number of SNP rows written.
    """
    if not os.path.isfile(in_path):
        raise FileNotFoundError(f"Input TSV not found: {in_path}")

    os.makedirs(os.path.dirname(snp_out), exist_ok=True)
    os.makedirs(os.path.dirname(snpsloc_out), exist_ok=True)

    written = 0

    with open(in_path, "r", encoding="utf-8") as fin, \
         open(snp_out, "w", encoding="utf-8") as fout_snp, \
         open(snpsloc_out, "w", encoding="utf-8") as fout_loc:

        header = fin.readline().rstrip("\n").split("\t")

        try:
            chrom_idx = header.index("#CHROM")
            pos_idx = header.index("POS")
            id_idx = header.index("ID")
            fmt_idx = header.index("FORMAT")
        except ValueError as e:
            raise SystemExit(f"ERROR: Missing required column in header: {e}")

        sample_ids = header[fmt_idx + 1:]
        if not sample_ids:
            raise SystemExit("ERROR: No sample columns found after FORMAT.")

        # Headers
        fout_snp.write("snpid\t" + "\t".join(sample_ids) + "\n")
        fout_loc.write("snpid\tchr\tpos\n")

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < fmt_idx + 2:
                continue

            chrom = fields[chrom_idx]
            pos = fields[pos_idx]
            snpid = fields[id_idx]

            fmt_keys = fields[fmt_idx].split(":")
            gt_i = fmt_keys.index("GT") if "GT" in fmt_keys else None

            # snpsloc row
            fout_loc.write(f"{snpid}\t{chrom}\t{pos}\n")

            # SNP matrix row
            row = [snpid]
            for cell in fields[fmt_idx + 1:fmt_idx + 1 + len(sample_ids)]:
                parts = cell.split(":")
                gt = parts[gt_i] if gt_i is not None and gt_i < len(parts) else "."
                row.append(alt_count_from_gt(gt))

            fout_snp.write("\t".join(row) + "\n")
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
    DEFAULT_SNP_OUT = os.path.join(
        REPO_ROOT,
        "output/genotype_run1/eqtl/SNP.txt",
    )
    DEFAULT_SNPLOC_OUT = os.path.join(
        REPO_ROOT,
        "output/genotype_run1/eqtl/snpsloc.txt",
    )

    in_arg, subset_n = parse_args(sys.argv[1:])
    in_path = in_arg or DEFAULT_IN

    if subset_n is None:
        snp_out = DEFAULT_SNP_OUT
        snpsloc_out = DEFAULT_SNPLOC_OUT
        max_rows = None
        mode = "full"
    else:
        snp_out = subset_out_path(DEFAULT_SNP_OUT, subset_n)
        snpsloc_out = subset_out_path(DEFAULT_SNPLOC_OUT, subset_n)
        max_rows = subset_n
        mode = f"subset ({subset_n} rows)"

    print(f"[INFO] Input    : {in_path}")
    print(f"[INFO] SNP out  : {snp_out}")
    print(f"[INFO] snpsloc  : {snpsloc_out}")
    print(f"[INFO] Mode     : {mode}")

    n = write_matrices_and_snpsloc(in_path, snp_out, snpsloc_out, max_rows)
    print(f"[INFO] Done     : wrote {n} SNP rows")


if __name__ == "__main__":
    main()

"""
========================
How to run
========================

From the repo root (illumina-genotyping-pipeline/):

1) Load env vars (sets REPO_ROOT)
   source scripts/00_config.sh

2) Full outputs (ALL rows)
   python3 scripts/eQTL/make_eqtl_inputs.py

   Outputs:
     output/genotype_run1/eqtl/SNP.txt
     output/genotype_run1/eqtl/snpsloc.txt

3) Subset-only outputs (FIRST N rows, e.g. 100)
   python3 scripts/eQTL/make_eqtl_inputs.py --subset 100

   Outputs:
     output/genotype_run1/eqtl/SNP.subset100.txt
     output/genotype_run1/eqtl/snpsloc.subset100.txt
"""
