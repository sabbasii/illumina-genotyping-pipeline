#!/usr/bin/env python3
import pandas as pd

IN_TSV = "/home/sima/git_projects/illumina-genotyping-pipeline/input_data/expr_array/expr_table.tsv"
OUT_TXT = "/home/sima/git_projects/illumina-genotyping-pipeline/output/genotype_run1/eqtl/GE.txt"

df = pd.read_csv(IN_TSV, sep="\t")

# --- sanity checks ---
required = {"SYMBOL"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

uasg_cols = [c for c in df.columns if c.startswith("UASG_")]
if not uasg_cols:
    raise ValueError("No columns starting with 'UASG_' were found.")

# --- build GE table ---
ge = df[["SYMBOL"] + uasg_cols].copy()
ge = ge.rename(columns={"SYMBOL": "geneid"})

# Drop rows without geneid
ge = ge.dropna(subset=["geneid"])
ge = ge[ge["geneid"].astype(str).str.strip() != ""]

# OPTIONAL (recommended): collapse duplicates by mean
# Comment this out if you truly want duplicates kept
ge = ge.groupby("geneid", as_index=False)[uasg_cols].mean()

# Write GE.txt as TSV
ge.to_csv(OUT_TXT, sep="\t", index=False)
print(f"Wrote {OUT_TXT} with shape {ge.shape}")
