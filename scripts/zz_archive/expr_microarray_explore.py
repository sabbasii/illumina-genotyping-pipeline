#!/usr/bin/env python3
import os, sys, re, gzip, io, json
from datetime import datetime

import pandas as pd

# ---------- config from environment (00_config.sh sets RUN/REPO_ROOT/OUT_DIR/META_DIR) ----------
REPO_ROOT = os.environ.get("REPO_ROOT", os.getcwd())
RUN       = os.environ.get("RUN", "genotype_run1")
OUT_DIR   = os.environ.get("OUT_DIR", os.path.join(REPO_ROOT, "output", RUN))
META_DIR  = os.environ.get("META_DIR", os.path.join(REPO_ROOT, "metadata"))
EXP_INDIR = os.path.join(REPO_ROOT, "input_data", "expression_microarray")
EXP_OUT   = os.path.join(OUT_DIR, "expr", "explore")
os.makedirs(EXP_OUT, exist_ok=True)

FILES = {
    "tga_txt": os.path.join(EXP_INDIR, "TGA-cohort.txt"),
    "tga_csv": os.path.join(EXP_INDIR, "TGA-cohort-BackUp.csv"),
    "transpose_csv": os.path.join(EXP_INDIR, "transpose_numbers.csv"),
}

# ---------- helpers ----------
def sniff_and_read(path):
    """Best-effort read: detect encoding, delimiter; return (DataFrame, info_dict)."""
    info = {"path": path, "exists": os.path.exists(path), "rows": 0, "cols": 0, "delimiter": None, "encoding": None}
    if not info["exists"]:
        return None, info
    # read small head to sniff
    with open(path, "rb") as fh:
        head = fh.read(65536)
    # encoding guess (very light)
    for enc in ("utf-8-sig", "utf-8", "iso-8859-1", "cp1252"):
        try:
            head.decode(enc)
            info["encoding"] = enc
            break
        except UnicodeDecodeError:
            continue
    if info["encoding"] is None:
        info["encoding"] = "utf-8"

    # delimiter guess
    preview = head.decode(info["encoding"], errors="replace")
    candidates = [",", "\t", ";", "|"]
    counts = {d: preview.count(d) for d in candidates}
    delim = max(counts, key=counts.get)
    info["delimiter"] = delim if counts[delim] > 0 else ","

    df = pd.read_csv(path, encoding=info["encoding"], sep=info["delimiter"], engine="python")
    # normalize column names
    df.columns = [re.sub(r"[^0-9a-zA-Z]+", "_", c.strip()).strip("_").lower() for c in df.columns]
    info["rows"], info["cols"] = df.shape
    return df, info

def guess_id_col(cols):
    pats = [
        r"^iid$", r"^sample_?id$", r"^id$", r"^subject_?id$",
        r"^sentrix.*(barcode|position)", r"^gsm\d+$"
    ]
    for c in cols:
        lc = c.lower()
        if any(re.search(p, lc) for p in pats):
            return c
    # fallback: first column
    return cols[0] if cols else None

def guess_pheno_col(cols):
    pats = [r"^pheno", r"^phenotype", r"^status$", r"^case[_-]?control$", r"^group$", r"^diagnosis$", r"^stroke$"]
    for c in cols:
        lc = c.lower()
        if any(re.search(p, lc) for p in pats):
            return c
    return None

def tidy_preview(df, name, out_dir=EXP_OUT, n=10):
    out = os.path.join(out_dir, f"{name}.head.csv")
    df.head(n).to_csv(out, index=False)
    return out

def load_optional_psam_and_vcf_samples():
    paths = {
        "psam": os.path.join(META_DIR, "cohort.sex.psam"),
        "vcf_samples": os.path.join(META_DIR, "vcf.samples"),
    }
    psam = None
    if os.path.exists(paths["psam"]):
        psam = pd.read_csv(paths["psam"], sep="\t", comment="#")
        psam.columns = [c.lower() for c in psam.columns]
    vcf = None
    if os.path.exists(paths["vcf_samples"]):
        vcf = pd.read_csv(paths["vcf_samples"], header=None, names=["iid"])
    return psam, vcf, paths

# ---------- main ----------
report_lines = []
report_lines.append(f"# Microarray expression exploration ({RUN})")
report_lines.append(f"_Generated: {datetime.utcnow().isoformat()}Z_")
report_lines.append("")
report_lines.append(f"Input dir: `{EXP_INDIR}`  \nOutput dir: `{EXP_OUT}`")
report_lines.append("")

dfs = {}
infos = {}
for key, path in FILES.items():
    df, info = sniff_and_read(path)
    infos[key] = info
    if df is not None:
        dfs[key] = df
        tidy_path = tidy_preview(df, f"{key}")
        report_lines.append(f"**{key}**: {info['rows']} rows ?? {info['cols']} cols  (enc=`{info['encoding']}`, sep=`{info['delimiter']}`)  \nPreview ??? `{os.path.relpath(tidy_path, REPO_ROOT)}`")
    else:
        report_lines.append(f"**{key}**: _missing_ at `{path}`")
report_lines.append("")

# guesses & summaries
for name, df in dfs.items():
    id_col = guess_id_col(df.columns.tolist())
    ph_col = guess_pheno_col(df.columns.tolist())
    report_lines.append(f"### {name}")
    report_lines.append(f"- guessed ID column: `{id_col}`")
    report_lines.append(f"- guessed phenotype column: `{ph_col}`")
    # simple value counts
    if ph_col and ph_col in df.columns:
        vc = df[ph_col].astype(str).str.strip().str.lower().value_counts(dropna=False)
        vc_path = os.path.join(EXP_OUT, f"{name}.phenotype_value_counts.csv")
        vc.to_csv(vc_path, header=["count"])
        report_lines.append(f"- phenotype value counts ??? `{os.path.relpath(vc_path, REPO_ROOT)}`")
    # types
    dtypes_path = os.path.join(EXP_OUT, f"{name}.dtypes.csv")
    pd.Series(df.dtypes.astype(str)).to_csv(dtypes_path, header=["dtype"])
    report_lines.append(f"- dtypes ??? `{os.path.relpath(dtypes_path, REPO_ROOT)}`")
    report_lines.append("")

# cross-file ID reconciliation
def normalize_id_series(s):
    return s.astype(str).str.strip().str.replace(r"\s+", "", regex=True)

id_maps = {}
for name, df in dfs.items():
    id_col = guess_id_col(df.columns.tolist())
    if id_col and id_col in df.columns:
        ids = normalize_id_series(df[id_col]).rename("iid")
        id_maps[name] = ids

if id_maps:
    all_ids = pd.Series(dtype=str)
    for name, s in id_maps.items():
        tag = pd.Series([name]*len(s), index=s.index, dtype=str, name="source")
        tmp = pd.concat([s.rename("iid"), tag], axis=1)
        all_ids = pd.concat([all_ids, s], ignore_index=True)
    # pairwise overlaps
    names = list(id_maps.keys())
    venn_rows = []
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            a = set(id_maps[names[i]])
            b = set(id_maps[names[j]])
            venn_rows.append({"A": names[i], "B": names[j], "A_only": len(a-b), "B_only": len(b-a), "A???B": len(a&b)})
    if venn_rows:
        venn_df = pd.DataFrame(venn_rows)
        venn_path = os.path.join(EXP_OUT, "id_overlap_pairs.csv")
        venn_df.to_csv(venn_path, index=False)
        report_lines.append(f"**ID overlaps** ??? `{os.path.relpath(venn_path, REPO_ROOT)}`")
        report_lines.append("")

# compare to genotyping identifiers if available
psam, vcf, paths = load_optional_psam_and_vcf_samples()
if vcf is not None:
    vcf["iid_norm"] = normalize_id_series(vcf["iid"])
if psam is not None and "iid" in psam.columns:
    psam["iid_norm"] = normalize_id_series(psam["iid"])

for name, ids in id_maps.items():
    s = normalize_id_series(ids).to_frame()
    s.columns = ["iid_norm"]
    if vcf is not None:
        inter = pd.Series(list(set(s["iid_norm"]) & set(vcf["iid_norm"])))
        out = os.path.join(EXP_OUT, f"{name}.overlap_with_vcf_samples.txt")
        inter.sort_values().to_csv(out, index=False, header=False)
    if psam is not None:
        inter = pd.Series(list(set(s["iid_norm"]) & set(psam["iid_norm"])))
        out = os.path.join(EXP_OUT, f"{name}.overlap_with_psam.txt")
        inter.sort_values().to_csv(out, index=False, header=False)

# write a simple markdown report
report_path = os.path.join(EXP_OUT, "exploration_report.md")
with open(report_path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(report_lines))

print(f"[OK] Wrote exploration outputs to: {EXP_OUT}")
print(f" - report: {report_path}")
for f in sorted(os.listdir(EXP_OUT)):
    print("   ", f)
