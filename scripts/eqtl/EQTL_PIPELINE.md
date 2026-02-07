# eQTL Pipeline — Illumina Genotyping + Expression Integration

This document explains how to run the eQTL pipeline step-by-step using the scripts under:

scripts/eqtl/


The pipeline performs:

- SNP matrix generation (from BCF or TSV)
- Gene expression matrix preparation
- Sample overlap harmonization
- Covariate construction + encoding
- Matrix eQTL association testing (ALL + CIS + TRANS)

---

# 📂 Pipeline Structure

```text
scripts/
├── 00_config.sh
│
└── eqtl/
    ├── 01_build_inputs/
    │   ├── build_GE_from_expr.py
    │   ├── make_geneloc_from_gtf.R
    │   ├── make_snp_and_snpsloc_from_tsv.py
    │   └── make_snp_and_snpsloc_from_bcf.py
    │
    ├── 02_prepare_overlap/
    │   └── prepare_overlap_matrices.py
    │      # builds overlap list + writes GE_overlap.txt + SNP_overlap.txt + overlap_summary.txt
    │
    ├── 03_covariates/
    │   ├── make_covariates.py
    │   └── preprocess_covariates.py
    │      # outputs Covariates_numeric.txt + covariates_preprocess_summary.txt
    │
    ├── 04_run_eqtl/
    │   └── eqtl_run.R
    │
    ├── 05_inspect_results/
    │




    ├── utils/
    │   ├── genotype_helpers.py
    │   ├── matrix_helpers.py
    │   └── covariate_helpers.py
    │
    └── EQTL_PIPELINE.md
```

---

# Before You Start

## Activate environment

```sh
conda activate bio_work
```

## Source config (REQUIRED)

```sh
source scripts/00_config.sh
```


This sets:
- REPO_ROOT
- RUN directories
- Output paths

---

# Pipeline Execution Order

---

## STEP 1 — Build SNP + Location Matrices

### Option A — From BCF (recommended for imputed data)

```sh
python3 scripts/eqtl/01_build_inputs/make_snp_and_snpsloc_from_bcf.py
```

### Option B — From TSV genotype table

```sh
python3 scripts/eqtl/01_build_inputs/make_snp_and_snpsloc_from_tsv.py
```

### Outputs  

```text
eqtl/   
├ SNP.txt  
└ snpsloc.txt  
```


---

## STEP 2 — Build Gene Expression Matrix

```sh
python3 scripts/eqtl/01_build_inputs/build_GE_from_expr.py
```


### Notes
- Uses SYMBOL → geneid
- Duplicate SYMBOL rows → **max expression value kept**

### Output
```text
eqtl/GE.txt
```


---

## STEP 3 — Build Gene Location Table

```sh
Rscript scripts/eqtl/01_build_inputs/make_geneloc_from_gtf.R
```

### Output
```text
eqtl/geneloc.txt
```


---

## STEP 4 — Harmonize Sample Overlap

```sh
python3 scripts/eqtl/02_prepare_overlap/prepare_overlap_matrices.py
```


### What this does
- Finds sample overlap between:
  - SNP.txt
  - GE.txt
  - metadata
- Subsets matrices to shared samples

### Outputs
```text
GE_overlap.txt    
SNP_overlap.txt  
uasg_ge_snp_meta_overlap.txt  
```


---

## STEP 5 — Build Covariate Table

```sh
python3 scripts/eqtl/03_covariates/make_covariates.py
```

### Output
```text
Covariates.txt
```


---

## STEP 6 — Encode Covariates (Numeric + One-Hot)

```sh
python3 scripts/eqtl/03_covariates/preprocess_covariates.py
```


### Converts
| Variable | Encoding |
|---|---|
Sex | Female=0, Male=1  
Yes/No | No=0, Yes=1  
Diagnosis | One-hot (control reference)  
Ancestry | One-hot (optional)  

### Outputs
```text
Covariates_numeric.txt
covariates_preprocess_summary.txt
```


---

## STEP 7 — Run Matrix eQTL

```sh
Rscript scripts/eqtl/04_run_eqtl/eqtl_run.R
```


### Runs
- ALL eQTL
- CIS eQTL
- TRANS eQTL

### Outputs
```text
results/
├ eqtl_all.tsv
├ eqtl_cis.tsv
├ eqtl_trans.tsv
└ eqtl_run_summary.txt
```


---

# Key Design Decisions

## SNP Encoding
Uses **ALT allele counts (0 / 1 / 2)**

## Gene Expression Collapsing
If multiple probes map to one gene:
- Keep **maximum expression value** per SYMBOL

## Chromosome Format
All outputs use:
```text
chr1 chr2 chr3 ...
```


---

# ALL SNP Workflow (Optional)

You can run the same pipeline pointing to:
```text
output/genotype_run1/eqtl_allSNPs/
```

The script `02_prepare_overlap/prepare_overlap_matrices.py` is fully generic and can subset any `GE.txt` and `SNP.txt` pair.

Example:

```sh
python prepare_overlap_matrices.py \
  --ge   output/genotype_run1/eqtl/GE.txt \
  --snp  output/genotype_run1/eqtl_allSNPs/SNP.txt \
  --meta metadata/clinical_data.csv \
  --outdir output/genotype_run1/eqtl_allSNPs
```

**Notes:**
- No need to copy GE.txt into eqtl_allSNPs
- You can keep a single GE.txt source


No special scripts needed.

---

# Utilities

Located in:
```text
scripts/eqtl/utils/
```


Includes helpers for:
- genotype parsing
- matrix subsetting
- covariate encoding

---

# Troubleshooting

## REPO_ROOT not set
source scripts/00_config.sh


## Module import errors
Run from repo root:
```text
python3 scripts/eqtl/...
```


---

# Summary

**Pipeline order:**

1. Build SNP + snpsloc  
2. Build GE  
3. Build geneloc  
4. Prepare overlap matrices  
5. Build covariates  
6. Encode covariates  
7. Run eQTL  

---

# Maintainer Notes

If adding new genotype source:  
→ Update only **01_build_inputs**

If adding new covariates:  
→ Update **03_covariates**

---

# References

### Matrix eQTL
`Shabalin AA. Matrix eQTL: ultra fast eQTL analysis via large matrix operations. Bioinformatics (2012)`

GitHub:
https://github.com/andreyshabalin/MatrixEQTL

Documentation:
https://github.com/andreyshabalin/MatrixEQTL/blob/master/README.md