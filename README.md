# Illumina Genotyping Pipeline

Linux-native, POSIX-compatible pipeline intended to run on WSL2, native Linux, and macOS.  
This repository provides reproducible scripts and configuration templates for processing Illumina Infinium genotyping array data (IDAT → GTC → VCF) using **DRAGEN Array CLI** and **bcftools** plugins.

👉 For environment setup and installation instructions, see [docs/SETUP.md](docs/SETUP.md).  
👉 For detailed requirements, see [docs/requirements.md](docs/requirements.md) (to be added).

---

## Features

- **Cross-platform:** works on WSL2 (Ubuntu), Linux, and macOS (zsh/bash).
- __Reproducible:__ config-driven (`scripts/00_config.sh`).
- **Automated pipeline:** supports IDAT → GTC conversion, GTC → VCF, normalization, QC, and downstream analysis.
- **Environment isolation:** Conda/mamba-based reproducible environment with bcftools plugins.
- **Scalable:** multi-threading supported in DRAGEN and bcftools steps.

---

## Directory structure

```text
illumina-genotyping-pipeline/
├─ docs/                      
│  ├─ SETUP.md
│  ├─ sexcheck.md
│  ├─ qc_filters.md
│  └─ CHANGELOG.md
├─ env/                       # environment & reproducibility
│  ├─ environment.yml
│  └─ environment.lock.yml
├─ reference/                 # genomes/manifests that are versioned/immutable
│  ├─ GRCh37/...
│  └─ manifests/...
├─ input_data/                # raw inputs (not tracked or via .gitignore)
│  ├─ idat/
│  ├─ manifest/
│  ├─ cluster/
│  └─ sample_sheet/
├─ metadata/                  # small, text metadata that *is* tracked
│  ├─ cohort.sex.psam
│  ├─ sexmap.txt
│  ├─ overrides/             # manual curation lives here
│  │  └─ sex_overrides.txt
│  └─ runs/                  # run manifests for provenance
│     └─ genotype_run1.yaml  # parameters used for the run
├─ scripts/                   # numbered, composable CLI steps
│  ├─ 00_config.sh
│  ├─ 01_verify_inputs.sh
│  ├─ 02_idat_to_gtc_dragena.sh
│  ├─ 03_gtc_to_vcf_bcftools.sh
│  ├─ 10_qc_vcf.sh
│  ├─ 11_build_psam_from_barcode.sh
│  ├─ 12_export_sexcheck_reports.sh
│  ├─ 20_qc_filters.sh        # (new)
│  ├─ plot_pca.py             # (optional helper)
│  └─ utils/                  # tiny helpers if needed
├─ output/                    # per-run sandboxes (not tracked)
│  └─ genotype_run1/
│     ├─ logs/                # all logs, by stage, timestamped
│     ├─ tmp/                 # scratch, safe to nuke
│     ├─ gtc/
│     ├─ vcf/
│     ├─ qc/
│     │  ├─ summaries/        # txt/tsv.gz snapshots (af, hardy, missing)
│     │  ├─ sexcheck/         # sexcheck specific outputs
│     │  │  ├─ chrX.{pgen,pvar,psam}
│     │  │  ├─ cohort.sexcheck.sexcheck
│     │  │  └─ reports/       # tidy: slim.tsv, problems.tsv, histograms
│     │  ├─ plink/            # working pfiles (autosomes.*; analysis.*)
│     │  │  ├─ plink_tmp/     # throwaway intermediates from 10_qc_vcf.sh
│     │  │  ├─ autosomes.{pgen,pvar,psam}
│     │  │  ├─ analysis.clean.{pgen,pvar,psam}
│     │  │  ├─ analysis.clean.prune.in/out
│     │  │  ├─ analysis.clean.pca.{eigenvec,eigenval}
│     │  │  ├─ analysis.clean.unrel.{pgen,pvar,psam}
│     │  │  └─ analysis.clean.unrel.pca.{eigenvec,eigenval}
│     │  └─ reports/          # human-readable summaries + PNGs (PCA plots)
│     └─ cnv/                 # (if used)
├─ tests/                     # tiny fixtures + CI sanity checks (optional)
└─ README.md
```

---

## Requirements (tools)

- **Conda/Mamba** (tested with mamba/conda 24.x)
- **bcftools = 1.22**, **bcftools-gtc2vcf-plugin = 1.22**, **htslib = 1.22.1**
- **samtools** (for FASTA indexing)
- **plink2** (for downstream QC)
- **Optional:** Illumina DRAGEN Array CLI (`dragena`) — or `dragena.exe` via WSL

---

## Quick start

1. Clone the repository:

```bash
git clone https://github.com/<your-username>/illumina-genotyping-pipeline.git
cd illumina-genotyping-pipeline
```

2. Set up the environment (see docs/SETUP.md):

```bash
conda env create -f environment.yml
conda activate array-pipeline
./scripts/env_check.sh
```

3. Configure your run:

```bash
cp scripts/00_config.example.sh scripts/00_config.sh
nano scripts/00_config.sh   # edit paths, RUN label, reference build
```

4. Run the pipeline:

```bash
./scripts/02_idat_to_gtc_dragena.sh
./scripts/03_gtc_to_vcf_bcftools.sh
```

---

### Inputs and outputs

- **IDATs:** per-sample intensity files (required).
- **BPM + CSV manifest:** assay definition.
- **EGT cluster file:** cluster definitions for genotyping.
- **Reference FASTA:** genome build (GRCh37 or GRCh38).
- **Sample sheet:** defines sample metadata.

## Reference Genome (Required)

This pipeline requires a **human reference genome** that matches the **genome build used by your Illumina genotyping array and downstream tools** (DRAGEN / gtc→VCF / PLINK).

<strong> Important</strong><br>
A mismatched reference genome (wrong build or chromosome naming) will cause:
<ul>
  <li>REF/ALT allele mismatches</li>
  <li>Missing contigs</li>
  <li>gtc→VCF or DRAGEN failures</li>
</ul>

### Quick start (GRCh37 / g1k_v37)

For most **Illumina GSAMD-24v3** genotyping arrays, the manifest specifies:
- **Genome build:** GRCh37 / g1k_v37

Download the recommended reference:

```bash
ref_url=https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh37_genome.zip
wget $ref_url
unzip GRCh37_genome.zip
```

👉 Before proceeding, verify integrity and indexing (see detailed steps below).

Detailed reference documentation

For:
- Genome build background (hg19 vs GRCh37 vs GRCh38)
- 1000 Genomes & Illumina conventions
- Reference integrity checks
- Chromosome naming pitfalls

➡️ See reference/README.md.

## Outputs

See [Output Files](docs/OUTPUTS.md) for details about generated files.

- **GTC files →** genotype calls from DRAGEN.
- **VCF/BCF →** normalized variant calls.
- **CNV calls** (if enabled).
- **QC logs/reports →** sanity checks, concordance metrics, sample sex checks.
- **PLINK files →** optional export for downstream GWAS/QC.

---

## Troubleshooting

- **Performance in WSL:** always work inside `/home/<user>`; avoid `/mnt/c/...`.
- __Plugin errors:__ check that `BCFTOOLS_PLUGINS` is set by `activate_env` in `00_config.sh`.
- __Input validation:__ use `check_inputs_exist` from `00_config.sh` to confirm BPM, EGT, FASTA, and sample sheet exist.
- **Common error logs:** stored under `output/<run>/logs/`.

---

## License

MIT (code).
Data and third-party tools (e.g., Illumina DRAGEN, manifests, cluster files, reference genomes) are not included.