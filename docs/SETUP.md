# Setup Guide

This guide explains how to prepare and verify the environment for the **array-pipeline** project.  
It covers installation of dependencies, creation of a reproducible conda environment, plugin setup, and verification with a helper script.

---

## 1. Install Conda (and Mamba)

You need [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (Anaconda or Miniconda).  
Check if it’s already installed:

```bash
conda --version
```

If not, install Miniconda from the official website.

---

## 2. Create the Environment

## Platform-specific conda environments

This repository defines all required tools and versions using **two platform-specific conda environment files**:

- `environment.mac.yml`
- `environment.linux.yml`

This split is intentional and necessary to ensure reliable setup across operating systems.



### --- Why this matters ---
The pipeline relies on `bcftools` for variant processing. One specific step — **Illumina GTC → VCF conversion** — uses the `bcftools +gtc2vcf` plugin.

While `bcftools` itself runs on both macOS and Linux, the `gtc2vcf` plugin is distributed and supported only on **Linux** through Bioconda. While `bcftools` runs on macOS, the `gtc2vcf` plugin is not provided as a macOS package and is available only as source code ([gtc2vcf GitHub repository](https://github.com/freeseek/gtc2vcf)).




### --- What `gtc2vcf` does ---
`gtc2vcf` converts **Illumina genotyping array output** (`.gtc` files) into standard **VCF/BCF** format. This is the only Illumina-specific step in the pipeline.

After this conversion, all downstream analyses operate on standard VCF/BCF files and use tools that are fully supported on both macOS and Linux.



### --- Platform behavior ---

**macOS**
- Supports `bcftools`, `htslib`, `samtools`, and `plink2`
- Suitable for all downstream VCF/BCF processing, QC, and analysis
- Does not support `bcftools +gtc2vcf` via conda

**Linux**
- Supports the full toolchain, including `bcftools +gtc2vcf`
- Required for Illumina GTC → VCF conversion

**WSL (Windows Subsystem for Linux)**
- Treated the same as native Linux
- Fully supported for `bcftools +gtc2vcf`
- Recommended for Windows users



### --- Choosing the correct environment ---
Select the environment file that matches your platform:

```bash
```bash
conda env create -f environment.mac.yml
conda env create -f environment.linux.yml

```

+ Both environments provide the same downstream analysis tools.
+ Only the Linux/WSL environment includes support for Illumina GTC → VCF conversion.

---

## 3. Activate the Environment

Once created, activate it:

```bash
conda activate array-pipeline
```

Verify a few tools:

```bash
bcftools --version
samtools --version
plink2 --version
```

---

## 4. Ensure Plugins are Available

This pipeline uses the **bcftools gtc2vcf plugin**.

Install (if not already available):

```bash
conda install -c bioconda bcftools=1.22 bcftools-gtc2vcf-plugin=1.22
```

On activation, the environment should automatically set the `BCFTOOLS_PLUGINS` path.  
You can check:

```bash
echo $BCFTOOLS_PLUGINS
bcftools plugin -l | head
```

```bash
bcftools --version
bcftools plugin -l | grep gtc2vcf
```

If you see plugins listed (and ideally `gtc2vcf`), you are good.

---

## 5. Verify with `01_verify_inputs.sh`

A helper script `scripts/01_verify_inputs.sh` is included to quickly confirm your setup.

### What it does

- Checks that the current conda environment is `array-pipeline`
- Prints the `BCFTOOLS_PLUGINS` path
- Confirms `bcftools` is available
- Prints version info
- Lists available plugins

### Why this matters

Instead of manually checking every tool, this script gives a one-command smoke test.  
It reduces mistakes and ensures you don’t waste time debugging later.

### How to run it

From the repo root:

```bash
conda activate array-pipeline
./scripts/01_verify_inputs.sh
```

If everything is configured, it will print versions and plugin list without errors.

---

## 6. Where to Put Data Files

The repository is structured like this:

```ini
input_data/
  ├─ idat/         # IDAT files
  ├─ manifest/     # BPM/CSV manifests
  ├─ sample_sheet/ # sample sheets
  └─ cluster/      # EGT cluster files

reference/
  ├─ GRCh37/
  └─ GRCh38/

output/
  ├─ genotype_run1/
  │  ├─ gtc/
  │  ├─ vcf/
  │  ├─ cnv/
  │  ├─ qc/
  │  └─ logs/
  └─ genotype_run2/   # future runs
```

- __Input data__ (IDAT, manifests, cluster files) → `input_data/`
- **Reference genomes** (FASTA + index) → `reference/<BUILD>/`
- **Outputs** from runs → `output/`

⚠️ These are ignored by Git (`.gitignore`), so they stay local.

---

## 7. Common Problems & Fixes

**Problem A: SSH error when cloning (`Permission denied (publickey)`)**

- Ensure your SSH key is added to your GitHub account.

- Test:

```bash
ssh -T git@github.com
```

**Problem B: bcftools plugin not found**

- Install plugin:

```bash
conda install -c bioconda bcftools-gtc2vcf-plugin=1.22
```

- Confirm `BCFTOOLS_PLUGINS` is set (run `01_verify_inputs.sh`).

**Problem C: Slow file I/O on WSL**

- Always work in `/home/<user>` (Linux filesystem).
- Avoid running jobs directly on `/mnt/c/...`.

**Problem D: DRAGEN CLI missing**

- Install Illumina DRAGEN Array CLI separately (not provided here).

- Verify with:

```bash
dragena --help
```

---

## 8. Updating the Environment

If you change `environment.yml`, update your environment with:

```bash
conda env update -f environment.yml --prune
```

To capture exact versions (lockfile):

```bash
conda env export --no-builds > environment.lock.yml
```

---

## Minimal Quick Start

```bash
git clone git@github.com:sabbasii/illumina-genotyping-pipeline.git
cd illumina-genotyping-pipeline
conda env create -f environment.yml
conda activate array-pipeline
./scripts/01_verify_inputs.sh
```

If the check passes, you are ready to run the pipeline.
