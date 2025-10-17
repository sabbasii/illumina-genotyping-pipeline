# Requirements

The **array-pipeline** is Linux-native and POSIX-compatible (works on Linux, WSL2, macOS).

## Core tools (installed via conda)
- bcftools = 1.22
- bcftools-gtc2vcf-plugin = 1.22
- htslib = 1.22.1
- vcftools = 0.1.17
- samtools (for FASTA index)
- plink2
- wget, pigz (helpers)

## Illumina DRAGEN Array (proprietary)
- Install Illumina **DRAGEN Array** CLI (the `dragena` executable).
- Not included in this repository. Obtain from Illumina, install, and ensure `dragena` is on PATH.
- Verify installation with: `dragena --help`

## Reference genomes
- Place FASTA files under `reference/<BUILD>/` (e.g., `reference/GRCh37/human_g1k_v37.fasta`).
- Index with `samtools faidx`. Ensure contig names match the array manifest (often without a "chr" prefix).

## Manifests and cluster files
- Place array **BPM** + **CSV** manifests under `input_data/manifest/`.
- Place array **EGT** cluster files under `input_data/cluster/`.
- These files are proprietary/large and must not be committed (see `.gitignore`).

## Conda environment (recommended)
Recreate the software environment used by this pipeline:
```bash
conda env create -f environment.yml
conda activate array-pipeline
# bcftools plugin path inside this environment:
export BCFTOOLS_PLUGINS="${CONDA_PREFIX}/libexec/bcftools"

