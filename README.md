# Illumina Genotyping Pipeline

Linux-native, POSIX-compatible pipeline intended to run on WSL2, native Linux, and macOS.  
This repository provides reproducible scripts and configuration templates for processing Illumina Infinium genotyping array data (IDAT → GTC → VCF) using **DRAGEN Array CLI** and **bcftools** plugins.

👉 For environment setup and installation instructions, see [docs/SETUP.md](docs/SETUP.md).  
👉 For detailed requirements, see [docs/requirements.md](docs/requirements.md) (to be added).

---

## Features
- **Cross-platform:** works on WSL2 (Ubuntu), Linux, and macOS (zsh/bash).
- **Reproducible:** config-driven (`scripts/00_config.sh`).
- **Automated pipeline:** supports IDAT → GTC conversion, GTC → VCF, normalization, QC, and downstream analysis.
- **Environment isolation:** Conda/mamba-based reproducible environment with bcftools plugins.
- **Scalable:** multi-threading supported in DRAGEN and bcftools steps.
---
```text
## Directory structure
illumina-genotyping-pipeline/
├── input_data/
│ ├── idat/                       # Raw IDAT files
│ ├── manifest/                   # BPM + CSV manifest files
│ ├── sample_sheet/               # Illumina sample sheet(s)
│ └── cluster/                    # EGT cluster files
├── reference/
│ ├── GRCh37/                     # Reference FASTA
│ └── GRCh38/                     # Alternative reference build
├── output/
│ ├── genotype_run1/              # Example run outputs
│ │ ├── gtc/                      # GTC output files
│ │ ├── vcf/                      # DRAGEN VCFs
│ │ ├── cnv/                      # Copy number variation calls
│ │ ├── qc/                       # QC reports
│ │ ├── logs/                     # Run logs
│ │ └── tmp/                      # Temporary files
├── scripts/
│ ├── 00_config.sh                # Central config template
│ ├── 02_idat_to_gtc_dragena.sh
│ ├── 03_gtc_to_vcf_bcftools.sh
│ └── env_check.sh                # Checks environment & plugins
└── docs/
├── SETUP.md                      # Full environment setup guide
└── requirements.md
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

2. Set up the environment (see docs/SETUP.md):

mamba env create -f environment.yml
conda activate array-pipeline
./scripts/env_check.sh

3. Configure your run:

cp scripts/00_config.example.sh scripts/00_config.sh
nano scripts/00_config.sh   # edit paths, RUN label, reference build

4. Run the pipeline:

./scripts/02_idat_to_gtc_dragena.sh
./scripts/03_gtc_to_vcf_bcftools.sh
---
## Inputs
•   **IDATs:** per-sample intensity files (required).
•   **BPM + CSV manifest:** assay definition.
•   **EGT cluster file:** cluster definitions for genotyping.
•   **Reference FASTA:** genome build (GRCh37 or GRCh38).
•   **Sample sheet:** defines sample metadata.
## Outputs 
See [Output Files](docs/OUTPUTS.md) for details about generated files.
•   **GTC files →** genotype calls from DRAGEN.
•   **VCF/BCF →** normalized variant calls.
•   **CNV calls** (if enabled).
•   **QC logs/reports →** sanity checks, concordance metrics, sample sex checks.
•   **PLINK files →** optional export for downstream GWAS/QC.
---
## Troubleshooting
•   **Performance in WSL:** always work inside /home/<user>; avoid /mnt/c/....
•   **Plugin errors:** check that BCFTOOLS_PLUGINS is set by activate_env in 00_config.sh.
•   **Input validation:** use check_inputs_exist from 00_config.sh to confirm BPM, EGT, FASTA, and sample sheet exist.
•   **Common error logs:** stored under output/<run>/logs/.
---
## License
MIT (code).
Data and third-party tools (e.g., Illumina DRAGEN, manifests, cluster files, reference genomes) are not included.