# Pipeline Scripts: Execution Flow

This document explains **how to initialize and run the pipeline scripts**.
Follow the steps in order. Do not skip the initialization step.

---

## Before You Start

This pipeline assumes:
- Required tools are installed and available in your environment
- You are working from the **repository root**
- You have access to the input data referenced by the scripts

Installation and environment setup are documented elsewhere and are **not repeated here**.

---

## Step 1: Initialize the Pipeline Environment

All pipeline runs **start by sourcing the configuration script**:

```bash
source scripts/00_config.sh

This step defines:

+ Repository root and directory layout
+ Run identifier (e.g. genotype_run1)
+ Reference genome build (GRCh37 or GRCh38)
+ Input and output paths shared by all scripts
+ Genome-build–specific settings required by downstream tools

**Important**

`scripts/00_config.sh` must be sourced, not executed.

Correct:
`source scripts/00_config.sh`

Incorrect:
`bash scripts/00_config.sh`

If this step is skipped or executed incorrectly, downstream scripts will fail.

---

## Step 2: Verify Inputs and Environment

Before running any pipeline step, validate that all required inputs, tools, and directories are in place.

Run the verification script:

```bash
bash scripts/01_verify_inputs.sh

This script performs a non-destructive sanity check of the pipeline setup. It does not modify input data and does not run any analysis.

Specifically, it:

+ Loads scripts/00_config.sh internally for its own execution
+ Confirms required tools are available
+ Verifies presence of required genotyping inputs (manifests, cluster file, reference FASTA, sample sheet)
+ Checks for IDAT files under the expected directory
+ Creates a FASTA index (.fai) if missing and if samtools is available
+ Reports expression/microarray inputs if present
+ Ensures all expected output directories exist

---