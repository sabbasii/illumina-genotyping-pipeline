# Setup

These steps create the `array-pipeline` Conda environment and ensure bcftools plugins are available.

## 1) Create the environment
```bash
# optional, faster solver
conda install -n base -c conda-forge mamba -y
mamba env create -f environment.yml
conda activate array-pipeline

