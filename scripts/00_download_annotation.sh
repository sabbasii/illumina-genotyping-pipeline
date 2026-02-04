#!/bin/bash
set -euo pipefail

# Make sure REPO_ROOT is available
: "${REPO_ROOT:?REPO_ROOT is not set. Did you source 00_config.sh?}"

ANNOT_DIR="${REPO_ROOT}/reference/annotation"
mkdir -p "$ANNOT_DIR"

wget -c \
  -P "$ANNOT_DIR" \
  ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

## How to Run
# source scripts/00_config.sh
# bash scripts/00_download_annotation.sh