#!/usr/bin/env bash
set -euo pipefail

# 1) Ensure the correct conda env is active
if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then
  echo "Please 'conda activate array-pipeline' first." >&2
  exit 1
fi

echo "=== Environment check ==="
echo "Env: ${CONDA_DEFAULT_ENV}   (CONDA_PREFIX=${CONDA_PREFIX:-<unset>})"

# 2) bcftools present?
if ! command -v bcftools >/dev/null 2>&1; then
  echo "[ERROR] bcftools not found on PATH"; exit 1
fi

# 3) Ensure BCFTOOLS_PLUGINS is set and usable; try common fallbacks if unset
if [ -z "${BCFTOOLS_PLUGINS:-}" ]; then
  for p in "$CONDA_PREFIX/libexec/bcftools" "$CONDA_PREFIX/lib/bcftools/plugins" "$CONDA_PREFIX/share/bcftools/plugins"; do
    if [ -d "$p" ]; then export BCFTOOLS_PLUGINS="$p"; break; fi
  done
fi

echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"
if [ -z "${BCFTOOLS_PLUGINS:-}" ] || [ ! -d "$BCFTOOLS_PLUGINS" ]; then
  echo "[ERROR] BCFTOOLS_PLUGINS not set to a valid directory."; exit 1
fi
if [ ! -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then
  echo "[ERROR] gtc2vcf.so not found in $BCFTOOLS_PLUGINS"; ls -1 "$BCFTOOLS_PLUGINS" || true; exit 1
fi

# 4) Core checks
echo "[OK] bcftools: $(bcftools --version | head -n 1)"
echo "[OK] plugin dir contains gtc2vcf.so"
echo "[OK] plugin list (first few):"
bcftools plugin -l | head || true

echo
echo "[OK] bcftools +gtc2vcf smoke test:"
bcftools +gtc2vcf -h | head -n 12 || { echo "[ERROR] plugin help failed"; exit 1; }

# 5) Optional helpers
if command -v samtools >/dev/null 2>&1; then
  echo "[OK] samtools: $(samtools --version | head -n 1)"
else
  echo "[WARN] samtools not found (needed to index FASTA)."
fi
if command -v plink2 >/dev/null 2>&1; then
  echo "[OK] plink2: $(plink2 --version | head -n 1)"
else
  echo "[WARN] plink2 not found."
fi

echo
echo "=== Environment looks good ==="

