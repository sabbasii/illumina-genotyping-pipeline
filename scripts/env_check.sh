#!/usr/bin/env bash
set -euo pipefail

# ensure we're in the right env
if [ "${CONDA_DEFAULT_ENV:-}" != "array-pipeline" ]; then
  echo "Please 'conda activate array-pipeline' first." >&2
  exit 1
fi

echo "Env: $CONDA_DEFAULT_ENV"
echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"
command -v bcftools >/dev/null
bcftools --version
bcftools plugin -l | head || true
