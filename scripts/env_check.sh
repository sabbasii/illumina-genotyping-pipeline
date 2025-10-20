#!/usr/bin/env bash
set -euo pipefail

# --- Settings ---
EXPECTED_ENV="${EXPECTED_ENV:-array-pipeline}"  # override if you use a different env name

# --- Helpers ---
info()  { printf '%s\n' "[INFO] $*"; }
ok()    { printf '%s\n' "[OK]   $*"; }
warn()  { printf '%s\n' "[WARN] $*" >&2; }
err()   { printf '%s\n' "[ERROR] $*" >&2; }

# --- Optional: ensure the right conda env ---
if [ "${CONDA_DEFAULT_ENV:-}" != "$EXPECTED_ENV" ]; then
  warn "Conda env is '${CONDA_DEFAULT_ENV:-<unset>}', expected '$EXPECTED_ENV'. Continuing, but tools must be on PATH."
fi

echo "=== Environment check ==="
echo "Env: ${CONDA_DEFAULT_ENV:-<unset>}   (CONDA_PREFIX=${CONDA_PREFIX:-<unset>})"

# --- bcftools present? ---
if ! command -v bcftools >/dev/null 2>&1; then
  err "'bcftools' not found on PATH. Try: 'conda activate $EXPECTED_ENV' or install bcftools."
  exit 1
fi
ok "bcftools: $(bcftools --version | head -n 1)"
info "bcftools binary: $(command -v bcftools)"

# --- WSL note (sometimes helpful for path issues) ---
if grep -qi microsoft /proc/version 2>/dev/null; then
  ok "WSL detected"
fi

# --- Ensure BCFTOOLS_PLUGINS points to a directory that contains gtc2vcf.so ---
# If it's already valid, keep it. Otherwise, search common locations and set it.
need_plugin_check=1
if [ -n "${BCFTOOLS_PLUGINS:-}" ] && [ -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then
  ok "BCFTOOLS_PLUGINS already set and contains gtc2vcf.so → $BCFTOOLS_PLUGINS"
  need_plugin_check=0
fi

if [ $need_plugin_check -eq 1 ]; then
  candidates=()
  # Conda-style locations
  if [ -n "${CONDA_PREFIX:-}" ]; then
    candidates+=("$CONDA_PREFIX/libexec/bcftools")
    candidates+=("$CONDA_PREFIX/lib/bcftools/plugins")
    candidates+=("$CONDA_PREFIX/share/bcftools/plugins")
  fi
  # Paths relative to bcftools binary (Homebrew/system installs)
  bcftools_bin="$(command -v bcftools)"
  bcftools_root="$(cd -- "$(dirname -- "$bcftools_bin")/.." && pwd -P 2>/dev/null || true)"
  if [ -n "$bcftools_root" ]; then
    candidates+=("$bcftools_root/libexec/bcftools")
    candidates+=("$bcftools_root/lib/bcftools/plugins")
    candidates+=("$bcftools_root/share/bcftools/plugins")
  fi
  # System fallbacks
  candidates+=("/usr/local/libexec/bcftools" "/usr/local/lib/bcftools/plugins" "/usr/local/share/bcftools/plugins")
  candidates+=("/usr/libexec/bcftools" "/usr/lib/bcftools/plugins" "/usr/share/bcftools/plugins")

  found=""
  for p in "${candidates[@]}"; do
    if [ -f "$p/gtc2vcf.so" ]; then found="$p"; break; fi
  done

  if [ -n "$found" ]; then
    export BCFTOOLS_PLUGINS="$found"
    ok "Found gtc2vcf.so at: $BCFTOOLS_PLUGINS"
  else
    warn "Could not locate gtc2vcf.so via common paths."
  fi
fi

echo "BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"

# --- Smoke tests for the plugin ---
# 1) If BCFTOOLS_PLUGINS is set and has the .so, list plugins and show help
if [ -n "${BCFTOOLS_PLUGINS:-}" ] && [ -f "$BCFTOOLS_PLUGINS/gtc2vcf.so" ]; then
  ok "plugin dir contains gtc2vcf.so"
  info "Installed plugins (first few):"
  (bcftools plugin -l | head || true)
  info "bcftools +gtc2vcf -h (first lines):"
  bcftools +gtc2vcf -h | head -n 12 || { err "Plugin help failed"; exit 1; }
else
  # 2) Maybe this build can load +gtc2vcf without BCFTOOLS_PLUGINS (some packages bake it in)
  if bcftools +gtc2vcf -h >/dev/null 2>&1; then
    ok "gtc2vcf works without BCFTOOLS_PLUGINS (package knows plugin path)"
  else
    err "bcftools plugin 'gtc2vcf' not available."
    err "Tried BCFTOOLS_PLUGINS='${BCFTOOLS_PLUGINS:-<unset>}' and implicit lookup."
    err "Tips:"
    err "  • If using conda:  conda install -c bioconda bcftools"
    err "  • Locate plugin:   find \"\${CONDA_PREFIX:-/}\" -name gtc2vcf.so 2>/dev/null"
    err "  • Then export:     export BCFTOOLS_PLUGINS=/path/to/plugins"
    exit 1
  fi
fi

# --- Optional helpers ---
if command -v samtools >/dev/null 2>&1; then
  ok "samtools: $(samtools --version | head -n 1)"
else
  warn "samtools not found (useful for FASTA indexing: samtools faidx)."
fi

if command -v plink2 >/dev/null 2>&1; then
  ok "plink2: $(plink2 --version | head -n 1)"
else
  warn "plink2 not found."
fi

echo "=== Environment looks good ==="
