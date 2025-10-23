#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
_SCRIPT_DIR="$(cd -- "$(dirname -- "$_SCRIPT")" && pwd -P)"
REPO_ROOT="$(cd -- "$_SCRIPT_DIR/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

ensure_dirs

LOG="$LOG_DIR/02_idat_to_gtc_dragena.log"
{
  echo "== DRAGEN IDAT->GTC =="
  echo "RUN=$RUN  REF_BUILD=$REF_BUILD"
} | tee "$LOG"

# Decide how to call DRAGEN:
# 1) native Linux 'dragena' if present
# 2) Windows 'dragena.exe' via WSL if present
# 3) allow explicit override

DRAGENA_BIN=""
MODE=""

# (A) Respect explicit override first
if [[ -n "${DRAGENA_BIN_OVERRIDE:-}" && -x "${DRAGENA_BIN_OVERRIDE}" ]]; then
  DRAGENA_BIN="$DRAGENA_BIN_OVERRIDE"
  # Guess mode from filename
  if [[ "$DRAGENA_BIN" == *.exe ]]; then MODE="windows"; else MODE="linux"; fi
fi

# (B) Auto-detect if not set by override
if [[ -z "$DRAGENA_BIN" ]]; then
  if command -v dragena >/dev/null 2>&1; then
    DRAGENA_BIN="dragena"
    MODE="linux"
  else
    for p in \
      "/mnt/c/Program Files/Illumina/DRAGEN Array/dragena.exe" \
      "/mnt/c/Program Files (x86)/Illumina/DRAGEN Array/dragena.exe"
    do
      if [[ -x "$p" ]]; then DRAGENA_BIN="$p"; MODE="windows"; break; fi
    done
  fi
fi

if [[ -z "$DRAGENA_BIN" ]]; then
  {
    echo "[ERROR] Could not find 'dragena' on PATH or 'dragena.exe' in Program Files."
    echo "        Install Illumina DRAGEN Array, or set DRAGENA_BIN_OVERRIDE to the binary."
  } | tee -a "$LOG"
  exit 1
fi

echo "[INFO] Using DRAGEN: $DRAGENA_BIN  (mode=$MODE)" | tee -a "$LOG"

# Normalize sample sheet line endings (safe no-op if already LF)
tmp_ss="$(mktemp)"
tr -d '\r' < "$SAMPLE_SHEET" > "$tmp_ss"
mv "$tmp_ss" "$SAMPLE_SHEET"

# Build command arguments; convert to Windows paths if invoking dragena.exe
if [[ "$MODE" == "windows" ]]; then
  if ! command -v wslpath >/dev/null 2>&1; then
    echo "[ERROR] 'wslpath' not found but Windows DRAGEN binary detected. Install WSL utilities." | tee -a "$LOG"
    exit 1
  fi

  IDAT_WIN=$(wslpath -w "$IDAT_DIR")
  BPM_WIN=$(wslpath -w "$BPM_MANIFEST")
  EGT_WIN=$(wslpath -w "$EGT_CLUSTER")
  SHEET_WIN=$(wslpath -w "$SAMPLE_SHEET")
  OUT_WIN=$(wslpath -w "$GTC_DIR")

  set -x
  "$DRAGENA_BIN" genotype call \
    --idat-folder "$IDAT_WIN" \
    --bpm-manifest "$BPM_WIN" \
    --cluster-file "$EGT_WIN" \
    --sample-sheet "$SHEET_WIN" \
    --num-threads "$THREADS" \
    --output-folder "$OUT_WIN" \
    2>&1 | tee -a "$LOG"
  set +x
else
  # Native Linux call
  set -x
  "$DRAGENA_BIN" genotype call \
    --idat-folder "$IDAT_DIR" \
    --bpm-manifest "$BPM_MANIFEST" \
    --cluster-file "$EGT_CLUSTER" \
    --sample-sheet "$SAMPLE_SHEET" \
    --num-threads "$THREADS" \
    --output-folder "$GTC_DIR" \
    2>&1 | tee -a "$LOG"
  set +x
fi

echo "[INFO] Checking for GTC outputs in $GTC_DIR" | tee -a "$LOG"
gtc_n=$(find "$GTC_DIR" -type f -name '*.gtc' | wc -l | tr -d ' ')
echo "GTC files: $gtc_n" | tee -a "$LOG"
if [[ "$gtc_n" -eq 0 ]]; then
  echo "[ERROR] No .gtc files produced. Check the log: $LOG" | tee -a "$LOG"
  exit 1
fi

echo "== Done: DRAGEN IDAT->GTC ==" | tee -a "$LOG"
