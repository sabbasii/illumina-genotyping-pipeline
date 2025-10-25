#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$HOME/git_projects/illumina-genotyping-pipeline"
cd "$REPO_ROOT"

# Load config for OUT_DIR/QC_DIR/PLINK_DIR
# shellcheck source=/dev/null
source scripts/00_config.sh

echo "Repo root: $REPO_ROOT"
echo "RUN=$RUN"
echo "OUT_DIR=$OUT_DIR"
echo "QC_DIR=$QC_DIR"
echo "PLINK_DIR=$PLINK_DIR"
echo

# Create destination folders
mkdir -p "$QC_DIR"/{summaries,sexcheck/reports,plink/plink_tmp,reports}

# helper: move-if-exists (no clobber)
mvn() { # src destdir
  local src="$1" dest="$2"
  if compgen -G "$src" > /dev/null; then
    # shellcheck disable=SC2086
    for f in $src; do
      [ -e "$f" ] || continue
      local base; base="$(basename "$f")"
      if [ -e "$dest/$base" ]; then
        echo "  [SKIP] Exists: $dest/$base"
      else
        echo "  [MOVE] $f -> $dest/"
        mv "$f" "$dest/"
      fi
    done
  else
    echo "  [MISS] $src"
  fi
}

echo "=== 1) Move chrX pfiles + sexcheck table into qc/sexcheck ==="
mvn "$PLINK_DIR/chrX.pgen"   "$QC_DIR/sexcheck"
mvn "$PLINK_DIR/chrX.pvar"   "$QC_DIR/sexcheck"
mvn "$PLINK_DIR/chrX.psam"   "$QC_DIR/sexcheck"
mvn "$PLINK_DIR/cohort.sexcheck.sexcheck" "$QC_DIR/sexcheck"

echo
echo "=== 2) Move sexcheck reports folder under qc/sexcheck/reports ==="
if [ -d "$QC_DIR/sexcheck_reports" ]; then
  # move contents, keep originals if names clash
  shopt -s nullglob
  for f in "$QC_DIR/sexcheck_reports"/*; do
    base="$(basename "$f")"
    if [ -e "$QC_DIR/sexcheck/reports/$base" ]; then
      echo "  [SKIP] Exists: $QC_DIR/sexcheck/reports/$base"
    else
      echo "  [MOVE] $f -> $QC_DIR/sexcheck/reports/"
      mv "$f" "$QC_DIR/sexcheck/reports/"
    fi
  done
  rmdir "$QC_DIR/sexcheck_reports" 2>/dev/null || true
else
  echo "  [MISS] $QC_DIR/sexcheck_reports"
fi

echo
echo "=== 3) Move general QC summaries into qc/summaries ==="
# cohort-wide summaries
mvn "$PLINK_DIR/cohort.afreq.gz"     "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.hardy.gz"     "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.smiss.gz"     "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.vmiss.gz"     "$QC_DIR/summaries"
# chrX-specific summaries (keep together in summaries)
mvn "$PLINK_DIR/cohort.chrX.afreq.gz"   "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.chrX.hardy.x.gz" "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.chrX.smiss.gz"   "$QC_DIR/summaries"
mvn "$PLINK_DIR/cohort.chrX.vmiss.gz"   "$QC_DIR/summaries"

# af.tsv.gz and other qc text from qc root → summaries
mvn "$QC_DIR/af.tsv.gz"             "$QC_DIR/summaries"
mvn "$QC_DIR/bcftools.stats.txt"    "$QC_DIR/summaries"
mvn "$QC_DIR/header.txt"            "$QC_DIR/summaries"
mvn "$QC_DIR/contigs.list"          "$QC_DIR/summaries"
mvn "$QC_DIR/primary_contigs.list"  "$QC_DIR/summaries"
mvn "$QC_DIR/ref_mismatch.summary"  "$QC_DIR/summaries"
mvn "$QC_DIR/refcheck.log"          "$QC_DIR/summaries"
mvn "$QC_DIR/samples.count"         "$QC_DIR/summaries"
mvn "$QC_DIR/samples.list"          "$QC_DIR/summaries"
mvn "$QC_DIR/variants.count"        "$QC_DIR/summaries"
mvn "$QC_DIR/variants_by_contig.txt" "$QC_DIR/summaries"

echo
echo "=== 4) Move plots into qc/reports ==="
if [ -d "$QC_DIR/vcfstats_plots" ]; then
  if [ ! -e "$QC_DIR/reports/vcfstats_plots" ]; then
    echo "  [MOVE] $QC_DIR/vcfstats_plots -> $QC_DIR/reports/"
    mv "$QC_DIR/vcfstats_plots" "$QC_DIR/reports/"
  else
    echo "  [SKIP] $QC_DIR/reports/vcfstats_plots already exists"
  fi
else
  echo "  [MISS] $QC_DIR/vcfstats_plots"
fi

echo
echo "=== 5) Ensure working pfiles remain under qc/plink ==="
# Keep autosomes.* and any analysis.* where they are; ensure plink_tmp exists
mkdir -p "$PLINK_DIR/plink_tmp"
# If autosomes.* are currently somewhere else (older layout), try to salvage:
if [ -e "$QC_DIR/plink_tmp/autosomes.pgen" ] || [ -e "$QC_DIR/plink_tmp/autosomes.pvar" ] || [ -e "$QC_DIR/plink_tmp/autosomes.psam" ]; then
  echo "  [OK] autosomes.* already in $QC_DIR/plink_tmp/"
else
  # check if autosomes.* are directly under $PLINK_DIR (rare)
  mvn "$PLINK_DIR/autosomes.pgen" "$QC_DIR/plink_tmp"
  mvn "$PLINK_DIR/autosomes.pvar" "$QC_DIR/plink_tmp"
  mvn "$PLINK_DIR/autosomes.psam" "$QC_DIR/plink_tmp"
fi

echo
echo "=== 6) Clean up stale symlinks (cohort.imiss/lmiss) if present ==="
for link in "$PLINK_DIR/cohort.imiss" "$PLINK_DIR/cohort.lmiss"; do
  if [ -L "$link" ]; then
    echo "  [RM] symlink $link"
    rm -f "$link"
  fi
done
echo "  [NOTE] Use summaries: cohort.smiss.gz and cohort.vmiss.gz in qc/summaries/"

echo
echo "=== 7) Final check ==="
echo "QC tree:"
tree -L 3 "$QC_DIR" || ls -R "$QC_DIR"

echo
echo "[DONE] Reorganization complete."
