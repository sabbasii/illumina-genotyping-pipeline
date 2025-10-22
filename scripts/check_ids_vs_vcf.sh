#!/usr/bin/env bash
set -euo pipefail
VCF="${1:?VCF path}"; SEXMAP="${2:?sexmap.txt path}"
tmpdir=$(mktemp -d)
bcftools query -l "$VCF" | sort > "$tmpdir/vcf.samples"
cut -f1 "$SEXMAP" | sort > "$tmpdir/sexmap.samples"
echo "In VCF but missing in sexmap (first 10):"
comm -23 "$tmpdir/vcf.samples" "$tmpdir/sexmap.samples" | head || true
echo "In sexmap but not in VCF (first 10):"
comm -13 "$tmpdir/vcf.samples" "$tmpdir/sexmap.samples" | head || true
rm -rf "$tmpdir"
