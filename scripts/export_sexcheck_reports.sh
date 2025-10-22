#!/usr/bin/env bash
# Export tidy sex-check reports: problems, borderlines, histograms.
# Writes to: $OUT_DIR/qc/sexcheck_reports/
set -euo pipefail

# --- Resolve repo root and load config
_SCRIPT="${BASH_SOURCE[0]:-$0}"
REPO_ROOT="$(cd -- "$(dirname -- "$_SCRIPT")/.." && pwd -P)"
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/00_config.sh"

# --- Inputs/Outputs
FILE="${SEXCHECK_FILE:-$PLINK_DIR/cohort.sexcheck.sexcheck}"   # PLINK --check-sex output
PSAM="$PSAM_SEX"                                               # recorded sex from SampleSheet
OUT="$QC_DIR/sexcheck_reports"
TMP="$QC_DIR/tmp"
mkdir -p "$OUT" "$TMP"

# --- Guardrails
[[ -s "$FILE" ]] || { echo "[ERR] Not found: $FILE" >&2; exit 1; }
[[ -s "$PSAM" ]] || { echo "[ERR] Not found: $PSAM" >&2; exit 1; }

# --- 1) Build slim table (header-aware): IID STATUS PEDSEX SNPSEX XF YRATE
awk -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++) H[$i]=i
  print "#IID","STATUS","PEDSEX","SNPSEX","XF","YRATE"
  next
}
{
  iid    = (H["IID"]    ? $(H["IID"])    : $2)
  status = (H["STATUS"] ? $(H["STATUS"]) : $6)
  ped    = (H["PEDSEX"] ? $(H["PEDSEX"]) : (H["SEX"]?$(H["SEX"]):""))
  snp    = (H["SNPSEX"] ? $(H["SNPSEX"]) : "")
  xf     = (H["F"]      ? $(H["F"])      : (H["XF"]?$(H["XF"]):""))
  yrate  = (H["YRATE"]  ? $(H["YRATE"])  : "")
  print iid,status,ped,snp,xf,yrate
}
' "$FILE" > "$OUT/sexcheck.slim.tsv"

# --- 2) Status summary (OK/PROBLEM counts)
awk 'NR>1{c[$2]++} END{for(k in c) printf "%s\t%d\n", k, c[k]}' "$OUT/sexcheck.slim.tsv" \
  | sort > "$OUT/sexcheck.summary.txt"

# --- 3) All PROBLEM rows (with header)
awk 'NR==1 || $2=="PROBLEM"' "$OUT/sexcheck.slim.tsv" > "$OUT/sexcheck.problems.tsv"

# --- 4) Borderline list: XF in [0.10,0.30] or [0.80,0.85)
awk -F'\t' 'NR>1{
  xf=$5+0;
  if ( (xf>=0.10 && xf<=0.30) || (xf>=0.80 && xf<0.85) )
    printf "%s\tXF=%s\tSTATUS=%s\tPEDSEX=%s\tSNPSEX=%s\tYRATE=%s\n", $1,$5,$2,$3,$4,$6
}' "$OUT/sexcheck.slim.tsv" > "$OUT/sexcheck.borderline.tsv"

# --- 5) Attach recorded sex to PROBLEM rows (left join by IID)
awk -F'\t' 'NR>1{print $2"\t"$5}' "$PSAM" | sort -t$'\t' -k1,1 > "$TMP/psam.sex"
sort -t$'\t' -k1,1 "$OUT/sexcheck.problems.tsv" > "$TMP/problems.sorted.tsv" || true
join -t $'\t' -1 1 -2 1 "$TMP/problems.sorted.tsv" "$TMP/psam.sex" \
  > "$OUT/sexcheck.problems_with_recorded.tsv" || true
# Add header if file not empty
if [[ -s "$OUT/sexcheck.problems_with_recorded.tsv" ]]; then
  (echo -e "#IID\tSTATUS\tPEDSEX\tSNPSEX\tXF\tYRATE\tRECORDED_SEX"; cat "$OUT/sexcheck.problems_with_recorded.tsv") \
    > "$OUT/.tmp" && mv "$OUT/.tmp" "$OUT/sexcheck.problems_with_recorded.tsv"
fi

# --- 6) Histograms (numeric + ASCII bars)
# XF (0.1 bins)
awk -F'\t' 'NR>1 && $5!=""{b=int(($5+0)*10)/10; c[b]++}
END{for(k in c) printf "%.1f\t%d\n", k, c[k]}' "$OUT/sexcheck.slim.tsv" \
  | sort -n > "$OUT/xf_hist.tsv"

awk -F'\t' '
NR>1 && $5!=""{b=int(($5+0)*10)/10; c[b]++}
END{
  max=0; for (k in c) if (c[k]>max) max=c[k];
  for (i=0;i<=10;i++){
    b=sprintf("%.1f", i/10);
    cnt = (b in c)? c[b]:0;
    bar = int((max? 50*cnt/max:0));
    printf "%s\t%4d\t|", b, cnt;
    for (j=0;j<bar;j++) printf "#";
    print ""
  }
}' "$OUT/sexcheck.slim.tsv" > "$OUT/xf_hist_ascii.txt"

# YRATE (0.01 bins up to 0.10; capped)
awk -F'\t' 'NR>1 && $6!=""{
  y=$6+0; if (y<0) y=0; if (y>0.10) y=0.10;
  b=sprintf("%.2f", int(y*100)/100); c[b]++
}
END{for(k in c) print k"\t"c[k]}' "$OUT/sexcheck.slim.tsv" \
  | sort -n > "$OUT/yrate_hist.tsv"

awk -F'\t' '
NR>1 && $6!=""{
  y=$6+0; if (y<0) y=0; if (y>0.10) y=0.10;
  b=sprintf("%.2f", int(y*100)/100); c[b]++
}
END{
  max=0; for (k in c) if (c[k]>max) max=c[k];
  for (i=0;i<=10;i++){
    b=sprintf("%.2f", i/100);
    cnt = (b in c)? c[b]:0;
    bar = int((max? 50*cnt/max:0));
    printf "%s\t%4d\t|", b, cnt;
    for (j=0;j<bar;j++) printf "#";
    print ""
  }
}' "$OUT/sexcheck.slim.tsv" > "$OUT/yrate_hist_ascii.txt"

echo "[OK] Wrote reports to: $OUT"
ls -lh "$OUT"
