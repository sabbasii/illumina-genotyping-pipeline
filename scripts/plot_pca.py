#!/usr/bin/env python3
"""
Minimal PCA plotter (PC1 vs PC2; PC3 vs PC4) for PLINK2 eigenvec files.
Uses matplotlib only; if not installed, prints a friendly message and exits 0.
"""
import sys, os
try:
    import matplotlib.pyplot as plt
except Exception as e:
    print("[HINT] matplotlib not available; skipping plots. Install via: pip install matplotlib")
    sys.exit(0)

import argparse
import csv

def load_eigenvec(path):
    # PLINK2 eigenvec: FID IID PC1 PC2 ...
    rows=[]
    with open(path, newline='') as fh:
        rdr=csv.reader(fh, delimiter=' ')
        for r in rdr:
            r=[c for c in r if c != '']
            if len(r) >= 6:
                rows.append((r[0], r[1], [float(x) for x in r[2:]]))
    return rows

def scatter_pcs(rows, i, j, outpng):
    x=[r[2][i] for r in rows]
    y=[r[2][j] for r in rows]
    plt.figure()
    plt.scatter(x,y,s=10)
    plt.xlabel(f"PC{i+1}")
    plt.ylabel(f"PC{j+1}")
    plt.tight_layout()
    plt.savefig(outpng, dpi=150)
    plt.close()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--eigenvec", required=True, help="Path to .eigenvec")
    ap.add_argument("--outdir", required=True, help="Output directory for PNGs")
    args=ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rows=load_eigenvec(args.eigenvec)
    if not rows:
        print("[WARN] No rows parsed from eigenvec; nothing to plot.")
        sys.exit(0)
    scatter_pcs(rows, 0, 1, os.path.join(args.outdir, "pca_pc1_pc2.png"))
    scatter_pcs(rows, 2, 3, os.path.join(args.outdir, "pca_pc3_pc4.png"))
    print("[OK] PCA plots written.")

if __name__ == "__main__":
    main()
