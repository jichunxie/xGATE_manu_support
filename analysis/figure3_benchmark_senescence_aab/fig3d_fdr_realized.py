#!/usr/bin/env python
"""REAL realized-FDR control on the labeled benchmark pathways (no simulation).

For each benchmark dataset we have ground-truth positive/negative labels (the same labels
used for precision/recall/MCC). Applying BH to xGATE's p-values, at each nominal FDR level
alpha the realized false-discovery proportion is FDP(alpha) = FP/max(rejections,1), where a
rejection is a benchmark pathway with BH q <= alpha and FP = a rejected NEGATIVE. Plotting
FDP (y) vs alpha (x) with the y=x diagonal shows whether BH controls FDR: on/below the
diagonal = controlled. At alpha=0.05 the point equals 1 - precision (reproduces the paper).
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path("/path/to/group/USER/xGATE/manu_v4/review")))
import build_assets as rra

OUTDIR = "/path/to/group/USER/xGATE/manu_v4/review"
COLORS = rra.DATASET_COLORS
alphas = np.linspace(0.0, 0.2, 201)

data = {}
for name, xgate_file, competing_file, _ in rra.DATASETS:
    df, truth, calls, scores = rra.load_benchmark(xgate_file, competing_file)
    q = df.q_BH.to_numpy(); truth = np.asarray(truth, bool)
    fdp = []
    for a in alphas:
        rej = q <= a
        r = int(rej.sum()); fp = int((rej & ~truth).sum())
        fdp.append(fp / r if r else 0.0)
    data[name] = dict(fdp=np.array(fdp), n=len(df), npos=int(truth.sum()),
                      fdp05=float(np.array(fdp)[np.argmin(abs(alphas - 0.05))]))

fig, ax = plt.subplots(figsize=(6.4, 5.2))
ax.plot([0, 0.2], [0, 0.2], ":", color="0.4", lw=1.3, label="nominal (y = x)")
for name in data:
    d = data[name]
    ax.plot(alphas, d["fdp"], "-", color=COLORS[name], lw=2.2,
            label=f"{name} ({d['npos']}+/{d['n'] - d['npos']}−)")
ax.axvline(0.05, ls="--", color="0.6", lw=1)
ax.set_xlabel(r"nominal BH FDR level $\alpha$")
ax.set_ylabel("realized false-discovery proportion (benchmark labels)")
ax.set_title("BH FDR control on real labeled benchmarks", fontweight="bold")
ax.set_xlim(0, 0.2); ax.set_ylim(0, 0.35)
ax.spines[["top", "right"]].set_visible(False); ax.grid(alpha=0.18)
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUTDIR}/fig3d_fdr_realized_preview.png", dpi=200, bbox_inches="tight")
print("realized FDP at alpha=0.05 (= 1 - precision, should match paper):")
for name in data:
    print(f"  {name}: FDP {data[name]['fdp05']:.3f}  (n={data[name]['n']}, pos={data[name]['npos']})")
print(f"-> {OUTDIR}/fig3d_fdr_realized_preview.png")
