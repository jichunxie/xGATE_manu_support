#!/usr/bin/env python
"""Proposed Fig S9 e/f: REAL FDR-control curves (no simulation).

For each real dataset, xGATE p-values are BH- and BY-adjusted. Panel e plots the number of
pathways called significant as the BH FDR level alpha is varied 0->0.2 (panel f: BY). This is
the empirical operating characteristic of the FDR-controlled procedure on the real data: at
each alpha the BH procedure guarantees expected FDR <= alpha, and the curve shows how many
discoveries that yields. The alpha=0.05 marker reproduces the reported counts.
(Realized FDR itself is not computable on real data without ground truth; this shows the
FDR-controlled decision rule applied across levels, using only the real q-values.)
"""
from xgate_paths import ROOT  # noqa: E402
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

RES = ROOT + "/results"
OUTDIR = "/path/to/group/USER/xGATE/manu_v4/review"
FILES = [("Pancreas", "supp_figS11_pancreas_ctrl_focused.csv", "#4C78A8"),
         ("Liver", "liver_xgate_sct.csv", "#F28E2B"),
         ("FUCCI U2OS", "fucci_xgate_lognorm.csv", "#59A14F"),
         ("TS Fibroblast", "ts_fibroblast_xgate_sct.csv", "#E15759")]

alphas = np.linspace(0, 0.2, 101)
data = {}
for name, fn, c in FILES:
    df = pd.read_csv(f"{RES}/{fn}").dropna(subset=["p_value"]).drop_duplicates("pathway")
    q_bh = multipletests(df.p_value, method="fdr_bh")[1]
    q_by = multipletests(df.p_value, method="fdr_by")[1]
    data[name] = dict(color=c, n=len(df),
                      bh=np.array([(q_bh <= a).sum() for a in alphas]),
                      by=np.array([(q_by <= a).sum() for a in alphas]),
                      bh05=int((q_bh <= 0.05).sum()), by05=int((q_by <= 0.05).sum()))

fig, (axe, axf) = plt.subplots(1, 2, figsize=(13, 4.6))
for ax, key, ttl, lab in [(axe, "bh", "BH: significant pathways vs FDR level", "e"),
                          (axf, "by", "BY: significant pathways vs FDR level", "f")]:
    for name in data:
        d = data[name]
        ax.plot(alphas, d[key], "-", color=d["color"], lw=2,
                label=f"{name} (n={d['n']})")
    ax.axvline(0.05, ls=":", color="0.4", lw=1.2)
    ax.text(0.052, ax.get_ylim()[1] * 0.02, r"$\alpha=0.05$", color="0.35", fontsize=9)
    ax.set_xlabel(r"FDR level $\alpha$ (q-value threshold)")
    ax.set_ylabel("pathways called significant")
    ax.set_title(ttl, fontweight="bold")
    ax.set_xlim(0, 0.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.18)
    ax.text(-0.08, 1.04, lab, transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
axe.legend(frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(f"{OUTDIR}/supp_figS14_bh_asymptotic_preview.png", dpi=200, bbox_inches="tight")
print("counts at alpha=0.05 (should match paper):")
for name in data:
    print(f"  {name}: BH {data[name]['bh05']}/{data[name]['n']}, BY {data[name]['by05']}/{data[name]['n']}")
print(f"-> {OUTDIR}/supp_figS14_bh_asymptotic_preview.png")
