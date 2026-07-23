#!/usr/bin/env python
"""R1.2 figure: (a) p-value histogram faceted by label (active vs inactive-biological);
(b) BH vs BY q-values; (c) gene-overlap (Jaccard) dependence diagnostic.
Shows inactive biological pathways are ~uniform/conservative -> BH assumption ok;
notes B=200 granularity (min p = 1/201)."""
from xgate_paths import ROOT  # noqa: E402
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway

OUT = ROOT + ""
df = pd.read_csv(f"{OUT}/results/supp_figS11_pancreas_ctrl_focused.csv")
sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
pal = {"active": "#2166ac", "inactive": "#b2182b"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# (a) p-value histograms faceted by label
bins = np.linspace(0, 1, 21)
for lab in ["inactive", "active"]:
    ax[0].hist(df[df.label==lab].p_value, bins=bins, alpha=.6, color=pal[lab],
               label=f"{lab} (n={ (df.label==lab).sum() })", density=True, edgecolor="white")
ax[0].axhline(1.0, ls="--", c="k", lw=1, label="uniform")
ax[0].axvline(1/201, ls=":", c="gray", lw=1)
ax[0].set_xlabel("raw xGATE p-value"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=8)
ax[0].set_title("p-value distribution by label\n(inactive ~uniform; active spike near 0)")

# (b) BH vs BY
for lab in ["active","inactive"]:
    s = df[df.label==lab]
    ax[1].scatter(s.q_BH, s.q_BY, c=pal[lab], s=55, edgecolor="k", linewidth=.4, alpha=.85, label=lab)
ax[1].plot([0,1],[0,1], ":", c="0.6"); ax[1].axhline(.05, ls="--", c="gray", lw=1); ax[1].axvline(.05, ls="--", c="gray", lw=1)
ax[1].set_xlabel("BH q-value"); ax[1].set_ylabel("BY q-value (conservative)"); ax[1].legend(fontsize=8)
ax[1].set_title("BH vs BY agreement")

# (c) gene-overlap (Jaccard) dependence among tested pathways
cats = get_categorized_pathways()
genesets = {}
for name in df.pathway:
    ids = gather_pathways_between(name, name, cats)
    try: genesets[name] = set(get_genes_in_pathway(ids)) if ids else set()
    except Exception: genesets[name] = set()
names = [n for n in df.pathway if genesets.get(n)]
J = np.zeros((len(names), len(names)))
for i,a in enumerate(names):
    for j,b in enumerate(names):
        u = genesets[a] | genesets[b]
        J[i,j] = len(genesets[a] & genesets[b]) / len(u) if u else 0
iu = np.triu_indices(len(names), k=1)
ax[2].hist(J[iu], bins=30, color="#4d4d4d", alpha=.8)
ax[2].set_xlabel("pairwise Jaccard gene overlap"); ax[2].set_ylabel("# pathway pairs")
ax[2].set_title(f"Tests are dependent via shared genes\n(median J={np.median(J[iu]):.3f}, max={J[iu].max():.2f})")

for a, lab in zip(ax, "abc"):
    a.text(-0.08, 1.06, f"({lab})", transform=a.transAxes, fontsize=15, fontweight="bold")
sns.despine(fig); plt.tight_layout()
for e in ("png","pdf"):
    fig.savefig(f"{OUT}/figures/supp_figS11_pvalue_diagnostics.{e}", dpi=300, bbox_inches="tight")
print("inactive p ~ uniform? KS vs Unif:")
from scipy.stats import kstest
ip = df[df.label=="inactive"].p_value.values
print("  inactive KS:", kstest(ip, "uniform"))
print(f"[done] -> {OUT}/figures/supp_figS11_pvalue_diagnostics.png")
