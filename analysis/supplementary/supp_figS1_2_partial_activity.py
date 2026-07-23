#!/usr/bin/env python
"""
R2.1 figure: partial-activity gene-mix sweep (v2). Run AFTER supp_figS1_2_build_genemix.py.

Panels:
  (a) MDS like Fig2c (fig2c_embedding_comparison.ipynb): for each alpha, project the
      synthetic-pathway embedding (row 0) together with its size-matched random
      null cloud (rows 1..) via MDS(n_components=3).fit_transform(data.T). Plot
      the random cloud (grey) + the pathway point colored by alpha on the SAME
      shared MDS so the pathway point moves active(far)->inactive(near) the cloud.
      We fit one MDS on the union of all per-alpha matrices for comparability.
  (b) distance-to-random vs alpha (mean +/- SD over seeds) -- the quantitative
      monotone continuum.
  (c) p-value and z-score vs alpha (mean +/- SD) -- significance attenuating with
      activity fraction.

Inputs:
  results/supp_figS1_2_genemix_sweep.csv
  results/supp_figS1_2_genemix_embeddings_alpha*.csv
Outputs:
  figures/supp_figS1_2_partial_activity.{png,pdf}
"""
from xgate_paths import ROOT  # noqa: E402
import sys, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.manifold import MDS

sys.path.insert(0, ROOT + "/analysis/shared")
from plot_style import apply_style, panel_label, save

OUT = ROOT + ""
RES = f"{OUT}/results"
# optional tag arg: e.g. "oxphos" -> reads supp_figS1_2_genemix_oxphos_* , writes _oxphos
TAG = sys.argv[1] if len(sys.argv) > 1 else ""
STEM = f"supp_figS1_2_genemix_{TAG}" if TAG else "supp_figS1_2_genemix"
FIGSTEM = f"supp_figS1_2_partial_activity_{TAG}" if TAG else "supp_figS1_2_partial_activity"
sweep = pd.read_csv(f"{RES}/{STEM}_sweep.csv")

agg = (sweep.groupby("alpha")
       .agg(p_mean=("p_value", "mean"), p_sd=("p_value", "std"),
            z_mean=("z_score", "mean"), z_sd=("z_score", "std"),
            dist_mean=("dist_to_random", "mean"), dist_sd=("dist_to_random", "std"))
       .reset_index().sort_values("alpha"))

# load per-alpha embeddings
emb_files = sorted(glob.glob(f"{RES}/{STEM}_embeddings_alpha*.csv"))
emb = {}
for f in emb_files:
    m = re.search(r"alpha(\dp\d+)\.csv", f)
    a = float(m.group(1).replace("p", "."))
    df = pd.read_csv(f)
    emb[a] = df
alphas = sorted(emb.keys(), reverse=True)   # 1.0 (active) -> 0.0 (inactive)

apply_style()
cmap = plt.get_cmap("viridis")
acolor = {a: cmap(1 - i / max(1, len(alphas) - 1)) for i, a in enumerate(alphas)}

fig = plt.figure(figsize=(16, 5))

# ---- panel (a): shared MDS of all alpha matrices ----
ax_a = fig.add_subplot(1, 3, 1, projection="3d")
if emb:
    # stack: keep track of which rows are pathway vs random and which alpha
    blocks, meta = [], []
    for a in alphas:
        df = emb[a]
        X = df.drop(columns=["row_type"]).values
        blocks.append(X)
        for rt in df["row_type"]:
            meta.append((a, rt))
    # embedding feature dim varies with subgraph size across alpha -> right-pad
    # every block with zeros to the global max width before stacking for one MDS.
    maxw = max(b.shape[1] for b in blocks)
    blocks = [np.pad(b, ((0, 0), (0, maxw - b.shape[1])), constant_values=0.0)
              if b.shape[1] < maxw else b for b in blocks]
    Xall = np.vstack(blocks)
    meta = np.array(meta, dtype=object)
    mds = MDS(n_components=3, random_state=42, normalized_stress="auto")
    Y = mds.fit_transform(Xall)
    # random clouds (grey, all alphas pooled visually)
    rmask = meta[:, 1] == "random"
    ax_a.scatter(Y[rmask, 0], Y[rmask, 1], Y[rmask, 2],
                 s=18, alpha=0.10, color="grey", marker="o")
    # pathway points, one per alpha, colored
    for a in alphas:
        pmask = (meta[:, 0] == a) & (meta[:, 1] == "pathway")
        ax_a.scatter(Y[pmask, 0], Y[pmask, 1], Y[pmask, 2],
                     s=260, marker="^", edgecolors="black", linewidths=1.5,
                     color=acolor[a], label=f"α={a:.2f}")
    ax_a.set_xticks([]); ax_a.set_yticks([]); ax_a.set_zticks([])
    ax_a.set_title("Embedding space (MDS)", pad=6, y=1.0)
    # alpha encoded by color -> use a colorbar instead of a 9-entry legend that
    # overlapped the cloud. Keep only a 2-entry shape legend (placed below the panel).
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    sm = ScalarMappable(cmap=cmap, norm=Normalize(0, 1)); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax_a, location="left", shrink=0.55, pad=0.02, aspect=18)
    cb.set_label("activity fraction α  (1=active → 0=inactive)", fontsize=10)
    shape_handles = [
        Line2D([0], [0], marker="^", color="w", markerfacecolor=cmap(0.85),
               markeredgecolor="black", markersize=13, label="pathway graph"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markersize=12, alpha=0.5, label="random graph embeddings"),
    ]
    ax_a.legend(handles=shape_handles, fontsize=9, loc="upper center",
                bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)
ax_a.text2D(-0.04, 1.12, "(a)", transform=ax_a.transAxes, fontsize=17, fontweight="bold")

# ---- panel (b): distance-to-random vs alpha ----
ax_b = fig.add_subplot(1, 3, 2)
ax_b.errorbar(agg["alpha"], agg["dist_mean"], yerr=agg["dist_sd"].fillna(0),
              marker="o", ms=8, lw=2, capsize=4, color="#d6604d")
ax_b.set_xlabel("activity fraction α  (1=active → 0=inactive)")
ax_b.set_ylabel("distance to random graph embeddings")
ax_b.set_title("Embedding distance vs activity fraction", pad=12)
ax_b.invert_xaxis()
panel_label(ax_b, "b", y=1.12)

# ---- panel (c): p-value + z-score vs alpha ----
ax_c = fig.add_subplot(1, 3, 3)
ax_c.errorbar(agg["alpha"], agg["p_mean"], yerr=agg["p_sd"].fillna(0),
              marker="s", ms=7, lw=2, capsize=4, color="#2166ac", label="p-value")
ax_c.axhline(0.05, ls="--", c="grey", lw=1)
ax_c.set_xlabel("activity fraction α")
ax_c.set_ylabel("p-value")
ax_c.set_title("Significance vs activity fraction", pad=12)
ax_c.invert_xaxis()
axc2 = ax_c.twinx()
axc2.errorbar(agg["alpha"], agg["z_mean"], yerr=agg["z_sd"].fillna(0),
              marker="^", ms=7, lw=2, capsize=4, color="#762a83", label="z-score")
axc2.set_ylabel("z-score", color="#762a83")
axc2.grid(False)
lines = [Line2D([0], [0], color="#2166ac", marker="s", label="p-value"),
         Line2D([0], [0], color="#762a83", marker="^", label="z-score")]
ax_c.legend(handles=lines, fontsize=10, loc="upper center", ncol=2)
panel_label(ax_c, "c", y=1.12)

plt.tight_layout(w_pad=2.5)
save(fig, f"{OUT}/figures/{FIGSTEM}")
print("[fig] alphas:", alphas)
print(agg.to_string(index=False))
print(f"[done] -> {OUT}/figures/{FIGSTEM}.png")
