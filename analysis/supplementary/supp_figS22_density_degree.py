#!/usr/bin/env python
"""
R1.4 (T4) figure: density/degree DISTRIBUTIONS. REVISED per PI:
  - more examples: all stored ACTIVE (insulin/AMPK/FoxO) + INACTIVE (bacterial
    invasion/Shigellosis/Salmonella) pancreas pathways (3 active + 3 inactive;
    was 1 active + 1 inactive);
  - degree panels drawn as SIDE-BY-SIDE (dodged) bars per bin, not overlaid, so
    the series are individually readable;
  - both nulls shown: the naive size-only random-node null (grey) AND the
    global-degree-matched null (orange) that matches each gene's global degree;
  - ONE panel letter per pathway (placed on the left/edge-density panel), so the
    density+degree pair for each pathway shares a single number (easy to trace).

Run AFTER supp_figS22_build_distributions.py produces supp_figS22_distributions.npz.

Per pathway (one row): (col 0) edge-density distribution -- pathway (line) vs
naive-random null vs degree-matched null (overlaid hists); (col 1) node-degree
distribution -- pathway vs naive vs degree-matched as grouped bars per degree bin.
The naive random-node null sits away from the pathway (does not reproduce its
density/degree); the degree-matched null tracks the pathway's density and degree
much more closely, yet active pathways still separate under xGATE -> the signal
is topological coordination, not a size or hubness artifact.

Output: figures/supp_figS22_density_degree.{png,pdf}  (= Supp_density_degree)
"""
from xgate_paths import ROOT  # noqa: E402
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, ROOT + "/analysis/shared")
from plot_style import apply_style, panel_label, save

OUT = ROOT + ""
RES = f"{OUT}/results"
z = np.load(f"{RES}/supp_figS22_distributions.npz", allow_pickle=True)
paths = list(z["__pathways"]); truth = list(z["__truth"])

def tag(n): return n.replace(" ", "_").replace("/", "_")
present = [(n, t) for n, t in zip(paths, truth) if f"{tag(n)}__path_dens" in z.files]
active = [(n, "active") for n, t in present if t == "active"]
inactive = [(n, "inactive") for n, t in present if t == "inactive"]
chosen = active + inactive           # all stored actives then inactives (pancreas)
nrow = len(chosen)

apply_style()
PC = "#2166ac"; NA = "#999999"; MA = "#f4a582"
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def render(rows, out):
    """Render one half (3 pathways) as taller rows with the legend at TOP-center.
    The two halves are placed as a two-page ContinuedFloat in the Supp so the long
    caption is never clipped and each subfigure has more height."""
    n = len(rows)
    fig, ax = plt.subplots(n, 2, figsize=(15, 4.0 * n),
                           gridspec_kw={"width_ratios": [1.0, 1.6]})
    if n == 1:
        ax = ax[None, :]
    for row, (name, tr, letter) in enumerate(rows):
        tg = tag(name)
        pd_ = float(z[f"{tg}__path_dens"][0])
        naive_d = z[f"{tg}__naive_dens"]; match_d = z[f"{tg}__match_dens"]
        path_deg = z[f"{tg}__path_deg"]; naive_deg = z[f"{tg}__naive_deg"]; match_deg = z[f"{tg}__match_deg"]

        # (col 0) edge-density distribution: pathway line + naive + degree-matched null hists
        a = ax[row, 0]
        lo = min(naive_d.min(), match_d.min(), pd_); hi = max(naive_d.max(), match_d.max(), pd_)
        bins = np.linspace(lo, hi, 40)
        a.hist(naive_d, bins=bins, color=NA, alpha=0.6, density=True, label="naive-random null")
        a.hist(match_d, bins=bins, color=MA, alpha=0.6, density=True, label="degree-matched null")
        a.axvline(pd_, color=PC, lw=2.5, label="pathway")
        a.set_xlabel("induced edge density"); a.set_ylabel("density")
        a.set_title(f"{name[:26]} ({tr}): edge density", pad=12)
        # ONE panel letter per pathway (row), on the left panel, tracing the density+degree pair
        panel_label(a, letter, x=-0.15, y=1.16)

        # (col 1) node-degree: 3 bars CENTERED on each integer node degree, every degree labelled
        b = ax[row, 1]
        alldeg = np.concatenate([path_deg, naive_deg, match_deg])
        Dcap = int(np.clip(np.percentile(alldeg, 99), 6, 15))
        degs = np.arange(0, Dcap + 1)
        def frac(a):
            r = np.round(a).astype(int)
            return np.array([(r == d).sum() for d in degs], dtype=float) / len(a)
        w = 0.27
        b.bar(degs - w, frac(path_deg),  width=w, color=PC, align="center", label="pathway")
        b.bar(degs,     frac(naive_deg), width=w, color=NA, align="center", label="naive-random")
        b.bar(degs + w, frac(match_deg), width=w, color=MA, align="center", label="degree-matched")
        b.set_xticks(degs); b.set_xlim(-0.6, Dcap + 0.6)
        b.set_xlabel("node degree (induced subgraph)"); b.set_ylabel("fraction of nodes")
        b.set_title(f"{name[:26]} ({tr}): node degree", pad=12)

    # shared legend at TOP-center (reserve top margin so it never overlaps the panels)
    fig.legend(handles=[Line2D([0], [0], color=PC, lw=2.5, label="pathway"),
                        Patch(facecolor=NA, alpha=0.7, label="naive-random null"),
                        Patch(facecolor=MA, alpha=0.7, label="degree-matched null")],
               ncol=3, loc="upper center", bbox_to_anchor=(0.5, 0.995), frameon=True,
               fontsize=13, edgecolor="0.8")
    plt.tight_layout(rect=[0, 0, 1, 0.945], h_pad=3.2, w_pad=2.2)
    save(fig, out)
    plt.close(fig)

rows_all = [(name, tr, chr(ord("a") + i)) for i, (name, tr) in enumerate(chosen)]
half = (len(rows_all) + 1) // 2
render(rows_all[:half], f"{OUT}/figures/supp_figS22_density_degree_1")   # actives  a-c
render(rows_all[half:], f"{OUT}/figures/supp_figS22_density_degree_2")   # inactives d-f

# console summary
for name, tr in chosen:
    tg = tag(name)
    pd_ = float(z[f"{tg}__path_dens"][0])
    print(f"[{tr}] {name}: pathDens={pd_:.4f} naiveDensMed={np.median(z[f'{tg}__naive_dens']):.4f} "
          f"matchDensMed={np.median(z[f'{tg}__match_dens']):.4f} "
          f"naiveDegMed={np.median(z[f'{tg}__naive_deg']):.1f} "
          f"matchDegMed={np.median(z[f'{tg}__match_deg']):.1f} "
          f"pathDegMed={np.median(z[f'{tg}__path_deg']):.1f}")
print(f"[done, {nrow} pathways] -> {OUT}/figures/supp_figS22_density_degree.png")
