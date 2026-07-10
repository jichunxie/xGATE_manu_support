#!/usr/bin/env python
"""
R1.5: reframe the computational-complexity analysis.
Show runtime/memory scale with the number of candidate gene pairs (|V|^2) and cell
count, NOT realized edge density. Adds an edge-density panel and visible panel labels.
Styling matches the original paper figure (seaborn whitegrid, hue/style=dataset_label).
"""
from xgate_paths import ROOT  # noqa: E402
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

SRC = ROOT + "/computational_benchmark_summary.csv"
ADDITIONS = ROOT + "/computational_benchmark_summary_additions.csv"
OUT = ROOT + ""

raw = pd.read_csv(SRC)
if os.path.exists(ADDITIONS):
    add = pd.read_csv(ADDITIONS)
    raw = pd.concat([raw, add], ignore_index=True)
    print(f"[input] base={len(raw) - len(add)} additions={len(add)} -> total={len(raw)}")
else:
    print(f"[input] base={len(raw)} additions=0 -> total={len(raw)}")
raw["dataset_label"] = raw["dataset_label"].fillna(raw["dataset"])
V = raw["n_vertices"].fillna(raw["n_genes"]).astype(float)
E = raw["n_edges"].astype(float)
df = pd.DataFrame({
    "dataset_label": raw["dataset_label"],
    "n_cells": pd.to_numeric(raw.get("n_cells"), errors="coerce"),
    "n_vertices": V,
    "n_edges": E,
    "cand_pairs": V * (V - 1) / 2.0,
    "edge_density": E / (V * (V - 1) / 2.0),
    "wall_time_s": pd.to_numeric(raw["wall_time_s"], errors="coerce"),
    "cpu_total_s": pd.to_numeric(raw["cpu_total_s"], errors="coerce"),
    "peak_mem_gb": pd.to_numeric(raw["peak_rss_mb"], errors="coerce") / 1024.0,
}).dropna(subset=["n_vertices", "n_edges", "cpu_total_s"]).reset_index(drop=True)
df.to_csv(f"{OUT}/results/table1_runtime_table.csv", index=False)

print(f"[loaded] {len(df)} clusters / {df.dataset_label.nunique()} studies")
print(df.groupby("dataset_label").agg(
    n=("n_vertices","size"), V_min=("n_vertices","min"), V_max=("n_vertices","max"),
    dens_min=("edge_density","min"), dens_max=("edge_density","max"),
    cpu_min=("cpu_total_s","min"), cpu_max=("cpu_total_s","max"),
    mem_max=("peak_mem_gb","max")).round(4).to_string())

def sp(a, b):
    r, p = spearmanr(a, b); return f"rho={r:+.3f} (p={p:.1e})"
print("\n[CPU-time drivers, Spearman]")
print("  cpu_total_s vs |V|^2       :", sp(df.cand_pairs, df.cpu_total_s))
print("  cpu_total_s vs |V|         :", sp(df.n_vertices, df.cpu_total_s))
print("  cpu_total_s vs n_edges     :", sp(df.n_edges, df.cpu_total_s))
print("  cpu_total_s vs edge_density:", sp(df.edge_density, df.cpu_total_s))
mc = df.dropna(subset=["n_cells"])
if len(mc) > 3:
    print(f"  cpu_total_s vs n_cells     : {sp(mc.n_cells, mc.cpu_total_s)} (n={len(mc)})")
print("[memory drivers]")
print("  peak_mem vs |V|^2          :", sp(df.cand_pairs, df.peak_mem_gb))
print("  peak_mem vs edge_density   :", sp(df.edge_density, df.peak_mem_gb))

# ---- figure: match original paper aesthetic ----
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
# fix study order + palette so colors are stable/consistent with the paper
order = ["Human CRC", "Human Prostate Cancer", "Pancreas", "Senescence", "Liver"]
order = [o for o in order if o in set(df.dataset_label)] + \
        [o for o in df.dataset_label.unique() if o not in order]
palette = dict(zip(order, sns.color_palette(n_colors=len(order))))

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

def scat(ax, x, y, xl, yl):
    sns.scatterplot(data=df, x=x, y=y, hue="dataset_label", style="dataset_label",
                    hue_order=order, style_order=order, palette=palette,
                    s=90, alpha=0.85, ax=ax, legend=(ax is axes[0,0]))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    if ax is axes[0,0]:
        ax.legend(title="Study", fontsize=8, title_fontsize=9, loc="lower right")

scat(axes[0,0], "cand_pairs", "cpu_total_s", r"Candidate gene pairs  $|V|(|V|-1)/2$", "CPU time (s)")
b, a0 = np.polyfit(np.log10(df.cand_pairs), np.log10(df.cpu_total_s), 1)
xs = np.array([df.cand_pairs.min(), df.cand_pairs.max()])
axes[0,0].plot(xs, 10**a0 * xs**b, "k--", lw=1.2, zorder=0)
axes[0,0].set_title(f"CPU time vs candidate pairs (|V|$^2$)  [slope={b:.2f}]", fontweight="bold")

scat(axes[0,1], "cand_pairs", "peak_mem_gb", r"Candidate gene pairs  $|V|(|V|-1)/2$", "Peak memory (GB)")
axes[0,1].set_title("Peak memory vs candidate pairs (|V|$^2$)", fontweight="bold")

scat(axes[1,0], "n_edges", "cpu_total_s", r"Realized edges  $|E|$", "CPU time (s)")
axes[1,0].set_title("CPU time vs realized edges (weaker)", fontweight="bold")

# edge density by study (box + strip), seaborn style
sns.boxplot(data=df, x="dataset_label", y="edge_density", order=order, hue="dataset_label",
            palette=palette, legend=False, ax=axes[1,1], showfliers=False)
sns.stripplot(data=df, x="dataset_label", y="edge_density", order=order, color=".25",
              size=3.5, alpha=0.5, ax=axes[1,1])
axes[1,1].set_yscale("log")
axes[1,1].set_xlabel("Study"); axes[1,1].set_ylabel(r"Edge density  $|E| / [|V|(|V|-1)/2]$")
axes[1,1].set_title("Edge density by study", fontweight="bold")
plt.setp(axes[1,1].get_xticklabels(), rotation=30, ha="right", fontsize=9)

for ax, lab in zip(axes.ravel(), "abcd"):
    ax.text(-0.10, 1.05, f"({lab})", transform=ax.transAxes, fontsize=16, fontweight="bold")

sns.despine(fig)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/figures/table1_complexity_analysis.{ext}", dpi=300, bbox_inches="tight")
print(f"\n[done] -> {OUT}/figures/table1_complexity_analysis.png + table")
