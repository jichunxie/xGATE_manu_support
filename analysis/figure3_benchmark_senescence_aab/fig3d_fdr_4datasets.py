#!/usr/bin/env python
"""
R1.2 Part 1 (4-dataset extension of fig3c_jaccard_pvalue.py).

Reviewer 1.2 dependence concern: pathway tests share genes -> p-values are
correlated -> is BH still controlling FDR? We report, for ALL FOUR benchmark
datasets (Pancreas beta-cells, Liver, FUCCI U2OS, TS Fibroblast):

  1. Pairwise Jaccard gene-overlap metrics + Spearman rho(Jaccard, |Delta p|).
     If shared genes drove the signal, |Delta p| would collapse to 0 as J -> 1.
  2. BH vs BY call-stability table: how many pathways are significant under BH
     (valid under positive dependence) vs BY (valid under ARBITRARY dependence).
     Stable calls BH<->BY => the dependence does not change the conclusions.

NO new xGATE runs: reads existing per-dataset score CSVs; KEGG gene sets are
resolved header-only (no big-graph read), the same gene space the scores used.

Outputs (all NEW names; does NOT overwrite the live 2-dataset figure):
  results/fig3d_fdr_jaccard_summary.csv     per-dataset Jaccard metrics
  results/bh_by_callstability.csv     per-dataset BH vs BY significance counts
  figures/fig3d_fdr_4datasets.{png,pdf}     2x4: (top) Jaccard vs |dp|, (bottom) p/BH/BY
"""
from xgate_paths import ROOT  # noqa: E402
import sys, itertools
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, ROOT + "/analysis/shared")
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
from plot_style import apply_style, panel_label, save

OUT = ROOT + ""
RES = f"{OUT}/results"
ALPHA = 0.05

# (display name, csv stem, p-value column already present?) -- each csv has
# `pathway` + `p_value`; truth/label kept only for reference, not used in stats.
DATASETS = [
    ("Pancreas",      "supp_figS11_pancreas_ctrl_focused"),
    ("Liver",         "liver_xgate_sct"),
    ("FUCCI U2OS",    "fucci_xgate_v2"),
    ("TS Fibroblast", "ts_fibroblast_xgate_final"),
]

def load(name, stem):
    df = pd.read_csv(f"{RES}/{stem}.csv")
    df = df[["pathway", "p_value"]].dropna().drop_duplicates("pathway").reset_index(drop=True)
    # recompute BOTH corrections from the raw p-values so all datasets are on the
    # exact same footing (pancreas csv ships q_BH/q_BY but we recompute for parity).
    df["q_BH"] = multipletests(df["p_value"], method="fdr_bh")[1]
    df["q_BY"] = multipletests(df["p_value"], method="fdr_by")[1]
    df["dataset"] = name
    return df

cats = get_categorized_pathways()
def resolve(name):
    ids = gather_pathways_between(name, name, cats)
    try:
        return set(get_genes_in_pathway(ids)) if ids else set()
    except Exception:
        return set()

def jaccard_table(df):
    gs = {n: resolve(n) for n in df["pathway"]}
    names = [n for n in df["pathway"] if gs[n]]
    pmap = dict(zip(df["pathway"], df["p_value"]))
    rows = []
    for a, b in itertools.combinations(names, 2):
        u = gs[a] | gs[b]
        if not u:
            continue
        rows.append(dict(jaccard=len(gs[a] & gs[b]) / len(u), dp=abs(pmap[a] - pmap[b])))
    return pd.DataFrame(rows), len(names)

apply_style()
fig, axes = plt.subplots(2, 4, figsize=(22, 9))
jac_rows, fdr_rows = [], []
panel = iter("abcdefgh")

for col, (name, stem) in enumerate(DATASETS):
    df = load(name, stem)
    jt, n_resolved = jaccard_table(df)
    rho, pv = (spearmanr(jt["jaccard"], jt["dp"]) if len(jt) > 2 else (np.nan, np.nan))

    # --- (top) Jaccard vs |delta p| ---
    ax = axes[0, col]
    hi = jt["jaccard"] >= 0.2
    ax.scatter(jt.loc[~hi, "jaccard"], jt.loc[~hi, "dp"], s=14, alpha=0.3, color="0.5", label="pair")
    ax.scatter(jt.loc[hi, "jaccard"], jt.loc[hi, "dp"], s=40, alpha=0.85, color="#d6604d",
               edgecolor="k", linewidth=0.3, label="high overlap (J≥0.2)")
    ax.set_xlabel("pairwise Jaccard gene overlap")
    ax.set_ylabel("|Δ p-value| between the pair")
    ax.set_title(f"{name}: ρ={rho:.2f} (n={len(jt)} pairs)", pad=12)
    ax.legend(fontsize=9, loc="upper right", framealpha=0.6)
    panel_label(ax, next(panel), y=1.13)

    # --- (bottom) raw p / BH q / BY q ---
    ax2 = axes[1, col]
    s = df.sort_values("p_value").reset_index(drop=True)
    x = np.arange(len(s))
    ax2.plot(x, s["p_value"], "o-", ms=4, lw=1, color="#2166ac", label="raw p")
    ax2.plot(x, s["q_BH"], "s-", ms=4, lw=1, color="#b2182b", label="BH q")
    ax2.plot(x, s["q_BY"], "^-", ms=4, lw=1, color="#1b7837", label="BY q (worst-case)")
    ax2.axhline(ALPHA, ls="--", c="grey", lw=1, label=f"α={ALPHA}")
    ax2.set_xlabel("pathway (ranked by raw p)")
    ax2.set_ylabel("p / q")
    n_raw = int((s["p_value"] < ALPHA).sum())
    n_bh  = int((s["q_BH"] < ALPHA).sum())
    n_by  = int((s["q_BY"] < ALPHA).sum())
    ax2.set_title(f"{name}: sig raw {n_raw} → BH {n_bh} → BY {n_by}", pad=12)
    ax2.legend(fontsize=9, loc="upper left", framealpha=0.6)
    panel_label(ax2, next(panel), y=1.13)

    jac_rows.append(dict(dataset=name, n_pathways=len(df), n_resolved=n_resolved,
                         n_pairs=len(jt), median_jaccard=round(float(jt["jaccard"].median()), 4),
                         max_jaccard=round(float(jt["jaccard"].max()), 4),
                         n_hi_overlap_pairs=int(hi.sum()),
                         spearman_rho=round(float(rho), 4), spearman_p=round(float(pv), 4)))
    # calls significant under BH but lost under the worst-case BY correction:
    flip = int(((df["q_BH"] < ALPHA) & (df["q_BY"] >= ALPHA)).sum())
    fdr_rows.append(dict(dataset=name, n_tested=len(df), n_sig_raw=n_raw,
                         n_sig_BH=n_bh, n_sig_BY=n_by, n_flip_BH_not_BY=flip,
                         calls_stable=(flip == 0)))

plt.tight_layout(h_pad=3.0, w_pad=2.2)
save(fig, f"{OUT}/figures/fig3d_fdr_4datasets")

jac_df = pd.DataFrame(jac_rows); fdr_df = pd.DataFrame(fdr_rows)
jac_df.to_csv(f"{RES}/fig3d_fdr_jaccard_summary.csv", index=False)
fdr_df.to_csv(f"{RES}/bh_by_callstability.csv", index=False)

pd.set_option("display.width", 160, "display.max_columns", 20)
print("\n=== Jaccard gene-overlap metrics (per dataset) ===")
print(jac_df.to_string(index=False))
print("\n=== BH vs BY call-stability (sig at q<0.05) ===")
print(fdr_df.to_string(index=False))
print(f"\n[done] figures/fig3d_fdr_4datasets.png + 2 summary CSVs in results/")
