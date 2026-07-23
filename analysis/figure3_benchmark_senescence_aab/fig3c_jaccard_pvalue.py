#!/usr/bin/env python
"""
R1.2 figure S9 (4-dataset dependence + direct FDR-control check). REVISED per PI:

  Row 1  (a-d, one panel per dataset) p-value CORRELATION binned by pairwise
         Jaccard overlap. Pathway pairs are grouped into Jaccard-overlap bins
         (0-0.05, 0.05-0.10, ...). Within each bin the paired p-values are
         collected (symmetrized over pair orientation) and their Pearson
         correlation is plotted at the bin center. If shared genes drove the
         signal, the p-value correlation would climb toward 1 in the high-overlap
         bins; instead it stays low across bins -> only mild dependence.
  Row 2  (e) BH and (f) BY DIRECT FDR-control check: a Gaussian-copula
         complete-null whose pairwise correlation is set by the gene-set Jaccard
         overlap. Null p-value vectors carrying that dependence are drawn, the
         correction (BH left, BY right) applied, and the empirical FDR reported
         against the nominal cutoff. On/below the diagonal = the correction
         controls FDR under the observed overlap dependence.

(The former Row 2, raw-p/BH/BY ranked per pathway, was dropped as uninformative.)

NO new xGATE runs: reads per-dataset score CSVs; KEGG gene sets resolved
header-only. Output -> figures/fig3c_jaccard_pvalue.{png,pdf} (= Supp_jaccard_pvalue, S9)
and results/bh_by_callstability.csv (Suppl. Table 6).
"""
from xgate_paths import ROOT  # noqa: E402
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr, pearsonr, norm
from numpy.random import default_rng
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, ROOT + "/analysis/shared")
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway

OUT = ROOT + ""; RES = f"{OUT}/results"; ALPHA = 0.05
DATASETS = [("Pancreas", "supp_figS11_pancreas_ctrl_focused"), ("Liver", "liver_xgate_sct"),
            ("FUCCI U2OS", "fucci_xgate_lognorm"), ("TS Fibroblast", "ts_fibroblast_xgate_sct")]
COLORS = {"Pancreas": "#2ca02c", "Liver": "#9467bd", "FUCCI U2OS": "#8c564b", "TS Fibroblast": "#e377c2"}
BINW = 0.05
BIN_EDGES = np.arange(0.0, 0.70 + BINW, BINW)   # Jaccard overlap bins
MIN_PAIRS = 8                                    # need enough pairs for a stable correlation

cats = get_categorized_pathways()
def resolve(name):
    ids = gather_pathways_between(name, name, cats)
    try: return set(get_genes_in_pathway(ids)) if ids else set()
    except Exception: return set()

def load(stem):
    df = pd.read_csv(f"{RES}/{stem}.csv")[["pathway", "p_value"]].dropna().drop_duplicates("pathway")
    df["q_BH"] = multipletests(df.p_value, method="fdr_bh")[1]
    df["q_BY"] = multipletests(df.p_value, method="fdr_by")[1]
    return df.reset_index(drop=True)

def jaccard_matrix(names):
    gs = {n: resolve(n) for n in names}
    names = [n for n in names if gs[n]]
    m = len(names); J = np.zeros((m, m))
    for i, a in enumerate(names):
        for j in range(i + 1, m):
            u = gs[a] | gs[names[j]]
            J[i, j] = J[j, i] = (len(gs[a] & gs[names[j]]) / len(u)) if u else 0
    return J, names

def nearest_psd_corr(J):
    C = J.copy(); np.fill_diagonal(C, 1.0)
    w, V = np.linalg.eigh(C); w = np.clip(w, 1e-6, None)
    C = (V * w) @ V.T
    d = np.sqrt(np.diag(C)); return C / np.outer(d, d)

def empirical_fdr(J, alphas, method="bh", N=3000, seed=0):
    """Complete-null FDR of BH/BY under the Jaccard-induced Gaussian copula.
    All pathways null -> any rejection is false -> FDR = P(>=1 rejection)."""
    C = nearest_psd_corr(J); L = np.linalg.cholesky(C); m = C.shape[0]
    P = norm.sf(L @ default_rng(seed).standard_normal((m, N)))     # m x N null p-values
    Ps = np.sort(P, axis=0)                                        # sort each column
    ranks = np.arange(1, m + 1)[:, None]
    cm = np.sum(1.0 / np.arange(1, m + 1)) if method == "by" else 1.0
    out = []
    for a in alphas:
        rej_any = (Ps <= a * ranks / (m * cm)).any(axis=0)        # BH/BY step-up
        out.append(rej_any.mean())
    return np.array(out)

def binned_pvalue_corr(jac, pi, pj):
    """Pearson correlation of paired p-values within each Jaccard-overlap bin
    (symmetrized over pair orientation). Returns bin centers, corr, n_pairs."""
    centers, corrs, ns = [], [], []
    for lo, hi in zip(BIN_EDGES[:-1], BIN_EDGES[1:]):
        m = (jac >= lo) & (jac < hi)
        n = int(m.sum())
        centers.append((lo + hi) / 2); ns.append(n)
        if n >= MIN_PAIRS:
            a = np.concatenate([pi[m], pj[m]]); b = np.concatenate([pj[m], pi[m]])
            if a.std() > 0 and b.std() > 0:
                corrs.append(pearsonr(a, b)[0])
            else:
                corrs.append(np.nan)
        else:
            corrs.append(np.nan)
    return np.array(centers), np.array(corrs), np.array(ns)

plt.rcParams.update({"axes.titlesize": 11, "axes.labelsize": 11, "figure.dpi": 120})
fig = plt.figure(figsize=(20, 10))
gs = GridSpec(2, 4, figure=fig, hspace=0.38, wspace=0.30)
alphas = np.linspace(0.005, 0.20, 25)
rows = []
ax_bh = fig.add_subplot(gs[1, 0:2])
ax_by = fig.add_subplot(gs[1, 2:4])

for col, (name, stem) in enumerate(DATASETS):
    df = load(stem)
    J, names = jaccard_matrix(list(df.pathway))
    pmap = dict(zip(df.pathway, df.p_value))
    iu = np.triu_indices(len(names), k=1)
    jac = J[iu]
    pi = np.array([pmap[names[i]] for i in iu[0]])
    pj = np.array([pmap[names[j]] for j in iu[1]])
    dp = np.abs(pi - pj)
    rho, _ = (spearmanr(jac, dp) if len(jac) > 2 else (np.nan, np.nan))

    # Row 1: p-value correlation binned by Jaccard overlap
    a1 = fig.add_subplot(gs[0, col])
    centers, corrs, ns = binned_pvalue_corr(jac, pi, pj)
    ok = ~np.isnan(corrs)
    a1.axhline(0, color="0.7", lw=0.8, ls="--")
    a1.plot(centers[ok], corrs[ok], "o-", color=COLORS[name], ms=4.5, lw=1.4, zorder=3)
    for k, (c, r, n) in enumerate(zip(centers[ok], corrs[ok], ns[ok])):
        dy = 20 if k % 2 == 0 else -22     # large alternating offset so labels clear the line
        a1.annotate(f"n={n}", (c, r), textcoords="offset points", xytext=(0, dy),
                    ha="center", va="center", fontsize=7, color="0.10", zorder=8, clip_on=False,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", lw=0.4, alpha=0.95))
    a1.set_ylim(-1.05, 1.05); a1.set_xlim(0, BIN_EDGES[-1])
    a1.set_xlabel("pairwise Jaccard overlap (bin)"); a1.set_ylabel("p-value correlation (Pearson)")
    a1.set_title(f"{name}")
    a1.text(-0.06, 1.08, f"({'abcd'[col]})", transform=a1.transAxes, fontsize=13, fontweight="bold")

    # copula empirical-FDR curves (accumulated onto the two Row-2 panels)
    efdr_bh = empirical_fdr(J, alphas, method="bh")
    efdr_by = empirical_fdr(J, alphas, method="by")
    ax_bh.plot(alphas, efdr_bh, "-o", ms=3, color=COLORS[name], label=name)
    ax_by.plot(alphas, efdr_by, "-o", ms=3, color=COLORS[name], label=name)

    n_bh = int((df.q_BH < ALPHA).sum()); n_by = int((df.q_BY < ALPHA).sum())
    # correlation in the highest populated bin (for the caption / text)
    hi_corr = corrs[ok][-1] if ok.any() else np.nan
    rows.append(dict(dataset=name, n_tested=len(df), n_sig_BH=n_bh, n_sig_BY=n_by,
                     spearman_rho=round(float(rho), 3),
                     max_bin_pcorr=round(float(hi_corr), 3) if not np.isnan(hi_corr) else np.nan,
                     emp_FDR_BH_0p05=round(float(efdr_bh[np.argmin(abs(alphas - .05))]), 3),
                     emp_FDR_BY_0p05=round(float(efdr_by[np.argmin(abs(alphas - .05))]), 3)))

for ax, ttl in ((ax_bh, "BH"), (ax_by, "BY")):
    ax.plot([0, .2], [0, .2], ":", c="k", lw=1, label="nominal (y=x)")
    ax.set_xlabel("nominal cutoff α"); ax.set_ylabel(f"empirical FDR ({ttl}, copula null)")
    ax.set_title(f"{ttl}: empirical FDR under the observed gene-overlap dependence")
    ax.set_ylim(0, 0.22); ax.legend(fontsize=9)

ax_bh.text(-0.03, 1.08, "(e)", transform=ax_bh.transAxes, fontsize=13, fontweight="bold")
ax_by.text(-0.03, 1.08, "(f)", transform=ax_by.transAxes, fontsize=13, fontweight="bold")
for e in ("png", "pdf"):
    fig.savefig(f"{OUT}/figures/fig3c_jaccard_pvalue.{e}", dpi=300, bbox_inches="tight")

tab = pd.DataFrame(rows); tab.to_csv(f"{RES}/bh_by_callstability.csv", index=False)
print(tab.to_string(index=False))
print(f"[done] -> {OUT}/figures/fig3c_jaccard_pvalue.png  (= Supp_jaccard_pvalue, S9)")
