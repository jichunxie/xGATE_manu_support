#!/usr/bin/env python
"""R1.3 batch figures (per user 2026-07-03): benchmark-restricted, ground-truth state + donor.

The batch analysis is now restricted to the SAME pathways used in the method benchmark, which
carry ground-truth state (positive=active / negative=inactive). This lets us run the SAME
two-factor decomposition as Suppl. Fig. 1c but with BOTH pathway state and donor, for the xGATE
read-out AND a matched gene-expression read-out, and compare the signal(state)-to-noise(donor)
split. Because the xGATE reconstruction-error read-out is heavy-tailed, the decomposition is on a
rank transform of the read-out (nonparametric; sign-log agrees, raw disclosed in the CSV).

Reads (results):
  supp_figS3_5_readout_anova_<ds>.csv            state/donor variance shares per method x transform
  supp_figS3_5_readout_<ds>_<donor>_emb.npz      per-(donor,pathway) xGATE + GE embeddings, state, donor
  supp_figS3_5_batch_manyset_<ds>.csv            full-panel donor recoverability (Suppl. Fig. 1 a,b)

Writes (figures):
  supp_figS3_batch_metrics.pdf/png    -> Suppl. Fig. 1 (fig:batch_metrics)
  supp_figS4_batch_pancreas.pdf/png   -> Suppl. Fig. 2 (fig:batch_pancreas)
  supp_figS5_batch_ts.pdf/png         -> Suppl. Fig. 3 (fig:batch_ts)
"""
from xgate_paths import ROOT  # noqa: E402
import glob
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

RES = Path(ROOT + "/results")
FIG = Path(ROOT + "/figures")
DS = [("pancreas", "Pancreatic β-cells"), ("ts", "TS fibroblasts")]
XG, GE = "#4C78A8", "#E15759"                 # xGATE, gene expression
STATE_C = {"positive": "#59A14F", "negative": "#B07AA1"}
DONOR_CMAP = "tab10"
TRANSFORM = "rank"                             # primary decomposition


def anova(ds):
    a = pd.read_csv(RES / f"supp_figS3_5_readout_anova_{ds}.csv")
    a = a[a["transform"] == TRANSFORM]
    return (a[a.method == "xGATE"].iloc[0], a[a.method == "GeneExpression"].iloc[0])


def load_emb(ds):
    Xx, Xg, st, don, pw = [], [], [], [], []
    for f in sorted(glob.glob(str(RES / f"supp_figS3_5_readout_{ds}_*_emb.npz"))):
        d = np.load(f, allow_pickle=True)
        Xx.append(d["xgate"]); Xg.append(d["geneexpr"])
        st += list(d["state"]); don += list(d["donor"]); pw += list(d["pathway"])
    return (np.vstack(Xx), np.vstack(Xg), np.array(st), np.array(don), np.array(pw))


def manyset(ds):
    m = pd.read_csv(RES / f"supp_figS3_5_batch_manyset_{ds}.csv")
    return (m[m.method == "xGATE"].iloc[0], m[m.method == "GeneExpression"].iloc[0])


# ---------- Suppl. Fig. 1: donor recoverability from the embeddings (a,b) ----------
# The ground-truth read-out signal-vs-batch decomposition (formerly panel c) now lives per
# dataset in Suppl. Figs. 2c and 3c, alongside the PCA maps, so it is not repeated here.
def fig_metrics():
    labels = [d[1] for d in DS]
    xs = np.arange(len(DS)); w = 0.36
    acc_x, acc_g, ch_x, ch_g, chance = [], [], [], [], []
    for ds, _ in DS:
        mx, mg = manyset(ds)
        acc_x.append(mx.donor_logreg_acc); acc_g.append(mg.donor_logreg_acc)
        ch_x.append(mx.CH_donor_raw); ch_g.append(mg.CH_donor_raw); chance.append(mx.chance)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
    ax[0].bar(xs - w / 2, acc_x, w, color=XG, label="xGATE")
    ax[0].bar(xs + w / 2, acc_g, w, color=GE, label="gene expression")
    for i, c in enumerate(chance):
        ax[0].hlines(c, i - 0.5, i + 0.5, color="0.35", ls="--", lw=1.3, label="chance" if i == 0 else None)
    ax[0].set_ylim(0, 0.6); ax[0].set_ylabel("donor-classification accuracy")
    ax[0].set_title("Donor recoverability", fontweight="bold"); ax[0].legend(frameon=False, fontsize=9)
    ax[1].bar(xs - w / 2, ch_x, w, color=XG); ax[1].bar(xs + w / 2, ch_g, w, color=GE)
    ax[1].set_yscale("log"); ax[1].set_ylabel("Calinski–Harabasz index (donor)")
    ax[1].set_title("Donor separation (CH)", fontweight="bold")
    for a, lab in zip(ax, "ab"):
        a.set_xticks(xs); a.set_xticklabels(labels)
        a.spines[["top", "right"]].set_visible(False); a.grid(axis="y", alpha=0.18)
        a.text(-0.14, 1.06, lab, transform=a.transAxes, fontsize=17, fontweight="bold", va="top")
    fig.tight_layout()
    fig.savefig(FIG / "supp_figS3_batch_metrics.pdf", bbox_inches="tight")
    fig.savefig(FIG / "supp_figS3_batch_metrics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("[metrics] acc", list(zip(labels, acc_x, acc_g)))


# ---------- Suppl. Figs. 2/3: benchmark pathways, donor + state, and read-out SNR ----------
def scatter(a, X, state, don, title, chance):
    donors = sorted(set(don.tolist()))
    cmap = plt.get_cmap(DONOR_CMAP); cidx = {d: cmap(i % 10) for i, d in enumerate(donors)}
    Xs = StandardScaler().fit_transform(np.asarray(X, float))
    p = PCA(n_components=2, random_state=0).fit_transform(Xs)
    for s, mk, fill in [("positive", "o", True), ("negative", "s", False)]:
        m = state == s
        a.scatter(p[m, 0], p[m, 1], marker=mk, s=42,
                  c=[cidx[d] for d in don[m]] if fill else "none",
                  edgecolors="none" if fill else [cidx[d] for d in don[m]],
                  linewidths=0 if fill else 1.3, alpha=0.85)
    a.set_title(title, fontsize=11, fontweight="bold")
    a.set_xlabel("PC1"); a.set_ylabel("PC2"); a.spines[["top", "right"]].set_visible(False)
    return donors, cidx


def fig_scatter(ds, title, out):
    Xx, Xg, state, don, pw = load_emb(ds)
    ax_, ag_ = anova(ds)
    mx, mg = manyset(ds)
    fig = plt.figure(figsize=(15, 5.0))
    a0 = fig.add_subplot(1, 3, 1); a1 = fig.add_subplot(1, 3, 2); a2 = fig.add_subplot(1, 3, 3)
    donors, cidx = scatter(a0, Xx, state, don, "xGATE graph embedding", mx.chance)
    scatter(a1, Xg, state, don, "gene-expression embedding", mg.chance)
    # (c) read-out variance decomposition (rank): state (signal) vs donor (batch), xGATE vs GE
    xs = np.arange(2); w = 0.36
    st = [ax_.state_pct, ag_.state_pct]; dn = [ax_.donor_pct, ag_.donor_pct]
    a2.bar(xs - w / 2, st, w, color="#59A14F", label="pathway state (signal)")
    a2.bar(xs + w / 2, dn, w, color="#B07AA1", label="donor (batch)")
    for i, (s, dd, rr) in enumerate(zip(st, dn, [ax_.signal_over_donor, ag_.signal_over_donor])):
        a2.text(i - w / 2, s + 1, f"{s:.0f}%", ha="center", fontsize=9)
        a2.text(i + w / 2, dd + 1, f"{dd:.1f}%", ha="center", fontsize=9)
        a2.text(i, max(s, dd) + 5, f"{rr:.0f}:1", ha="center", fontsize=10, fontweight="bold", color="0.25")
    a2.set_xticks(xs); a2.set_xticklabels(["xGATE", "gene expr."])
    a2.set_ylabel("% variance of read-out (rank)"); a2.set_ylim(0, max(st) * 1.3)
    a2.set_title("Read-out signal vs batch", fontsize=11, fontweight="bold")
    a2.legend(frameon=False, fontsize=9, loc="upper right")
    a2.spines[["top", "right"]].set_visible(False); a2.grid(axis="y", alpha=0.18)
    for a, lab in zip([a0, a1, a2], "abc"):
        a.text(-0.10, 1.07, lab, transform=a.transAxes, fontsize=17, fontweight="bold", va="top")
    # legends: donor colors + state markers, below the two scatter panels
    dh = [mlines.Line2D([], [], color=cidx[d], marker="o", ls="", ms=7, label=d) for d in donors]
    sh = [mlines.Line2D([], [], color="0.3", marker="o", ls="", ms=7, label="active (positive)"),
          mlines.Line2D([], [], color="0.3", marker="s", ls="", ms=7, mfc="none", label="inactive (negative)")]
    leg1 = fig.legend(handles=dh, title="Donor", frameon=False, fontsize=8, title_fontsize=9,
                      loc="lower center", bbox_to_anchor=(0.33, -0.10), ncol=min(len(donors), 7))
    fig.add_artist(leg1)
    fig.legend(handles=sh, title="Pathway state", frameon=False, fontsize=8, title_fontsize=9,
               loc="lower center", bbox_to_anchor=(0.72, -0.10), ncol=2)
    fig.suptitle(f"{title}: benchmark pathways, colored by donor and marked by ground-truth state",
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / f"{out}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{out}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[{ds}] {len(don)} points = {len(set(pw))} pathways x {len(set(don))} donors; "
          f"state xGATE {ax_.state_pct:.0f}% donor {ax_.donor_pct:.1f}% ({ax_.signal_over_donor:.0f}:1)")


if __name__ == "__main__":
    fig_metrics()
    fig_scatter("pancreas", "Pancreatic β-cells (5 donors)", "supp_figS4_batch_pancreas")
    fig_scatter("ts", "Tabula Sapiens fibroblasts (7 donors)", "supp_figS5_batch_ts")
    print("done")
