#!/usr/bin/env python
"""
R1.3: two-factor ANOVA of the per-donor benchmark read-outs (per user 2026-07-03).

Combines the per-donor read-out CSVs (supp_figS3_5_readout_<ds>_<donor>.csv) and runs, PER METHOD
(xGATE and matched gene expression), the same decomposition as Suppl. Fig. 1c:
    z ~ C(state) + C(donor)
with state = benchmark ground truth (positive=active / negative=inactive) and donor = batch.
Reports the signal(state)-to-noise(donor) split: variance shares, the F statistics, and the
state:donor ratio (the signal-to-noise ratio). xGATE vs gene expression are compared directly.

Writes results/supp_figS3_5_readout_anova_<ds>.csv (one row per method) and prints a summary table.
"""
from xgate_paths import ROOT  # noqa: E402
import sys
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import rankdata

RES = ROOT + "/results"
DATASETS = ["pancreas", "ts"]
METHODS = {"z_xgate": "xGATE", "z_ge": "GeneExpression"}
# The read-out z is a VAE-reconstruction-error ratio: heavy-tailed and unbounded (xGATE
# positives reach ~10^4). A raw-scale OLS is dominated by a few extreme values, so we
# decompose a monotone transform of the read-out. "rank" (nonparametric, primary) and
# "signlog" = sign(z)*log1p(|z|) both preserve the state separation while taming the tail;
# "raw" is kept for transparency.
TRANSFORMS = ["rank", "signlog", "raw"]


def transform(z, kind):
    z = np.asarray(z, float)
    if kind == "raw":
        return z
    if kind == "signlog":
        return np.sign(z) * np.log1p(np.abs(z))
    if kind == "rank":
        return rankdata(z)
    raise ValueError(kind)


def anova(df, zcol, kind):
    d = df[["donor", "state", zcol]].rename(columns={zcol: "z"}).dropna()
    d = d[np.isfinite(d["z"])].copy()
    d["z"] = transform(d["z"].values, kind)
    aov = sm.stats.anova_lm(smf.ols("z ~ C(state) + C(donor)", data=d).fit(), typ=2)
    ss = aov["sum_sq"]; tot = ss.sum()
    s_ss, d_ss, r_ss = ss["C(state)"], ss["C(donor)"], ss["Residual"]
    return dict(
        transform=kind, n=len(d), n_donor=d.donor.nunique(),
        state_pct=100 * s_ss / tot, donor_pct=100 * d_ss / tot, resid_pct=100 * r_ss / tot,
        F_state=aov["F"]["C(state)"], p_state=aov["PR(>F)"]["C(state)"],
        F_donor=aov["F"]["C(donor)"], p_donor=aov["PR(>F)"]["C(donor)"],
        signal_over_donor=s_ss / d_ss if d_ss > 0 else np.inf)


def main():
    import glob
    for ds in DATASETS:
        files = sorted(glob.glob(f"{RES}/supp_figS3_5_readout_{ds}_*.csv"))
        if not files:
            print(f"[{ds}] no per-donor files yet; skip")
            continue
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        rows = []
        for kind in TRANSFORMS:
            for zc, mname in METHODS.items():
                r = anova(df, zc, kind); r["dataset"] = ds; r["method"] = mname
                rows.append(r)
        out = pd.DataFrame(rows)[["dataset", "method", "transform", "n", "n_donor",
                                  "state_pct", "donor_pct", "resid_pct", "F_state", "p_state",
                                  "F_donor", "p_donor", "signal_over_donor"]]
        out.to_csv(f"{RES}/supp_figS3_5_readout_anova_{ds}.csv", index=False)
        print(f"\n===== {ds}  ({df.donor.nunique()} donors x {df.pathway.nunique()} pathways) =====")
        pd.set_option("display.width", 240, "display.max_columns", 30)
        print(out.round(3).to_string(index=False))
        for kind in ["rank", "signlog"]:
            x = out[(out["method"] == "xGATE") & (out["transform"] == kind)].iloc[0]
            g = out[(out["method"] == "GeneExpression") & (out["transform"] == kind)].iloc[0]
            print(f"  [{kind}] SIGNAL(state): xGATE {x.state_pct:.0f}% (F={x.F_state:.0f}) "
                  f"vs GE {g.state_pct:.0f}% (F={g.F_state:.0f}); "
                  f"BATCH(donor): xGATE {x.donor_pct:.1f}% (p={x.p_donor:.2f}) vs GE {g.donor_pct:.1f}%; "
                  f"xGATE signal:donor = {x.signal_over_donor:.0f}:1")
    print("\n[done]")


if __name__ == "__main__":
    main()
