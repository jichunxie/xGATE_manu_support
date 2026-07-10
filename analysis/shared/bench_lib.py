"""Benchmark metrics + precision-recall plotting for the xGATE method comparison.
Reusable across datasets (FUCCI, TS fibroblast, and optionally liver/pancreas).

Each method gives, per pathway: a binary call (active) and a continuous score
(higher = more active) used for AUCPR. Metrics: precision, recall, F1, specificity,
accuracy, balanced accuracy, MCC, AUCPR.
"""
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score, matthews_corrcoef
from statsmodels.stats.multitest import multipletests

def metrics_from_calls(truth, call, score=None):
    y = np.asarray(truth, bool); c = np.asarray(call, bool)
    tp=int((y&c).sum()); fp=int((~y&c).sum()); fn=int((y&~c).sum()); tn=int((~y&~c).sum())
    P = tp/(tp+fp) if tp+fp else np.nan
    R = tp/(tp+fn) if tp+fn else np.nan
    F1 = 2*P*R/(P+R) if (P and R) else 0.0
    spec = tn/(tn+fp) if tn+fp else np.nan
    acc = (tp+tn)/len(y) if len(y) else np.nan
    bal = np.nanmean([R, spec])
    mcc = matthews_corrcoef(y, c) if (len(set(c))>1 and len(set(y))>1) else np.nan
    out = dict(precision=P, recall=R, F1=F1, specificity=spec, accuracy=acc,
               balanced_acc=bal, MCC=mcc, TP=tp, FP=fp, FN=fn, TN=tn)
    if score is not None and len(set(y))>1:
        s = np.asarray(score, float)
        finite = np.isfinite(s)
        if finite.any():
            if not finite.all():                 # NA score (method could not score a
                s = s.copy()                     # pathway) -> treat as least-active
                s[~finite] = np.nanmin(s[finite]) - 1.0
            out["AUCPR"] = average_precision_score(y, s)
    return out

def load_benchmark(xgate_csv, competing_csv, clean_neg, aucell_col="AUCell_active_20"):
    """Return per-pathway df with truth + each method's (call, score), restricted to the
    positives (from xGATE truth) + the clean negative list."""
    xg = pd.read_csv(xgate_csv)
    xg["truth"] = xg["truth"].replace({"active":"positive","inactive":"negative"})
    pos = set(xg[xg.truth=="positive"].pathway)
    keep = lambda p: (p in pos) or (p in set(clean_neg))
    xg = xg[xg.pathway.apply(keep)].copy()
    cm = pd.read_csv(competing_csv)
    cm.columns = [c.strip().strip('"') for c in cm.columns]
    for c in cm.columns:
        if cm[c].dtype==object: cm[c]=cm[c].astype(str).str.strip().str.strip('"')
    cm = cm[cm.pathway.apply(keep)].copy()
    df = xg[["pathway","truth","p_value"]].merge(cm, on="pathway", how="left")
    df["is_pos"] = df.truth=="positive"
    def b(x): return str(x).upper() in ("TRUE","1","1.0")
    methods = {}
    # xGATE call = BH (Benjamini-Hochberg) FDR-corrected p_value < 0.05, per dataset
    # (matches assemble_benchmark.py, the canonical benchmark generator).
    _p=df.p_value.values; _fin=np.isfinite(_p); _q=np.full(len(_p),1.0)
    if _fin.any(): _q[_fin]=multipletests(_p[_fin],method="fdr_bh")[1]
    methods["xGATE"]  = (_q<0.05,                    1.0-df.p_value)
    methods["AUCell"] = (df[aucell_col].map(b),      pd.to_numeric(df.get("AUCell_frac_enriched"),errors="coerce"))
    methods["ORA"]    = (df["ORA_active"].map(b),     1.0-pd.to_numeric(df["ORA_padj"],errors="coerce"))
    methods["GESECA"] = (df["GESECA_active"].map(b),  1.0-pd.to_numeric(df["GESECA_padj"],errors="coerce"))
    return df, methods

def benchmark_metrics(df, methods, dataset):
    rows=[]
    for m,(call,score) in methods.items():
        r = metrics_from_calls(df.is_pos.values, np.asarray(call), np.asarray(score))
        r.update(method=m, dataset=dataset); rows.append(r)
    return pd.DataFrame(rows)

# ---------- plots ----------
def pr_scatter(metrics_long, ax, method_palette, dataset_markers=None):
    """x=recall, y=precision; color=method, marker=dataset (Fig-3a style)."""
    import matplotlib.lines as mlines
    datasets = list(dict.fromkeys(metrics_long.dataset))
    markers = dataset_markers or dict(zip(datasets, ["o","^","s","D","v","P","X"]))
    for _,r in metrics_long.iterrows():
        ax.scatter(r["recall"], r["precision"], s=150,
                   color=method_palette.get(r["method"],"#444"),
                   marker=markers.get(r["dataset"],"o"),
                   edgecolor="k", linewidth=0.6, alpha=0.9, zorder=3)
    # F1 isolines
    import numpy as np
    rr=np.linspace(0.01,1,100)
    for f in (0.5,0.7,0.9):
        pp=f*rr/(2*rr-f); pp[(pp<0)|(pp>1)]=np.nan
        ax.plot(rr,pp,color="gray",lw=0.7,ls=":",alpha=0.6)
        ax.annotate(f"F1={f}", xy=(0.98, f*0.98/(2*0.98-f)), fontsize=8, color="gray")
    ax.set_xlim(0,1.04); ax.set_ylim(0,1.06)
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    meth_handles=[mlines.Line2D([],[],color=method_palette.get(m,"#444"),marker="o",ls="",ms=9,label=m)
                  for m in dict.fromkeys(metrics_long.method)]
    ds_handles=[mlines.Line2D([],[],color="0.3",marker=markers.get(d,"o"),ls="",ms=9,label=d)
                for d in datasets]
    return meth_handles, ds_handles
