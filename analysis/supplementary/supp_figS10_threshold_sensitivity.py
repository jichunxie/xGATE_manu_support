#!/usr/bin/env python
"""
(1) Threshold-sensitivity figure: AUCell at 10/20/30/50% cutoffs + scGSEA fixed vs
    F1-optimized cutoff -> recall/precision/F1 per cutoff, per dataset. Shows xGATE
    needs no tuning while threshold-based competitors are cutoff-sensitive.
(2) Extra-metrics figure: specificity, MCC, AUCPR (beyond precision/recall/F1) by method,
    per dataset, from fig3_benchmark_metrics_bh.csv.
"""
from xgate_paths import ROOT  # noqa: E402
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
sys.path.insert(0, os.path.dirname(__file__)); import plot_style as ps
R=ROOT + ""; RES,FIG=f"{R}/results",f"{R}/figures"; ps.apply_style()
# normalized percall files (SCTransform UMI; LogNormalize SMART-seq2 FUCCI)
DS=[("Liver","liver_sct_competing_percall"),("Pancreas","pancreas_sct_competing_percall"),
    ("FUCCI U2OS","fucci_lognorm_competing_percall"),("TS Fibroblast","ts_sct_competing_percall")]

def prf(active,is_pos):
    a=np.asarray(pd.Series(active).fillna(False),bool); y=np.asarray(is_pos,bool)
    tp=(a&y).sum();fp=(a&~y).sum();fn=(~a&y).sum()
    P=tp/(tp+fp) if tp+fp else 0; Rc=tp/(tp+fn) if tp+fn else 0
    return P,Rc,(2*P*Rc/(P+Rc) if P+Rc else 0)

# ---- (1) threshold sensitivity ----
rows=[]
for name,pc in DS:
    d=pd.read_csv(f"{RES}/{pc}.csv"); y=(d.label=="positive").values
    for q in [10,20,30,50]:
        col=f"AUCell_{q}_active"
        if col in d: P,Rc,F=prf(d[col].astype(str).str.upper().isin(["TRUE","1","1.0"]),y); rows.append(dict(dataset=name,method=f"AUCell @{q}%",F1=F,recall=Rc,precision=P))
    for col,lab in [("scGSEA_fixed_active","scGSEA fixed"),("scGSEA_active","scGSEA F1-opt")]:
        if col in d: P,Rc,F=prf(d[col].astype(str).str.upper().isin(["TRUE","1","1.0"]),y); rows.append(dict(dataset=name,method=lab,F1=F,recall=Rc,precision=P))
sens=pd.DataFrame(rows); sens.to_csv(f"{RES}/supp_figS10_threshold_sensitivity.csv",index=False)
fig,axes=plt.subplots(2,2,figsize=(14,9),sharey=True)
for ax,(name,_) in zip(axes.flat,DS):
    s=sens[sens.dataset==name]; x=np.arange(len(s)); w=0.27
    ax.bar(x-w,s.precision,w,label="precision",color="#4393c3")
    ax.bar(x,s.recall,w,label="recall",color="#f4a582")
    ax.bar(x+w,s.F1,w,label="F1",color="#d6604d")
    ax.set_xticks(x); ax.set_xticklabels(s.method,rotation=30,ha="right"); ax.set_ylim(0,1.05)
    ax.set_title(name,fontweight="bold"); ax.set_ylabel("score")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
h,l=axes.flat[0].get_legend_handles_labels()   # one legend at right-middle: applies to every panel
fig.legend(h,l,loc="center left",bbox_to_anchor=(0.99,0.5),frameon=False,title="metric")
fig.tight_layout(rect=[0,0,0.93,1]); ps.save(fig,f"{FIG}/supp_figS10_threshold_sensitivity"); plt.close(fig)
print("[sensitivity] ->", f"{FIG}/supp_figS10_threshold_sensitivity.png")
print(sens.round(3).to_string(index=False))

# ---- (2) extra metrics (specificity, MCC, AUCPR) ----
if os.path.exists(f"{RES}/fig3_benchmark_metrics_bh.csv"):
    M=pd.read_csv(f"{RES}/fig3_benchmark_metrics_bh.csv")
    extra=["specificity","MCC","AUCPR"]; order=["xGATE","ORA","AUCell","scGSEA","GESECA","PAGODA"]
    fig,axes=plt.subplots(1,len(extra),figsize=(6*len(extra),5),sharey=False)
    for ax,met in zip(axes,extra):
        sns.barplot(data=M,x="method",y=met,hue="dataset",order=order,ax=ax)
        ax.set_title(met,fontweight="bold"); ax.set_xlabel(""); ax.set_ylabel(met)
        ax.axhline(0,c="k",lw=.6); plt.setp(ax.get_xticklabels(),rotation=35,ha="right")
        for sp in ("top","right"): ax.spines[sp].set_visible(False)
        ax.legend_.set_visible(False)
    axes[-1].legend(title="dataset",frameon=False,loc="upper left",bbox_to_anchor=(1.02,1.0))
    fig.suptitle("Benchmark metrics beyond precision/recall/F1",fontweight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); ps.save(fig,f"{FIG}/supp_figS10_extra_metrics"); plt.close(fig)
    print("[extra-metrics] ->", f"{FIG}/supp_figS10_extra_metrics.png")
