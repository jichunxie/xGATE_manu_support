#!/usr/bin/env python
"""
Assemble benchmark figures/tables from xGATE (Python) + 5 competing methods (R outputs),
on the NORMALIZED inputs (SCTransform v2 for UMI datasets; LogNormalize for SMART-seq2
FUCCI). All 6 methods consume one identical normalized matrix per dataset.

Per dataset (4: liver, pancreas, FUCCI, TS): 6-method check/cross Supp table.
Across datasets: precision-recall scatter + extended-metrics bars
(precision/recall/F1/specificity/MCC/AUCPR) -> Fig 3a.

Ground truth = the `label` column in each *_competing_percall.csv (positive/negative);
xGATE call = p_value < 0.05, xGATE score = 1 - p_value. Competing calls/scores come from
the *_active / *_score columns. Everything computed from one source (percall) so xGATE and
the 5 competitors are scored identically.
"""
from xgate_paths import ROOT  # noqa: E402
import os, sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, seaborn as sns
sys.path.insert(0, os.path.dirname(__file__)); import plot_style as ps
from sklearn.metrics import average_precision_score, matthews_corrcoef
from statsmodels.stats.multitest import multipletests
R=ROOT + ""; RES=f"{R}/results"; FIG=f"{R}/figures"; ROUT=f"{R}/results"
METHOD_ORDER=["xGATE","ORA","AUCell","scGSEA","GESECA","PAGODA"]

# (display name, xGATE result csv, competing percall csv, marker)
DATASETS=[
 ("Liver",        "liver_xgate_sct",        "liver_sct_competing_percall",    "s"),
 ("Pancreas",     "pancreas_xgate_bench",   "pancreas_sct_competing_percall", "D"),
 ("FUCCI U2OS",   "fucci_xgate_lognorm",    "fucci_lognorm_competing_percall","o"),
 ("TS Fibroblast","ts_fibroblast_xgate_sct","ts_sct_competing_percall",       "^"),
]

def truthy(s):  # competing _active col -> bool
    return pd.Series(s).astype(str).str.upper().isin(["TRUE","1","1.0"]).values

def load_ds(xgs, pc):
    pcl=pd.read_csv(f"{RES}/{pc}.csv")
    pcl.columns=[c.strip().strip('"') for c in pcl.columns]
    pcl["pathway"]=pcl["pathway"].astype(str).str.strip().str.strip('"')
    pcl["label"]=pcl["label"].astype(str).str.strip().str.strip('"')
    xg=pd.read_csv(f"{RES}/{xgs}.csv")
    xg["pathway"]=xg["pathway"].astype(str).str.strip()
    df=pcl.merge(xg[["pathway","p_value"]], on="pathway", how="left")
    miss=df.p_value.isna().sum()
    if miss: print(f"  [warn] {xgs}: {miss}/{len(df)} pathways missing xGATE p_value")
    truth=(df.label=="positive").values
    # xGATE call = BH (Benjamini-Hochberg) FDR-corrected p_value < 0.05 (per dataset).
    _p=df.p_value.values; _fin=np.isfinite(_p); _q=np.full(len(_p),1.0)
    if _fin.any(): _q[_fin]=multipletests(_p[_fin],method="fdr_bh")[1]
    calls={"xGATE":(_q<0.05)}
    sc={"xGATE":(1-df.p_value).values}
    for m,col in [("ORA","ORA_active"),("AUCell","AUCell_active"),("scGSEA","scGSEA_active"),
                  ("GESECA","GESECA_active"),("PAGODA","PAGODA_active")]:
        calls[m]=truthy(df[col]) if col in df else np.full(len(df),np.nan)
        scol=col.replace("_active","_score")
        sc[m]=pd.to_numeric(df[scol],errors="coerce").values if scol in df else None
    return df, truth, calls, sc

def metrics(y,c,s):
    y=np.asarray(y,bool); c=np.asarray(pd.Series(c).fillna(False),bool)
    tp=int((y&c).sum());fp=int((~y&c).sum());fn=int((y&~c).sum());tn=int((~y&~c).sum())
    P=tp/(tp+fp) if tp+fp else np.nan; Rc=tp/(tp+fn) if tp+fn else np.nan
    F1=2*P*Rc/(P+Rc) if (P and Rc) else 0.0; spec=tn/(tn+fp) if tn+fp else np.nan
    mcc=matthews_corrcoef(y,c) if len(set(c))>1 and len(set(y))>1 else np.nan
    out=dict(precision=P,recall=Rc,F1=F1,specificity=spec,accuracy=(tp+tn)/len(y),MCC=mcc,
             TP=tp,FP=fp,FN=fn,TN=tn)
    sa=np.asarray(s,float) if s is not None else None
    # A competing method (e.g. scGSEA/gficf) can return NA for pathways whose genes fall
    # outside its feature space. Score such pathways as least-active (worst rank) instead
    # of voiding the whole metric, so AUCPR stays defined and comparable across methods on
    # the full labeled pathway set. Only genuinely no-signal cases (all-NA scores or a
    # single-class label vector) remain undefined.
    if sa is not None and len(set(y))>1 and np.isfinite(sa).any():
        finite=np.isfinite(sa)
        if not finite.all():
            sa=sa.copy(); sa[~finite]=np.nanmin(sa[finite])-1.0
        out["AUCPR"]=average_precision_score(y,sa)
    else: out["AUCPR"]=np.nan
    return out

# ---- per-dataset: metrics + check/cross table ----
allmet=[]
for name,xgs,pc,_ in DATASETS:
    df,truth,calls,sc=load_ds(xgs,pc)
    for m in METHOD_ORDER:
        r=metrics(truth,calls[m],sc.get(m)); r.update(method=m,dataset=name); allmet.append(r)
    # ---- check/cross Supp table ----
    ps.apply_style()
    pos=df[df.label=="positive"].pathway.tolist(); negl=df[df.label=="negative"].pathway.tolist()
    callmap={m:dict(zip(df.pathway,calls[m])) for m in METHOD_ORDER}
    def correct(p,m,is_pos):
        v=callmap[m].get(p)
        if v is None or (isinstance(v,float) and np.isnan(v)): return None
        return bool(v)==is_pos
    nm=len(METHOD_ORDER); maxlen=max(len(pos),len(negl))
    fig,axes=plt.subplots(1,2,figsize=(7+1.0*nm,0.5*maxlen+3.2),
                          gridspec_kw={"width_ratios":[max(1,len(pos)),max(1,len(negl))]})
    # SAME ylim on both panels so the Active/Inactive header boxes and the method
    # headers sit on the same physical row (shorter list just leaves blank rows below).
    for ax,(ttl,plist,isp) in zip(axes,[("Active / positive",pos,True),("Inactive / negative",negl,False)]):
        ax.set_xlim(0,nm); ax.set_ylim(0,maxlen); ax.invert_yaxis(); ax.axis("off")
        # Positive/Negative box sits ABOVE the 45-deg rotated method headers (no overlap)
        ax.text(nm/2,-2.7,ttl,ha="center",va="center",fontsize=15,fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.4",fc="#cfe2f3",ec="#9fc5e8"),clip_on=False)
        for j,m in enumerate(METHOD_ORDER):
            ax.text(j+0.5,-0.1,m,ha="left",va="bottom",fontsize=12,fontweight="bold",rotation=45)
        for i,p in enumerate(plist):
            ax.text(-0.1,i+0.5,p[:38],ha="right",va="center",fontsize=12)
            for j,m in enumerate(METHOD_ORDER):
                ok=correct(p,m,isp)
                ax.text(j+0.5,i+0.5,"✓" if ok else ("✗" if ok is False else "–"),
                        ha="center",va="center",fontsize=17,
                        color=("#2e7d32" if ok else ("#c62828" if ok is False else "#999")))
    fig.suptitle(f"{name} benchmark: correct (✓) / incorrect (✗) by method",fontsize=15,fontweight="bold",y=0.99)
    fig.tight_layout(rect=[0,0,1,0.95]); ps.save(fig,f"{FIG}/supp_figS6_9_benchmark_{name.replace(' ','_')}"); plt.close(fig)
    print(f"[table] {name} -> supp_figS6_9_benchmark_{name.replace(' ','_')}.png")

M=pd.DataFrame(allmet)[["dataset","method","precision","recall","F1","specificity","accuracy","MCC","AUCPR","TP","FP","FN","TN"]]
M.to_csv(f"{ROUT}/fig3_benchmark_metrics_bh.csv",index=False)
print("\n[combined metrics — 6 methods x 4 datasets]"); print(M.round(3).to_string(index=False))

# ---- Fig 3a: PR scatter + extended-metrics bars ----
ps.apply_style()
fig=plt.figure(figsize=(15,5.4)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.3],wspace=0.3)
axa=fig.add_subplot(gs[0]); import matplotlib.lines as ml
mk={d[0]:d[3] for d in DATASETS}
for _,r in M.iterrows():
    if not (np.isfinite(r.recall) and np.isfinite(r.precision)): continue
    axa.scatter(r.recall,r.precision,s=160,color=ps.METHOD_PALETTE.get(r.method,"#444"),
                marker=mk.get(r.dataset,"o"),edgecolor="k",linewidth=.6,alpha=.9,zorder=3)
rr=np.linspace(.01,1,100)
for f in (.5,.7,.9):
    pp=f*rr/(2*rr-f); pp[(pp<0)|(pp>1)]=np.nan; axa.plot(rr,pp,":",c="gray",lw=.7,alpha=.6)
    ylab=f/(2-f)   # precision where the constant-F1 curve meets recall=1
    axa.text(1.015, min(ylab,1.0), f"F1={f:.1f}", color="0.30", fontsize=12, fontweight="bold",
             ha="left", va="center", clip_on=False)   # bold, in right margin, clear of markers
axa.set_xlim(0,1.15);axa.set_ylim(0,1.06);axa.set_xlabel("Recall");axa.set_ylabel("Precision")
axa.set_title("Precision-recall (4 benchmarks)")
for s in ("top","right"): axa.spines[s].set_visible(False)
ps.panel_label(axa,"a",x=-0.12)
axb=fig.add_subplot(gs[1]); met6=["precision","recall","F1","specificity","MCC","AUCPR"]
long=M.melt(id_vars=["dataset","method"],value_vars=met6,var_name="metric",value_name="score")
sns.barplot(data=long,x="metric",y="score",hue="method",hue_order=METHOD_ORDER,palette=ps.METHOD_PALETTE,
            ax=axb,errorbar=("ci",0))
if axb.get_legend() is not None: axb.get_legend().remove()   # drop redundant inner legend (single right legend kept)
axb.set_ylim(-0.3,1.05);axb.set_xlabel("");axb.set_ylabel("score (mean of 4 datasets)")
axb.set_title("Extended metrics by method");plt.setp(axb.get_xticklabels(),rotation=20,ha="right")
for s in ("top","right"): axb.spines[s].set_visible(False)
ps.panel_label(axb,"b",x=-0.1)
mh=[ml.Line2D([],[],color=ps.METHOD_PALETTE.get(m,"#444"),marker="o",ls="",ms=10,label=m) for m in METHOD_ORDER]
dh=[ml.Line2D([],[],color="#444",marker=mk[d[0]],ls="",ms=10,label=d[0]) for d in DATASETS]
leg1=fig.legend(handles=mh,title="Method",loc="center left",bbox_to_anchor=(0.99,0.62),frameon=False)
fig.legend(handles=dh,title="Dataset (panel a marker)",loc="center left",bbox_to_anchor=(0.99,0.25),frameon=False)
fig.add_artist(leg1)
ps.save(fig,f"{FIG}/fig3a_pr_metrics"); plt.close(fig)
print("\n[done] -> 4 suppl tables + fig3_benchmark_metrics_bh.csv + fig3a_pr_metrics.png")
