#!/usr/bin/env python
"""R1.1 FUCCI benchmark figure: xGATE vs AUCell/ORA/GESECA precision/recall/F1 on the
clean 12-positive (cell-cycle/replication/repair) / 10-negative (immune+endocrine) set."""
from xgate_paths import ROOT  # noqa: E402
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
OUT=ROOT + ""
CLEAN_NEG=["JAK-STAT signaling pathway","Toll-like receptor signaling pathway",
 "Cytokine-cytokine receptor interaction","B cell receptor signaling pathway",
 "T cell receptor signaling pathway","Insulin secretion",
 "Maturity onset diabetes of the young","Bile secretion","Salivary secretion","Taste transduction"]

xg=pd.read_csv(f"{OUT}/results/fucci_xgate_v2.csv")
cm=pd.read_csv(f"{OUT}/results/fucci_competing.csv")
cm.columns=[c.strip().strip('"') for c in cm.columns]
for c in cm.columns:
    if cm[c].dtype==object: cm[c]=cm[c].astype(str).str.strip().str.strip('"')
pos=set(xg[xg.truth=="positive"].pathway)
keep=lambda p: (p in pos) or (p in CLEAN_NEG)
xg=xg[xg.pathway.apply(keep)].copy(); cm=cm[cm.pathway.apply(keep)].copy()
truth=dict(zip(xg.pathway, xg.truth))

def prf(active_set, allp):
    tp=sum(truth[p]=="positive" and p in active_set for p in allp)
    fp=sum(truth[p]=="negative" and p in active_set for p in allp)
    fn=sum(truth[p]=="positive" and p not in active_set for p in allp)
    P=tp/(tp+fp) if tp+fp else 0.0; R=tp/(tp+fn) if tp+fn else 0.0
    F=2*P*R/(P+R) if P+R else 0.0
    return P,R,F

allp=list(truth)
calls={}
calls["xGATE"]=set(xg[xg.p_value<0.05].pathway)
for col,name in [("AUCell_active_20","AUCell"),("ORA_active","ORA"),("GESECA_active","GESECA")]:
    act=cm[cm[col].astype(str).str.upper().isin(["TRUE","1","1.0"])].pathway
    calls[name]=set(act)

rows=[]
for m in ["xGATE","AUCell","ORA","GESECA"]:
    P,R,F=prf(calls[m],allp); rows.append(dict(method=m,Precision=P,Recall=R,F1=F))
res=pd.DataFrame(rows); res.to_csv(f"{OUT}/results/fucci_benchmark_PRF.csv",index=False)
print(res.to_string(index=False))

sns.set_theme(style="whitegrid",context="paper",font_scale=1.2)
m=res.melt(id_vars="method",var_name="metric",value_name="score")
order=["xGATE","GESECA","ORA","AUCell"]
fig,ax=plt.subplots(figsize=(7.5,4.8))
sns.barplot(data=m,x="metric",y="score",hue="method",hue_order=order,
            palette={"xGATE":"#d6604d","GESECA":"#4393c3","ORA":"#92c5de","AUCell":"#bababa"},ax=ax)
ax.set_ylim(0,1.05); ax.set_ylabel("score"); ax.set_xlabel("")
ax.set_title("FUCCI U2OS cell-cycle benchmark (12 positive / 10 negative pathways)",fontsize=11,fontweight="bold")
ax.legend(title="",ncol=4,loc="lower center",fontsize=9)
for c in ax.containers: ax.bar_label(c,fmt="%.2f",fontsize=7,padding=1)
sns.despine(fig); plt.tight_layout()
for e in ("png","pdf"): fig.savefig(f"{OUT}/figures/fig3_fucci_benchmark.{e}",dpi=300,bbox_inches="tight")
print(f"[done] -> {OUT}/figures/fig3_fucci_benchmark.png")
