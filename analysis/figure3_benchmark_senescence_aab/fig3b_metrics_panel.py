"""Main Fig.3 add-on: per-dataset benchmark metrics across all 6 methods, xGATE
highlighted. 2x2 grid = F1, MCC, Specificity, AUCPR (precision/recall already in the Fig.3a PR scatter).
Source: results/fig3_benchmark_metrics_bh.csv.

Notes baked in:
- Liver/PAGODA F1 = 0 (and MCC < 0): PAGODA made no correct active calls on liver
  (TP=0), so the bar is genuinely zero, not missing data.
- scGSEA AUCPR is undefined (its score distribution is degenerate -> average-precision
  not computable); shown as 'n/a'."""
from xgate_paths import ROOT  # noqa: E402
import pandas as pd, numpy as np, matplotlib.pyplot as plt
RES=ROOT + "/results"; OUT=ROOT + "/figures"
df=pd.read_csv(f"{RES}/fig3_benchmark_metrics_bh.csv")
methods=["xGATE","scGSEA","GESECA","ORA","AUCell","PAGODA"]          # xGATE first
datasets=["Liver","Pancreas","FUCCI U2OS","TS Fibroblast"]
colors={"xGATE":"#d1495b","scGSEA":"#3a7ca5","GESECA":"#2a9d8f","ORA":"#8d6cab",
        "AUCell":"#e9a14b","PAGODA":"#9aa0a6"}
panels=[("F1","F1 score"),("MCC","MCC"),("specificity","Specificity"),("AUCPR","AUCPR")]
fig,axes=plt.subplots(2,2,figsize=(11,7.2)); w=0.13
for ax,(metric,lab) in zip(axes.flat,panels):
    x=np.arange(len(datasets))
    for j,m in enumerate(methods):
        vals=[df[(df.dataset==d)&(df.method==m)][metric].values[0] for d in datasets]
        for k,v in enumerate(vals):                       # mark undefined (NaN) values
            if not np.isfinite(v):
                ax.text(k+(j-2.5)*w,0.02,"n/a",ha="center",va="bottom",fontsize=6.5,
                        rotation=90,color=colors[m])
        ax.bar(x+(j-2.5)*w,np.nan_to_num(vals,nan=0.0),w,label=m,color=colors[m],
               edgecolor="black" if m=="xGATE" else "none",
               linewidth=1.4 if m=="xGATE" else 0)
    ax.set_xticks(x); ax.set_xticklabels(datasets,fontsize=9)
    ax.set_ylabel(lab)
    ax.axhline(0,c="k",lw=0.7)
    lo=-0.6 if metric=="MCC" else 0.0
    ax.set_ylim(lo,1.18)                                  # headroom so annotation clears bars
    mean_m=df[df.method=="xGATE"][metric].mean()
    ax.text(0.0,1.02,f"xGATE mean {metric.replace('specificity','Specificity')} = {mean_m:.2f}",
            transform=ax.transAxes,fontsize=9.5,fontweight="bold",color=colors["xGATE"])
    for sp in ("top","right"): ax.spines[sp].set_visible(False)   # drop upper/right border
handles=[plt.Rectangle((0,0),1,1,fc=colors[m],ec="black" if m=="xGATE" else "none",
         lw=1.2 if m=="xGATE" else 0) for m in methods]
fig.legend(handles,methods,ncol=6,fontsize=10,loc="lower center",
           bbox_to_anchor=(0.5,0.995),frameon=False)
fig.tight_layout(h_pad=2.4,w_pad=2.6,rect=[0,0,1,0.97])
for ext in ("pdf","png"):
    fig.savefig(f"{OUT}/fig3b_metrics_panel.{ext}",dpi=200,bbox_inches="tight")
print("wrote fig3b_metrics_panel.pdf/.png")
print(df[df.method=="xGATE"][["dataset","F1","MCC","specificity","AUCPR"]].to_string(index=False))
