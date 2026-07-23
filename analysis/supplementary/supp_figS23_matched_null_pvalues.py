#!/usr/bin/env python
"""
Supp_matched_null_pvalues.pdf  —  default (size-only) vs global-degree-matched null.

Produces the two-panel figure that quantifies how xGATE p-values change between the
default size-only null and the global-degree-matched null:
  (a) per-pathway p_default vs p_matched scatter, coloured by ground-truth activity;
  (b) inactive (null) pathways called active (BH q<0.05) per dataset, default vs matched.
Result: active detection ~unchanged (pooled 45->47 of 51); inactive false positives
more than double (pooled 9->21 of 50, concentrated in liver 1->11). Motivates keeping
the size-only null as xGATE's default.

PROVENANCE / HOW TO RUN (DCC):
  - Inputs (read-only, owned by USER):
        /path/to/xGATE/results/supp_figS23_{fucci,liver,pancreas,ts}_matched_null.csv
        produced by /path/to/xGATE/analysis/shared/matched_null.py
        (columns: pathway, truth, p_default, p_matched, ...).
  - Dependency: plot_style.py from /path/to/xGATE/analysis/shared/.
  - This script is version-controlled here because /path/to/work is read-only for USER;
    the runnable copy lives at ~/xgate_scratch/ on DCC. To reproduce:
        cd ~/xgate_scratch
        PYTHONPATH=/path/to/xGATE/analysis/shared \
        LD_LIBRARY_PATH=/path/to/conda/envs/xgate/lib \
        /path/to/conda/envs/xgate/bin/python supp_figS23_matched_null_pvalues.py
  - Output PDF is renamed to manu_v4/supp_figures/Supp_matched_null_pvalues.pdf on the laptop.
"""
from xgate_paths import ROOT  # noqa: E402
import glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import plot_style as ps
ps.apply_style()
RES=ROOT + '/results'
OUT=ROOT + '/figures/r14_pvalue_compare'
DSNAME={'liver':'Liver','pancreas':'Pancreas','fucci':'FUCCI U2OS','ts':'TS Fibroblast'}
def act(t): return str(t).lower() in ('active','positive')
FLOOR=0.004
rows=[]; bars=[]
for f in sorted(glob.glob(f'{RES}/supp_figS23_*_matched_null.csv')):
    ds=f.split('supp_figS23_')[1].split('_matched')[0]; d=pd.read_csv(f); d['a']=d.truth.map(act); rows.append(d)
    ina=~d.a
    idb=int((ina&(multipletests(d.p_default,method='fdr_bh')[1]<.05)).sum())
    imb=int((ina&(multipletests(d.p_matched,method='fdr_bh')[1]<.05)).sum())
    bars.append(dict(ds=DSNAME.get(ds,ds), n_inact=int(ina.sum()), default=idb, matched=imb))
D=pd.concat(rows,ignore_index=True); B=pd.DataFrame(bars)
fig,ax=plt.subplots(1,2,figsize=(13,5.6)); fig.subplots_adjust(wspace=0.32)
a0=ax[0]
for lab,col in [('active','#2166ac'),('inactive','#b2182b')]:
    m=(D.truth.str.lower().isin(['active','positive'])) if lab=='active' else (~D.truth.str.lower().isin(['active','positive']))
    x=np.clip(D.p_default[m],FLOOR,1); y=np.clip(D.p_matched[m],FLOOR,1)
    a0.scatter(x,y,s=34,c=col,alpha=.75,edgecolor='white',linewidth=.5,label=f'{lab} (n={int(m.sum())})',zorder=3)
a0.plot([FLOOR,1],[FLOOR,1],'--',c='0.5',lw=1.1,zorder=1)
a0.axhline(.05,ls=':',c='0.4',lw=1); a0.axvline(.05,ls=':',c='0.4',lw=1)
a0.set_xscale('log'); a0.set_yscale('log'); a0.set_xlim(FLOOR,1.05); a0.set_ylim(FLOOR,1.05)
a0.set_xlabel('p-value, size-only null (default)'); a0.set_ylabel('p-value, global-degree-matched null')
a0.set_title('Per-pathway p-value shift',pad=10); a0.legend(frameon=False,loc='lower right',fontsize=11); a0.grid(alpha=.25)
a1=ax[1]; xi=np.arange(len(B)); w=.38
a1.bar(xi-w/2,B['default'],w,color='#9ecae1',edgecolor='0.3',label='size-only null (default)')
a1.bar(xi+w/2,B['matched'],w,color='#fc9272',edgecolor='0.3',label='global-degree-matched null')
for i,r in B.iterrows():
    a1.text(i-w/2,r['default']+.15,str(r['default']),ha='center',fontsize=11)
    a1.text(i+w/2,r['matched']+.15,str(r['matched']),ha='center',fontsize=11)
a1.set_xticks(xi); a1.set_xticklabels([f"{r.ds}\n(n={r.n_inact} inact.)" for _,r in B.iterrows()],fontsize=10)
a1.set_ylabel('inactive pathways called active (BH q<0.05)')
a1.set_title('False positives among inactive pathways',pad=10); a1.legend(frameon=False,loc='upper left',fontsize=11); a1.grid(axis='y',alpha=.25)
for a,l in zip(ax,'ab'):
    a.text(-0.10,1.04,l,transform=a.transAxes,fontsize=17,fontweight='bold',va='bottom')
fig.tight_layout()
for e in ('pdf','png'): fig.savefig(f'{OUT}.{e}',dpi=300,bbox_inches='tight')
print('pooled inactive FP default->matched:',int(B['default'].sum()),'->',int(B['matched'].sum()))
print('[done] ->',OUT+'.pdf')
