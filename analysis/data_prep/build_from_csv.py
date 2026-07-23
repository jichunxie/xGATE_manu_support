#!/usr/bin/env python
"""Build an xGATE co-expression graph from a genes(ensembl) x cells raw-count CSV
(e.g., liver hepatocytes), via the published pipeline; save hsa:entrez adjacency."""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, os, numpy as np, pandas as pd
from xGATE.utilities import (create_sifinet_object, quantile_thres2, cal_coexp, create_network,
                       filter_lowexp, convert_gene_ids)
ap=argparse.ArgumentParser(); ap.add_argument("--counts",required=True); ap.add_argument("--out",required=True)
a=ap.parse_args(); t0=time.time()
df=pd.read_csv(a.counts,index_col=0)            # genes x cells, ensembl rownames
df.index=[str(g).split(".")[0] for g in df.index]; df=df[~df.index.duplicated(keep="first")]
df=df.loc[(df>0).sum(axis=1)>=0.05*df.shape[1],:]
df=df.loc[:,(df>0).sum(axis=0)>=0.05*df.shape[0]]
rm=df.mean(axis=1); df=df.loc[~((rm==0)|(rm==1)),:]
print(f"[filter] genes x cells={df.shape} t={time.time()-t0:.0f}s",flush=True)
so=create_sifinet_object(df,rowfeature=True); so=quantile_thres2(so)
so=cal_coexp(so,X=so.data_thres['dt'],X_full=so.data_thres['dt'])
so=create_network(so,alpha=0.05,manual=False,least_edge_prop=0.01); so=filter_lowexp(so,t1=10,t2=0.9,t3=0.9)
adj=pd.DataFrame(np.where(np.abs(so.coexp-so.est_ms['mean'])>so.thres,np.abs(so.coexp),0.0))
adj.index=df.index; adj.columns=df.index
print(f"[adj] |V|={adj.shape[0]} |E|={int((adj.values!=0).sum()/2)} t={time.time()-t0:.0f}s",flush=True)
adj=convert_gene_ids(adj,"ensembl"); os.makedirs(os.path.dirname(a.out),exist_ok=True); adj.to_csv(a.out)
print(f"[done] -> {a.out} t={time.time()-t0:.0f}s",flush=True)
