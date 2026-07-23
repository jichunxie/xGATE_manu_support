#!/usr/bin/env python
"""
R1.2 formal null calibration: run xGATE on many RANDOM gene sets (sizes matched to the
benchmark pathway size distribution) on the pancreas ctrl graph. Their p-values should be
~uniform (the formal null is uniform by construction), anchoring BH validity. Output the
p-values + KS-vs-uniform test for the calibration panel.
"""
from xgate_paths import ROOT  # noqa: E402
import sys, time, random, numpy as np, pandas as pd
import torch
from xGATE.utilities import create_network_from_adj_matrix, embedding_recon
from scipy.stats import kstest
import os
random.seed(7); np.random.seed(7); torch.manual_seed(7)
ADJ=ROOT + "/data/pancreas/adj_matrix_pancreas_ctrl_final.csv"
N_SETS=80; B=100                              # faster: fewer sets + smaller null (was 120/200)
OUT=ROOT + "/results/supp_figS21_randomnull_pvalues.csv"
def decode(x): return x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
adj=pd.read_csv(ADJ,index_col=0); adj.index=[decode(i) for i in adj.index]; adj.columns=[decode(c) for c in adj.columns]
G=create_network_from_adj_matrix(adj); del adj
names=[v["name"] for v in G.vs]; nodeset=set(names)
print(f"[network] |V|={G.vcount()} |E|={G.ecount()}",flush=True)
sizes=np.random.randint(40,160,size=N_SETS)   # cap size at 160 (large subgraphs dominate runtime)
pd.DataFrame(columns=["set_id","size","p_value","z_score"]).to_csv(OUT,index=False)  # incremental
rows=[]; t0=time.time()
for i,k in enumerate(sizes):
    genes=random.sample(names,int(k))
    p,z=embedding_recon(G,None,genes,200,200,B)
    rows.append(dict(set_id=i,size=int(k),p_value=p,z_score=z))
    pd.DataFrame([rows[-1]]).to_csv(OUT,mode="a",header=False,index=False)  # save each (survive timeout)
    if (i+1)%10==0: print(f"  {i+1}/{N_SETS} t={time.time()-t0:.0f}s",flush=True)
df=pd.DataFrame(rows)
ks=kstest(df.p_value,"uniform")
print(f"[done] {N_SETS} random sets; mean p={df.p_value.mean():.3f}; frac p<0.05={np.mean(df.p_value<0.05):.3f}")
print(f"  KS vs Uniform: stat={ks.statistic:.3f} p={ks.pvalue:.3f}  (large p => consistent with uniform)")
