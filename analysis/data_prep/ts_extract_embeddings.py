#!/usr/bin/env python
"""
TS fibroblast batch (R1.3 2nd dataset), step 2: from the 7 per-donor co-expression graphs,
extract xGATE graph embeddings (pathway_embedding + null_embeddings) for the active fibroblast
pathways, saving a pkl in the same structure as the pancreas all_embedding_data.pkl.
Run after the per-donor graphs (donor_adj/adj_TSP*.csv) are built.
"""
from xgate_paths import ROOT  # noqa: E402
import sys, glob, os, random, pickle, numpy as np, pandas as pd
import torch
from xGATE.utilities import create_network_from_adj_matrix, get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
from xGATE.utilities.embeddings import generate_embedding
random.seed(12); np.random.seed(12); torch.manual_seed(12)
DA=ROOT + "/data/ts_fibroblast/donor_adj"; B=200
ACTIVE=["ECM-receptor interaction","Focal adhesion","Protein digestion and absorption","PI3K-Akt signaling pathway",
 "Regulation of actin cytoskeleton","Proteoglycans in cancer","TGF-beta signaling pathway","Relaxin signaling pathway",
 "AGE-RAGE signaling pathway in diabetic complications","Hippo signaling pathway"]
def decode(x): return x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
cats=get_categorized_pathways()
pgenes={p:get_genes_in_pathway(gather_pathways_between(p,p,cats)) for p in ACTIVE}
out={}
for adjf in sorted(glob.glob(f"{DA}/adj_*.csv")):
    if adjf.endswith("_timing.csv"): continue          # skip timing summaries
    donor=os.path.basename(adjf).replace("adj_","").replace(".csv","")
    adj=pd.read_csv(adjf,index_col=0); adj.index=[decode(i) for i in adj.index]; adj.columns=[decode(c) for c in adj.columns]
    if adj.shape[0]!=adj.shape[1] or adj.shape[0]<500:
        print(f"[skip] {donor}: non-square/small adj {adj.shape}",flush=True); continue
    G=create_network_from_adj_matrix(adj); del adj
    if G.ecount()<10000:                                # drop degenerate graphs (e.g. TSP14)
        print(f"[skip] {donor}: degenerate graph |E|={G.ecount()}",flush=True); continue
    names=[v["name"] for v in G.vs]; nodeset=set(names)
    out[donor]={}
    for p,genes in pgenes.items():
        found=[g for g in genes if g in nodeset]
        if len(found)<5: continue
        sub=G.subgraph([v.index for v in G.vs if v["name"] in set(found)]); k=sub.vcount()
        pe=generate_embedding(sub)
        nulls=[]
        for _ in range(B):
            verts=set(random.sample(names,k)); nulls.append(generate_embedding(G.subgraph([v.index for v in G.vs if v["name"] in verts])))
        out[donor][p]={"pathway_embedding":np.array(pe),"null_embeddings":[np.array(n) for n in nulls],"size":k}
    print(f"[{donor}] |V|={G.vcount()} pathways={len(out[donor])}",flush=True)
pickle.dump(out,open(ROOT + "/data/ts_fibroblast/all_embedding_data_ts.pkl","wb"))
print(f"[done] {len(out)} donors -> all_embedding_data_ts.pkl")
