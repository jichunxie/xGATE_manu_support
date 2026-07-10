#!/usr/bin/env python
"""Prep competing-method inputs for the TS fibroblast benchmark: export entrez-rowname
counts (pooled fibroblasts, same 5000-cell subsample seed as the graph) + the SAME KEGG
gene sets xGATE used."""
from xgate_paths import ROOT  # noqa: E402
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
from xGATE.utilities.pathway_analysis import get_entrez_mapping
TS="/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad"
OUT=ROOT + "/data/ts_fibroblast"
POS=["ECM-receptor interaction","Focal adhesion","Protein digestion and absorption",
 "PI3K-Akt signaling pathway","Regulation of actin cytoskeleton","Proteoglycans in cancer",
 "TGF-beta signaling pathway","Relaxin signaling pathway",
 "AGE-RAGE signaling pathway in diabetic complications","Hippo signaling pathway"]
# clean specialized negatives (verified inactive in TS fibroblasts; non-fibroblast cell-type programs)
NEG=["Melanogenesis","Aldosterone synthesis and secretion","Insulin secretion",
 "Fat digestion and absorption","Parathyroid hormone synthesis, secretion and action",
 "Pancreatic secretion","Carbohydrate digestion and absorption","Type II diabetes mellitus",
 "Bile secretion","Cortisol synthesis and secretion"]

def main():
    np.random.seed(12)
    a=ad.read_h5ad(TS,backed="r")
    idx=np.where((a.obs.cell_type=="fibroblast").values)[0]
    idx=np.sort(np.random.choice(idx,5000,replace=False))   # same as graph build
    sub=a[idx].to_memory()
    X=sp.csr_matrix(sub.raw.X if sub.raw is not None else sub.X)   # cells x genes
    genes=[str(g).split(".")[0] for g in (sub.raw.var_names if sub.raw is not None else sub.var_names)]
    df=pd.DataFrame(X.T.toarray(),index=genes)
    df=df[~df.index.duplicated(keep="first")]
    df=df.loc[(df>0).sum(axis=1)>=0.01*df.shape[1],:]
    print(f"genes x cells={df.shape}; mapping->entrez",flush=True)
    ent=get_entrez_mapping(list(df.index),"ensembl.gene"); df.index=[str(e) for e in ent]
    df=df[df.index!="None"].groupby(level=0).sum()
    df.to_csv(f"{OUT}/ts_fibroblast_counts_entrez.csv")
    cats=get_categorized_pathways(); detected=set(df.index); gs={}
    for name,lab in [(p,"positive") for p in POS]+[(p,"negative") for p in NEG]:
        ids=gather_pathways_between(name,name,cats); g=get_genes_in_pathway(ids) if ids else []
        ent=[x.replace("hsa:","") for x in g]; gs[name]=dict(label=lab,genes=ent,
            n_detected=len([x for x in ent if x in detected]))
    json.dump(gs,open(f"{OUT}/ts_fibroblast_genesets.json","w"))
    print("[done] wrote counts + genesets",flush=True)

if __name__=="__main__": main()
