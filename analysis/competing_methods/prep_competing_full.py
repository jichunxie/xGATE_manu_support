#!/usr/bin/env python
"""
Prep ENSEMBL inputs for the comprehensive competing-methods R script (all 5 methods incl
scGSEA/gficf which needs ensembl IDs). For FUCCI and TS fibroblast, export:
  <ds>_counts_ensembl.csv  : raw counts, genes(ensembl) x cells
  <ds>_genesets_ensembl.json: {pathway: {label, genes(ensembl)}}  (same benchmark pathways as xGATE)
"""
from xgate_paths import ROOT  # noqa: E402
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
from xGATE.utilities.pathway_analysis import get_entrez_mapping
import mygene
D=ROOT + "/data"

FUCCI_POS=["Cell cycle","DNA replication","Mismatch repair","Homologous recombination","Base excision repair",
 "Nucleotide excision repair","Fanconi anemia pathway","p53 signaling pathway","Oocyte meiosis",
 "Cellular senescence","Pyrimidine metabolism","Purine metabolism"]
FUCCI_NEG=["JAK-STAT signaling pathway","Toll-like receptor signaling pathway","Cytokine-cytokine receptor interaction",
 "B cell receptor signaling pathway","T cell receptor signaling pathway","Insulin secretion",
 "Maturity onset diabetes of the young","Bile secretion","Salivary secretion","Taste transduction"]
TS_POS=["ECM-receptor interaction","Focal adhesion","Protein digestion and absorption","PI3K-Akt signaling pathway",
 "Regulation of actin cytoskeleton","Proteoglycans in cancer","TGF-beta signaling pathway","Relaxin signaling pathway",
 "AGE-RAGE signaling pathway in diabetic complications","Hippo signaling pathway"]
TS_NEG=["Melanogenesis","Aldosterone synthesis and secretion","Insulin secretion","Fat digestion and absorption",
 "Parathyroid hormone synthesis, secretion and action","Pancreatic secretion","Carbohydrate digestion and absorption",
 "Type II diabetes mellitus","Bile secretion","Cortisol synthesis and secretion"]

cats=get_categorized_pathways(); mg=mygene.MyGeneInfo()
def genesets_ensembl(pos,neg,detected):
    gs={}
    for name,lab in [(p,"positive") for p in pos]+[(p,"negative") for p in neg]:
        raw=get_genes_in_pathway(gather_pathways_between(name,name,cats))
        ent=[g.replace("hsa:","") for g in raw]
        res=mg.querymany(ent,scopes="entrezgene",fields="ensembl.gene",species="human",verbose=False)
        ens=[]
        for r in res:
            if "ensembl" in r: ens.extend([e["gene"] for e in r["ensembl"]] if isinstance(r["ensembl"],list) else [r["ensembl"]["gene"]])
        ens=sorted(set(g for g in ens if g in detected))
        gs[name]=dict(label=lab,genes=ens,n_detected=len(ens))
    return gs

def fucci():
    c=pd.read_csv(f"{D}/fucci_u2os/GSE146773_Counts.csv.gz",index_col=0)  # cells x genes(ensembl)
    df=c.T; df.index=[str(g).split(".")[0] for g in df.index]; df=df[~df.index.duplicated(keep="first")]
    df=df.loc[(df>0).sum(axis=1)>=0.01*df.shape[1],:]
    df.to_csv(f"{D}/fucci_u2os/fucci_counts_ensembl.csv")
    gs=genesets_ensembl(FUCCI_POS,FUCCI_NEG,set(df.index))
    json.dump(gs,open(f"{D}/fucci_u2os/fucci_genesets_ensembl.json","w"))
    print("FUCCI:",df.shape,"genes x cells;",len(gs),"pathways")

def ts():
    np.random.seed(12)
    a=ad.read_h5ad("/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad",backed="r")
    idx=np.where((a.obs.cell_type=="fibroblast").values)[0]
    idx=np.sort(np.random.choice(idx,5000,replace=False))
    sub=a[idx].to_memory()
    X=sp.csr_matrix(sub.raw.X if sub.raw is not None else sub.X)
    genes=[str(g).split(".")[0] for g in (sub.raw.var_names if sub.raw is not None else sub.var_names)]
    df=pd.DataFrame(X.T.toarray(),index=genes); df=df[~df.index.duplicated(keep="first")]
    df=df.loc[(df>0).sum(axis=1)>=0.01*df.shape[1],:]
    df.to_csv(f"{D}/ts_fibroblast/ts_counts_ensembl.csv")
    gs=genesets_ensembl(TS_POS,TS_NEG,set(df.index))
    json.dump(gs,open(f"{D}/ts_fibroblast/ts_genesets_ensembl.json","w"))
    print("TS:",df.shape,"genes x cells;",len(gs),"pathways")

if __name__=="__main__":
    fucci(); ts(); print("[done]")
