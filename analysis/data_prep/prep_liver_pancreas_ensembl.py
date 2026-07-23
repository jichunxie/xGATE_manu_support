#!/usr/bin/env python
"""
Prep ENSEMBL inputs so the SAME run_all_competing_metrics.R can be run on liver + pancreas
(for a consistent 4-dataset, 6-method Fig 3a). Generates gene sets (ensembl) for both, and
the PANCREAS count matrix (ctrl beta cells) from pancreas_human.h5ad. Liver counts come from
the .rds via extract_liver_counts.R (run in RStudio; rownames must be ensembl).
"""
from xgate_paths import ROOT  # noqa: E402
import sys, json, numpy as np, pandas as pd
import anndata as ad, scanpy as sc, mygene
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
D=ROOT + "/data"

# manuscript benchmark pathway sets (Supp Fig 3 pancreas; liver-benchmark figure)
PANC_POS=["AMPK signaling pathway","Autophagy","HIF-1 signaling pathway","Insulin signaling pathway",
 "Oxidative phosphorylation","PI3K-Akt signaling pathway","PPAR signaling pathway",
 "Protein digestion and absorption","Protein processing in endoplasmic reticulum",
 "cGMP-PKG signaling pathway","mTOR signaling pathway"]
PANC_NEG=["Apoptosis","Bacterial invasion of epithelial cells","Cytokine-cytokine receptor interaction",
 "Pathogenic Escherichia coli infection","Salmonella infection","Shigellosis","Thyroid cancer"]
LIVER_POS=["Carbon metabolism","Biosynthesis of amino acids","Biosynthesis of cofactors","Glycolysis / Gluconeogenesis",
 "Pentose and glucuronate interconversions","Ascorbate and aldarate metabolism","Pyruvate metabolism","Fatty acid degradation",
 "Primary bile acid biosynthesis","Steroid hormone biosynthesis","Arachidonic acid metabolism","Linoleic acid metabolism",
 "Biosynthesis of unsaturated fatty acids","Glycine, serine and threonine metabolism","Cysteine and methionine metabolism",
 "Tyrosine metabolism","Taurine and hypotaurine metabolism","Retinol metabolism","Porphyrin metabolism",
 "Metabolism of xenobiotics by cytochrome P450","Cholesterol metabolism","Caffeine metabolism","Drug metabolism"]
LIVER_NEG=["Thyroid cancer","Shigellosis","Colorectal cancer","Pancreatic cancer","Hepatocellular carcinoma","Gastric cancer",
 "Glioma","Acute myeloid leukemia","Chronic myeloid leukemia","Basal cell carcinoma","Melanoma","Renal cell carcinoma",
 "Bladder cancer","Prostate cancer","Endometrial cancer","Breast cancer","Small cell lung cancer","Non-small cell lung cancer"]

cats=get_categorized_pathways(); mg=mygene.MyGeneInfo()
def gsets(pos,neg,detected=None):
    out={}
    for name,lab in [(p,"positive") for p in pos]+[(p,"negative") for p in neg]:
        ids=gather_pathways_between(name,name,cats)
        if not ids: print("  [warn] no KEGG id:",name); continue
        raw=get_genes_in_pathway(ids); ent=[g.replace("hsa:","") for g in raw]
        res=mg.querymany(ent,scopes="entrezgene",fields="ensembl.gene",species="human",verbose=False)
        ens=[]
        for r in res:
            if "ensembl" in r: ens.extend([e["gene"] for e in r["ensembl"]] if isinstance(r["ensembl"],list) else [r["ensembl"]["gene"]])
        ens=sorted(set(ens if detected is None else [g for g in ens if g in detected]))
        out[name]=dict(label=lab,genes=ens,n_detected=len(ens))
    return out

# pancreas counts (ctrl beta cells, ensembl)
a=sc.read(f"{D}/pancreas/pancreas_human.h5ad")
fd=a[(a.obs.get('cell_type')=='type B pancreatic cell')&(a.obs.get('disease_state')=='Control')].copy() \
   if 'disease_state' in a.obs else a[a.obs.get('cell_type')=='type B pancreatic cell'].copy()
X=fd.X.toarray() if not isinstance(fd.X,np.ndarray) else fd.X
df=pd.DataFrame(X.T, index=[str(g).split('.')[0] for g in fd.var_names])
df=df[~df.index.duplicated(keep="first")]; df=df.loc[(df>0).sum(axis=1)>=0.01*df.shape[1],:]
df.to_csv(f"{D}/pancreas/pancreas_counts_ensembl.csv")
json.dump(gsets(PANC_POS,PANC_NEG,set(df.index)),open(f"{D}/pancreas/pancreas_genesets_ensembl.json","w"))
print(f"pancreas: counts {df.shape}; genesets written")
# liver genesets (counts extracted separately in R from the .rds)
json.dump(gsets(LIVER_POS,LIVER_NEG,None),open(f"{D}/liver/liver_genesets_ensembl.json","w"))
print("liver: genesets written (counts via extract_liver_counts.R)")
print("[done]")
