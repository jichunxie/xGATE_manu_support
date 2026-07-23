#!/usr/bin/env python
"""
Prep for competing-method comparison on FUCCI: export (1) an entrez-rowname count
matrix (genes x cells) and (2) the SAME KEGG pathway gene sets (entrez) xGATE used,
with positive/negative labels. Guarantees identical gene IDs across methods.
"""
from xgate_paths import ROOT  # noqa: E402
import sys, json, numpy as np, pandas as pd
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway
from xGATE.utilities.pathway_analysis import get_entrez_mapping

D = ROOT + "/data/fucci_u2os"
OUT_COUNTS = f"{D}/fucci_counts_entrez.csv"
OUT_GS = f"{D}/fucci_genesets.json"

POS = ["Cell cycle","DNA replication","Mismatch repair","Homologous recombination",
       "Base excision repair","Nucleotide excision repair","Fanconi anemia pathway",
       "p53 signaling pathway","Oocyte meiosis","Cellular senescence","Pyrimidine metabolism",
       "Purine metabolism"]
NEG = ["JAK-STAT signaling pathway","Toll-like receptor signaling pathway",
       "Cytokine-cytokine receptor interaction","B cell receptor signaling pathway",
       "T cell receptor signaling pathway","Insulin secretion",
       "Maturity onset diabetes of the young","Bile secretion","Salivary secretion",
       "Taste transduction"]   # clean negatives (immune + endocrine/secretory)

def main():
    print("[load] FUCCI counts ...", flush=True)
    c = pd.read_csv(f"{D}/GSE146773_Counts.csv.gz", index_col=0)  # cells x genes (ensembl)
    df = c.T                                                       # genes x cells
    df.index = [str(g).split(".")[0] for g in df.index]
    df = df[~df.index.duplicated(keep="first")]
    # keep genes detected in >=1% cells (AUCell ranking needs a real matrix, drop all-zero)
    df = df.loc[(df > 0).sum(axis=1) >= 0.01 * df.shape[1], :]
    print(f"  genes x cells = {df.shape}; mapping ensembl->entrez ...", flush=True)
    ent = get_entrez_mapping(list(df.index), "ensembl.gene")
    df.index = [str(e) for e in ent]
    df = df[df.index != "None"]
    df = df.groupby(level=0).sum()       # collapse dup entrez (sum counts)
    df.to_csv(OUT_COUNTS)
    print(f"  wrote {OUT_COUNTS} {df.shape}", flush=True)

    cats = get_categorized_pathways()
    detected = set(df.index)
    gs = {}
    for name, lab in [(p,"positive") for p in POS]+[(p,"negative") for p in NEG]:
        ids = gather_pathways_between(name, name, cats)
        genes = get_genes_in_pathway(ids) if ids else []
        ent = [g.replace("hsa:","") for g in genes]               # entrez
        det = [g for g in ent if g in detected]
        gs[name] = dict(label=lab, genes=ent, n_detected=len(det))
        print(f"  {lab:8s} {name:40s} detected {len(det)}/{len(ent)}", flush=True)
    json.dump(gs, open(OUT_GS,"w"))
    print(f"[done] wrote {OUT_GS}", flush=True)

if __name__ == "__main__":
    main()
