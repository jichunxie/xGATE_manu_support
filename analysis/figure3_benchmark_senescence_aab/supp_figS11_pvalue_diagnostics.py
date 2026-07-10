#!/usr/bin/env python
"""
R1.2 (efficient): p-value distribution + BH/BY on a curated, biologically meaningful
pathway set (NOT the full 353-pathway KEGG catalog). Positives = beta-cell pathways
(expect small p); negatives = disease/infection/cancer pathways = real INACTIVE
biological pathways (expect ~uniform p -> the calibration check the reviewer wants).

Outputs raw p, z, BH q, BY q + the labels for a faceted histogram.
"""
from xgate_paths import ROOT  # noqa: E402
import sys, time, os, random
import numpy as np, pandas as pd
import torch
from statsmodels.stats.multitest import multipletests
from xGATE.utilities import (create_network_from_adj_matrix, get_categorized_pathways,
                       gather_pathways_between, get_genes_in_pathway, embedding_recon)

ADJ = ROOT + "/data/pancreas/adj_matrix_pancreas_ctrl_final.csv"
OUT = ROOT + "/results/supp_figS11_pancreas_ctrl_focused.csv"

ACTIVE = [  # beta-cell / metabolic / housekeeping-active (expect small p)
 "Insulin signaling pathway","AMPK signaling pathway","mTOR signaling pathway",
 "FoxO signaling pathway","PPAR signaling pathway","HIF-1 signaling pathway",
 "Longevity regulating pathway","Pancreatic secretion","Insulin secretion",
 "Protein processing in endoplasmic reticulum","Oxidative phosphorylation","Ribosome",
 "PI3K-Akt signaling pathway","Rap1 signaling pathway","cGMP-PKG signaling pathway",
 "Type II diabetes mellitus","Maturity onset diabetes of the young",
 "Adipocytokine signaling pathway","Glucagon signaling pathway","Mineral absorption",
 "Antigen processing and presentation","Protein digestion and absorption",
]
INACTIVE = [  # disease/infection/cancer = real inactive biological pathways (expect ~uniform)
 "Bacterial invasion of epithelial cells","Shigellosis","Pathogenic Escherichia coli infection",
 "Salmonella infection","Vibrio cholerae infection","Legionellosis","Pertussis","Yersinia infection",
 "Tuberculosis","Leishmaniasis","Amoebiasis","Malaria","Toxoplasmosis","Chagas disease",
 "African trypanosomiasis","Influenza A","Measles","Hepatitis B","Hepatitis C",
 "Herpes simplex virus 1 infection","Human papillomavirus infection","Epstein-Barr virus infection",
 "Colorectal cancer","Pancreatic cancer","Hepatocellular carcinoma","Gastric cancer","Glioma",
 "Thyroid cancer","Acute myeloid leukemia","Chronic myeloid leukemia","Basal cell carcinoma",
 "Melanoma","Renal cell carcinoma","Bladder cancer","Prostate cancer","Endometrial cancer",
 "Breast cancer","Small cell lung cancer","Non-small cell lung cancer","Rheumatoid arthritis",
 "Systemic lupus erythematosus","Asthma","Inflammatory bowel disease","Alzheimer disease",
 "Parkinson disease","Huntington disease","Amyotrophic lateral sclerosis","Prion disease",
]

def decode(x): return x.decode() if isinstance(x,(bytes,bytearray)) else str(x)

def main():
    random.seed(12); np.random.seed(12); torch.manual_seed(12)
    t0 = time.time()
    adj = pd.read_csv(ADJ, index_col=0)
    adj.index=[decode(i) for i in adj.index]; adj.columns=[decode(c) for c in adj.columns]
    G = create_network_from_adj_matrix(adj); del adj
    node_set = set(v["name"] for v in G.vs)
    cats = get_categorized_pathways()
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    rows=[]
    for name, tr in [(p,"active") for p in ACTIVE]+[(p,"inactive") for p in INACTIVE]:
        ids = gather_pathways_between(name, name, cats)
        genes = get_genes_in_pathway(ids) if ids else []
        det = [g for g in genes if g in node_set]
        if len(det) < 15:
            print(f"  [skip] {name} ({len(det)} detected)", flush=True); continue
        p,z = embedding_recon(G, cats, genes, 200,200,200)
        rows.append(dict(pathway=name, label=tr, n_detected=len(det), p_value=p, z_score=z))
        print(f"  [{tr:8s}] {name:42s} p={p:.4f} z={z:8.2f} n={len(det)} t={time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    # multiple-testing corrections over the tested set
    p = df.p_value.values
    df["q_BH"] = multipletests(p, method="fdr_bh")[1]
    df["q_BY"] = multipletests(p, method="fdr_by")[1]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\n[done] {len(df)} pathways t={time.time()-t0:.0f}s -> {OUT}", flush=True)
    for lab in ["active","inactive"]:
        s=df[df.label==lab]
        print(f"  {lab}: n={len(s)} median p={s.p_value.median():.3f} "
              f"frac p<0.05={np.mean(s.p_value<0.05):.2f} BH-sig={np.sum(s.q_BH<0.05)} BY-sig={np.sum(s.q_BY<0.05)}")

if __name__ == "__main__":
    main()
