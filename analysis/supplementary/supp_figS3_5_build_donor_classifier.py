#!/usr/bin/env python
"""
R1.3 settle: full-dim donor-CLASSIFIER, xGATE topological embeddings vs gene-expression
SVD embeddings (control pancreatic beta cells, 5 donors). Builds the SAME embeddings as
the batch-effect / donor-silhouette analysis (see supp_figS3_5_batch_figures.py) but asks: how
recoverable is donor identity from
each embedding? Lower accuracy = donor less recoverable = better donor mixing. This checks
whether the donor-classifier agrees or disagrees with the reported donor silhouette
(xGATE -0.23 vs gene-expr -0.14), and whether the gene-expr embedding is degenerate.
"""
from xgate_paths import ROOT  # noqa: E402
import os, sys, pickle, numpy as np, pandas as pd
import scanpy as sc, mygene
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway

BE = ROOT + "/data/batch_effect"; OUT = ROOT + ""

# --- xGATE active embeddings per (donor, pathway) ---
d = pickle.load(open(f"{BE}/all_embedding_data.pkl", "rb"))
net_active = []; net_don = []
for donor, paths in d.items():
    for pname, dat in paths.items():
        net_active.append(np.array(dat["pathway_embedding"]).flatten()); net_don.append(donor)
net_don = np.array(net_don)

# --- gene-expression SVD embeddings (Suppl. Note 3 procedure) ---
ad_ = sc.read(f"{BE}/pancreas_human.h5ad")
fd = ad_[(ad_.obs['cell_type'] == 'type B pancreatic cell') & (ad_.obs['disease_state'] == 'Control')].copy()
keep = fd.obs['donor_id'].value_counts(); keep = keep[keep >= 200].index
fd = fd[fd.obs['donor_id'].isin(keep)].copy()
Xf = fd.X.toarray() if not isinstance(fd.X, np.ndarray) else fd.X
genes = sorted(fd.var_names.tolist()); g2i = {g: i for i, g in enumerate(genes)}
test_pathways = ["Autophagy", "Protein processing in endoplasmic reticulum", "mTOR signaling pathway",
 "HIF-1 signaling pathway", "AMPK signaling pathway", "PPAR signaling pathway", "Protein digestion and absorption",
 "PI3K-Akt signaling pathway", "Insulin signaling pathway", "Oxidative phosphorylation", "cGMP-PKG signaling pathway"]
cats = get_categorized_pathways(); mg = mygene.MyGeneInfo(); pgd = {}
for pn in test_pathways:
    raw = get_genes_in_pathway(gather_pathways_between(pn, pn, cats))
    if not raw: continue
    res = mg.querymany([g.replace('hsa:', '') for g in raw], scopes='entrezgene', fields='ensembl.gene', species='human', verbose=False)
    ens = []
    for r in res:
        if 'ensembl' in r: ens.extend([e['gene'] for e in r['ensembl']] if isinstance(r['ensembl'], list) else [r['ensembl']['gene']])
    v = sorted(set(g for g in ens if g in g2i))
    if v: pgd[pn] = v
k = min(len(v) for v in pgd.values())
def svd_emb(Xsub):
    _, _, Vt = np.linalg.svd(Xsub, full_matrices=False); return np.mean(Xsub @ Vt[:k, :].T, axis=0)
sc_active = []; sc_don = []
donor_idx = {dn: np.where(fd.obs['donor_id'] == dn)[0] for dn in keep}
for dn in keep:
    Xd = Xf[donor_idx[dn], :]
    for pn, ag in pgd.items():
        sc_active.append(svd_emb(Xd[:, [g2i[g] for g in ag]])); sc_don.append(dn)
sc_don = np.array(sc_don)
net_active = np.array([v[:k] for v in net_active]); sc_active = np.array(sc_active)

def donor_acc(X, don, knn=False):
    don = np.array(don); classes = sorted(set(don)); n = len(classes); chance = 1.0 / n
    counts = np.bincount([classes.index(v) for v in don]); cv = max(2, min(5, int(counts.min())))
    if knn:
        clf = make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=5))
    else:
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight='balanced'))
    pred = cross_val_predict(clf, X, don, cv=StratifiedKFold(cv, shuffle=True, random_state=0))
    return accuracy_score(don, pred), chance, n

axl, ch, n = donor_acc(net_active, net_don); agl, _, _ = donor_acc(sc_active, sc_don)
axk, _, _ = donor_acc(net_active, net_don, knn=True); agk, _, _ = donor_acc(sc_active, sc_don, knn=True)
# degeneracy diagnostic: how much variance does each embedding carry?
def spread(X): return float(np.mean(np.std(X, axis=0)))
print(f"[R1.3 donor-classifier full-dim k={k}] pancreas control, {n} donors, chance={ch:.3f}")
print(f"  logreg: xGATE acc={axl:.3f}  gene-expr acc={agl:.3f}")
print(f"  kNN   : xGATE acc={axk:.3f}  gene-expr acc={agk:.3f}")
print(f"  embedding spread (mean per-dim SD): xGATE={spread(net_active):.4g}  gene-expr={spread(sc_active):.4g}")
print(f"  rows: xGATE={len(net_don)} gene-expr={len(sc_don)}")
print("  lower acc = donor LESS recoverable = better mixing; compare to donor silhouette xGATE -0.23 vs GE -0.14")
pd.DataFrame([
    {"emb": "xGATE", "logreg_acc": axl, "knn_acc": axk, "chance": ch, "spread": spread(net_active)},
    {"emb": "GeneExpression", "logreg_acc": agl, "knn_acc": agk, "chance": ch, "spread": spread(sc_active)},
]).to_csv(f"{OUT}/results/supp_figS3_5_donor_classifier_pancreas.csv", index=False)
print("saved -> results/supp_figS3_5_donor_classifier_pancreas.csv")
