#!/usr/bin/env python
"""
R1.3 batch effect, CORRECTED. The batch variable is DONOR. We ask a single question:
in each embedding space, how separable are the donors? Lower separation = the embedding
is more resilient to batch (donor) effects.

Fixes vs the previous scripts:
  * NO random-null gene sets anywhere. Points are ONLY the real pathway embeddings,
    one per (donor, pathway).
  * EVERY metric is computed with the DONOR label (not active-vs-random state).
  * xGATE topological embedding vs gene-expression SVD embedding, same procedure.

Metrics (all on the donor label): Calinski-Harabasz, silhouette, Davies-Bouldin, plus a
cross-validated donor classifier (logreg + kNN) and an embedding-spread diagnostic. Because
there is no random draw, each metric is a single deterministic value (no 30-draw averaging).

We report metrics on the RAW embeddings AND on per-dimension standardized embeddings, because
the two methods live on very different scales (the gene-expression SVD embedding is near-
degenerate) and any single normalization choice can flip a scale-sensitive metric. Direction
for donor separability: CH lower=better, silhouette lower=better, Davies-Bouldin higher=better
(worse-separated), classifier accuracy closer to chance=better.
"""
from xgate_paths import ROOT  # noqa: E402
import os, sys, pickle, numpy as np, pandas as pd
import scanpy as sc, anndata as ad, mygene
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway

OUT = ROOT + ""
BE = ROOT + "/data/batch_effect"
TS_H5 = "/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad"
TS_PKL = ROOT + "/data/ts_fibroblast/all_embedding_data_ts.pkl"

PANC_PATHWAYS = ["Autophagy", "Protein processing in endoplasmic reticulum", "mTOR signaling pathway",
    "HIF-1 signaling pathway", "AMPK signaling pathway", "PPAR signaling pathway",
    "Protein digestion and absorption", "PI3K-Akt signaling pathway", "Insulin signaling pathway",
    "Oxidative phosphorylation", "cGMP-PKG signaling pathway"]
TS_PATHWAYS = ["ECM-receptor interaction", "Focal adhesion", "Protein digestion and absorption",
    "PI3K-Akt signaling pathway", "Regulation of actin cytoskeleton", "Proteoglycans in cancer",
    "TGF-beta signaling pathway", "Relaxin signaling pathway",
    "AGE-RAGE signaling pathway in diabetic complications", "Hippo signaling pathway"]

_CATS = None
def pathway_ensembl(pname, valid=None):
    """KEGG pathway -> ensembl gene ids present in `valid` (a gene->index map), via mygene."""
    global _CATS, _MG
    if _CATS is None:
        _CATS = get_categorized_pathways(); _MG = mygene.MyGeneInfo()
    raw = get_genes_in_pathway(gather_pathways_between(pname, pname, _CATS))
    if not raw:
        return []
    res = _MG.querymany([g.replace("hsa:", "") for g in raw], scopes="entrezgene",
                        fields="ensembl.gene", species="human", verbose=False)
    ens = []
    for r in res:
        if "ensembl" in r:
            ens.extend([e["gene"] for e in r["ensembl"]] if isinstance(r["ensembl"], list)
                       else [r["ensembl"]["gene"]])
    ens = sorted(set(ens))
    return [g for g in ens if (valid is None or g in valid)]


def xgate_points(pkl_path):
    """Real pathway embeddings only (NO nulls): X (points x dim), donor labels."""
    d = pickle.load(open(pkl_path, "rb"))
    X, don = [], []
    for donor, paths in d.items():
        for _pname, dat in paths.items():
            X.append(np.asarray(dat["pathway_embedding"]).flatten())
            don.append(donor)
    k = min(len(v) for v in X)
    return np.array([v[:k] for v in X]), np.array(don), k


def svd_emb(Xsub, k):
    """Mean top-k right-singular-vector projection; always length exactly k (zero-pad if short)."""
    _, _, Vt = np.linalg.svd(Xsub, full_matrices=False)
    kk = min(k, Vt.shape[0])
    e = np.mean(Xsub @ Vt[:kk, :].T, axis=0)
    if len(e) < k:
        e = np.concatenate([e, np.zeros(k - len(e))])
    return e[:k]


def geneexpr_points_pancreas():
    a = sc.read(f"{BE}/pancreas_human.h5ad")
    fd = a[(a.obs["cell_type"] == "type B pancreatic cell") & (a.obs["disease_state"] == "Control")].copy()
    keep = fd.obs["donor_id"].value_counts(); keep = keep[keep >= 200].index
    fd = fd[fd.obs["donor_id"].isin(keep)].copy()
    Xf = fd.X.toarray() if not isinstance(fd.X, np.ndarray) else fd.X
    genes = sorted(fd.var_names.tolist()); g2i = {g: i for i, g in enumerate(genes)}
    pgd = {p: v for p in PANC_PATHWAYS if (v := pathway_ensembl(p, g2i))}
    k = min(len(v) for v in pgd.values())           # smallest gene set -> homogeneous dim
    idx = {dn: np.where(fd.obs["donor_id"] == dn)[0] for dn in keep}
    X, don = [], []
    for dn in keep:
        Xd = Xf[idx[dn], :]
        for _p, ag in pgd.items():
            X.append(svd_emb(Xd[:, [g2i[g] for g in ag]], k)); don.append(dn)
    return np.array(X), np.array(don)


def geneexpr_points_ts(k=10):
    a = ad.read_h5ad(TS_H5, backed="r")
    donors = list(pickle.load(open(TS_PKL, "rb")).keys())
    X, don = [], []
    for dn in donors:
        sel = (a.obs.cell_type == "fibroblast").values & (a.obs.donor_id == dn).values
        ii = np.where(sel)[0]
        if len(ii) == 0:
            continue
        if len(ii) > 1000:
            ii = np.sort(np.random.RandomState(42).choice(ii, 1000, replace=False))
        sub = a[ii].to_memory(); Xf = sub.X
        Xf = Xf.toarray() if not isinstance(Xf, np.ndarray) else Xf
        vn = [str(g).split(".")[0] for g in sub.var_names]; g2i = {g: i for i, g in enumerate(vn)}
        for p in TS_PATHWAYS:
            cols = [g2i[g] for g in pathway_ensembl(p, g2i)]
            if len(cols) < 5:
                continue
            X.append(svd_emb(Xf[:, cols], k)); don.append(dn)
    return np.array(X), np.array(don)


def donor_acc(X, don, knn=False):
    don = np.array(don); classes = sorted(set(don)); n = len(classes)
    counts = np.bincount([classes.index(v) for v in don]); cv = max(2, min(5, int(counts.min())))
    clf = (make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=5)) if knn else
           make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight="balanced")))
    pred = cross_val_predict(clf, X, don, cv=StratifiedKFold(cv, shuffle=True, random_state=0))
    return accuracy_score(don, pred), 1.0 / n


def metrics_block(name, method, X, don):
    spread = float(np.mean(np.std(X, axis=0)))
    Xz = StandardScaler().fit_transform(X)
    lg, chance = donor_acc(X, don); kn, _ = donor_acc(X, don, knn=True)
    row = dict(dataset=name, method=method, n_points=len(X), n_donors=len(set(don)),
               chance=round(chance, 3), spread=spread,
               CH_donor_raw=calinski_harabasz_score(X, don),
               sil_donor_raw=silhouette_score(X, don),
               DB_donor_raw=davies_bouldin_score(X, don),
               CH_donor_std=calinski_harabasz_score(Xz, don),
               sil_donor_std=silhouette_score(Xz, don),
               DB_donor_std=davies_bouldin_score(Xz, don),
               donor_logreg_acc=lg, donor_knn_acc=kn)
    return row


def main():
    rows = []
    # ---- pancreas ----
    Xn, dn_n, _ = xgate_points(f"{BE}/all_embedding_data.pkl")
    Xg, dn_g = geneexpr_points_pancreas()
    rows.append(metrics_block("Pancreas", "xGATE", Xn, dn_n))
    rows.append(metrics_block("Pancreas", "GeneExpression", Xg, dn_g))
    # ---- TS fibroblast ----
    Xn2, dn_n2, _ = xgate_points(TS_PKL)
    Xg2, dn_g2 = geneexpr_points_ts()
    rows.append(metrics_block("TS Fibroblast", "xGATE", Xn2, dn_n2))
    rows.append(metrics_block("TS Fibroblast", "GeneExpression", Xg2, dn_g2))

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/results/supp_figS3_5_batch_donor_metrics.csv", index=False)
    pd.set_option("display.width", 200, "display.max_columns", 30)
    print(df.round(3).to_string(index=False))
    print("\nDONOR-SEPARABILITY (batch effect). Lower CH_donor / lower sil_donor / higher DB_donor /")
    print("classifier acc closer to chance => donors less separable => MORE batch-resilient.")
    print("Compare xGATE vs GeneExpression within each dataset, on raw AND std, plus classifier.")
    for nm in ["Pancreas", "TS Fibroblast"]:
        x = df[(df.dataset == nm) & (df.method == "xGATE")].iloc[0]
        g = df[(df.dataset == nm) & (df.method == "GeneExpression")].iloc[0]
        print(f"\n[{nm}] spread xGATE={x.spread:.3g} vs GE={g.spread:.3g}"
              f"  (GE near-degenerate if <<1 -> metrics on GE unreliable)")
        print(f"  sil_donor raw  xGATE {x.sil_donor_raw:+.3f} vs GE {g.sil_donor_raw:+.3f}"
              f"   | std xGATE {x.sil_donor_std:+.3f} vs GE {g.sil_donor_std:+.3f}")
        print(f"  CH_donor  raw  xGATE {x.CH_donor_raw:.2f} vs GE {g.CH_donor_raw:.2f}"
              f"   | std xGATE {x.CH_donor_std:.2f} vs GE {g.CH_donor_std:.2f}")
        print(f"  donor logreg   xGATE {x.donor_logreg_acc:.3f} vs GE {g.donor_logreg_acc:.3f}"
              f" (chance {x.chance})  | knn xGATE {x.donor_knn_acc:.3f} vs GE {g.donor_knn_acc:.3f}")
    print("\n[done] -> results/supp_figS3_5_batch_donor_metrics.csv")


if __name__ == "__main__":
    main()
