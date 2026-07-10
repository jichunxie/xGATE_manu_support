#!/usr/bin/env python
"""
R1.3 batch effect, MANY-GENE-SET donor analysis (per user 2026-07-02).

Batch variable = DONOR. We embed the SAME panel of REAL KEGG pathways with BOTH the xGATE
co-expression-graph embedding AND the matched gene-expression SVD embedding, per donor, then
ask how separable the DONORS are in each embedding space.

Key points (per user 2026-07-02):
  * The panel is REAL KEGG pathways (>= MIN_PATH_GENES detected genes), NOT random gene sets.
    There are no random gene sets and no nulls anywhere in this analysis.
  * The SAME pathway panel is used for xGATE and for gene expression, per donor.
  * EVERY metric uses the DONOR label. No active/inactive contrast, no nulls-as-contrast.
  * Deterministic (fixed pathway panel).

xGATE embedding  = generate_embedding(donor_graph.subgraph(set))  (topology descriptor).
Gene-expr embed. = mean top-k SVD projection of the donor's log-normalized expression on
                   the same genes (Suppl. Note 3 procedure).
Both truncated to a shared dimension k = min set size.

Usage: supp_figS3_5_build_batch_manyset.py --dataset {pancreas,ts} [--n_sets 200] [--seed 0]
Writes results/supp_figS3_5_batch_manyset_<dataset>.csv (one row per method).
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, os, pickle, random, json
import numpy as np, pandas as pd
import igraph as ig
from xGATE.utilities.embeddings import generate_embedding
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score

OUT = ROOT + ""
PANC_H5 = ROOT + "/data/batch_effect/pancreas_human.h5ad"
PANC_ADJ = ROOT + "/data/pancreas/donor_adj"
PANC_DONORS = ["HPAP022", "HPAP035", "HPAP036", "HPAP037", "HPAP040"]
TS_H5 = "/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad"
TS_ADJ = ROOT + "/data/ts_fibroblast/donor_adj"
TS_DONORS = ["TSP1", "TSP2", "TSP4", "TSP6", "TSP7", "TSP10", "TSP14"]
KEGG_PANEL = ROOT + "/data/kegg_pathways_hsa.json"  # {name: [hsa:entrez]}
MIN_PATH_GENES = 20               # keep pathways with >= this many detected genes; = embed dim


def load_graph(path):
    adf = pd.read_csv(path, index_col=0)
    names = [str(g).split(".")[0] for g in adf.index]
    A = adf.values.astype(float)
    iu = np.triu_indices_from(A, k=1)
    w = A[iu]
    nz = w != 0
    edges = list(zip(iu[0][nz].tolist(), iu[1][nz].tolist()))
    G = ig.Graph(n=len(names), edges=edges)
    G.vs["name"] = names
    G.es["weight"] = w[nz].tolist()
    return G, names


def emb_xgate(G, name2idx, geneset, k):
    idx = [name2idx[g] for g in geneset if g in name2idx]
    if len(idx) < 2:
        return None
    sub = G.subgraph(idx)
    e = np.asarray(generate_embedding(sub), dtype=float)
    if e.size < k:
        e = np.concatenate([e, np.zeros(k - e.size)])
    return e[:k]


def emb_geneexpr(Xdonor, g2i, geneset, k):
    cols = [g2i[g] for g in geneset if g in g2i]
    if len(cols) < 2:
        return None
    Xsub = Xdonor[:, cols]
    _, _, Vt = np.linalg.svd(Xsub, full_matrices=False)
    kk = min(k, Vt.shape[0])
    e = np.mean(Xsub @ Vt[:kk, :].T, axis=0)
    if e.size < k:
        e = np.concatenate([e, np.zeros(k - e.size)])
    return e[:k]


def donor_acc(X, don, knn=False):
    don = np.array(don); classes = sorted(set(don)); n = len(classes)
    counts = np.bincount([classes.index(v) for v in don]); cv = max(2, min(5, int(counts.min())))
    clf = (make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=5)) if knn else
           make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, class_weight="balanced")))
    pred = cross_val_predict(clf, X, don, cv=StratifiedKFold(cv, shuffle=True, random_state=0))
    return accuracy_score(don, pred), 1.0 / n


def metrics_block(dataset, method, X, don):
    X = np.asarray(X); spread = float(np.mean(np.std(X, axis=0)))
    Xz = StandardScaler().fit_transform(X)
    lg, chance = donor_acc(X, don); kn, _ = donor_acc(X, don, knn=True)
    return dict(dataset=dataset, method=method, n_points=len(X), n_donors=len(set(don)),
                chance=round(chance, 3), spread=spread,
                CH_donor_raw=calinski_harabasz_score(X, don),
                sil_donor_raw=silhouette_score(X, don),
                DB_donor_raw=davies_bouldin_score(X, don),
                CH_donor_std=calinski_harabasz_score(Xz, don),
                sil_donor_std=silhouette_score(Xz, don),
                DB_donor_std=davies_bouldin_score(Xz, don),
                donor_logreg_acc=lg, donor_knn_acc=kn)


def entrez_to_ensembl(hsa_names):
    import mygene
    mg = mygene.MyGeneInfo()
    ent = [n.replace("hsa:", "") for n in hsa_names]
    res = mg.querymany(ent, scopes="entrezgene", fields="ensembl.gene", species="human", verbose=False)
    m = {}
    for r in res:
        if "ensembl" in r:
            g = r["ensembl"]
            ens = g[0]["gene"] if isinstance(g, list) else g["gene"]
            m["hsa:" + str(r["query"])] = str(ens).split(".")[0]
    return m


def run(dataset, n_sets, seed):
    import scanpy as sc, anndata as ad
    rng = random.Random(seed)
    if dataset == "pancreas":
        donors = PANC_DONORS
        graphs = {d: load_graph(f"{PANC_ADJ}/adj_{d}.csv") for d in donors}
        # graph nodes are ensembl; expression is ensembl -> one ID space
        a = sc.read(PANC_H5)
        fd = a[(a.obs["cell_type"] == "type B pancreatic cell") & (a.obs["disease_state"] == "Control")].copy()
        fd = fd[fd.obs["donor_id"].isin(donors)].copy()
        Xf = fd.X.toarray() if not isinstance(fd.X, np.ndarray) else fd.X
        expr_genes = {str(g).split(".")[0]: i for i, g in enumerate(fd.var_names)}
        donor_expr = {d: (Xf[np.where(fd.obs["donor_id"].values == d)[0], :], expr_genes) for d in donors}
        graph_id = {d: {n: i for i, n in enumerate(graphs[d][1])} for d in donors}
        universe = set.intersection(*[set(graphs[d][1]) for d in donors]) & set(expr_genes)
    else:  # ts
        donors = TS_DONORS
        graphs = {d: load_graph(f"{TS_ADJ}/adj_{d}.csv") for d in donors}
        # graph nodes are hsa:entrez -> map to ensembl to align with expression
        all_hsa = set.union(*[set(pd.read_csv(f"{TS_ADJ}/adj_{d}.csv", index_col=0, nrows=0).columns)
                              for d in donors])
        h2e = entrez_to_ensembl(sorted(all_hsa))
        a = ad.read_h5ad(TS_H5, backed="r")
        expr_all = {str(g).split(".")[0]: i for i, g in enumerate(a.var_names)}
        donor_expr = {}
        for d in donors:
            sel = (a.obs.cell_type == "fibroblast").values & (a.obs.donor_id == d).values
            ii = np.where(sel)[0]
            if len(ii) > 1500:
                ii = np.sort(np.random.RandomState(42).choice(ii, 1500, replace=False))
            sub = a[ii].to_memory(); X = sub.X
            X = X.toarray() if not isinstance(X, np.ndarray) else X
            g2i = {str(g).split(".")[0]: i for i, g in enumerate(sub.var_names)}
            donor_expr[d] = (X, g2i)
        graph_id = {d: {n: i for i, n in enumerate(graphs[d][1])} for d in donors}
        # universe of ENSEMBL ids present as graph node (via map) in all donors AND in expression
        def ens_nodes(d):
            return {h2e[h] for h in graphs[d][1] if ("hsa:" + h.split(":")[-1]) in h2e or h in h2e}
        # graphs[d][1] already hsa:entrez strings; map directly
        per_donor_ens = []
        e2h = {}
        for d in donors:
            s = set()
            for h in graphs[d][1]:
                if h in h2e:
                    s.add(h2e[h]); e2h.setdefault(d, {})[h2e[h]] = h
            per_donor_ens.append(s)
        universe = set.intersection(*per_donor_ens) & set(expr_all)

    universe = sorted(universe); universe_set = set(universe)
    print(f"[{dataset}] shared gene universe = {len(universe)} genes; {len(donors)} donors", flush=True)
    # REAL KEGG pathways (no random gene sets, no nulls): map each pathway's genes into the
    # shared ensembl universe; keep pathways with >= MIN_PATH_GENES detected genes. The SAME
    # pathway panel is embedded by both methods for every donor; the only label is the donor.
    panel = json.load(open(KEGG_PANEL))                 # {name: [hsa:entrez]}
    all_ent = sorted({g for genes in panel.values() for g in genes})
    ent2ens = entrez_to_ensembl(all_ent)                # {"hsa:E": ensembl}
    gene_sets, set_names = [], []
    for name, genes in panel.items():
        ens = {ent2ens[g] for g in genes if g in ent2ens} & universe_set
        if len(ens) >= MIN_PATH_GENES:
            gene_sets.append(sorted(ens)); set_names.append(name)
    print(f"[{dataset}] real pathways used = {len(gene_sets)} (>= {MIN_PATH_GENES} detected genes)", flush=True)
    k = MIN_PATH_GENES  # shared embedding dim

    rows = []
    Xx, Xe, don = [], [], []
    for d in donors:
        G, names = graphs[d]; nid = graph_id[d]
        Xdon, g2i = donor_expr[d]
        if dataset == "ts":
            e2h_d = e2h[d]
        for gs in gene_sets:
            if dataset == "pancreas":
                gs_graph = gs                       # ensembl node names in graph
                gs_expr = gs
            else:
                gs_graph = [e2h_d[g] for g in gs if g in e2h_d]  # back to hsa for graph
                gs_expr = gs
            ex = emb_xgate(G, nid, gs_graph, k)
            eg = emb_geneexpr(Xdon, g2i, gs_expr, k)
            if ex is None or eg is None:
                continue
            Xx.append(ex); Xe.append(eg); don.append(d)
    print(f"[{dataset}] usable points per method = {len(don)}", flush=True)
    # save the per-point embeddings + donor labels so the same-latent-space scatter can be
    # drawn without re-running (each row = one (donor, gene-set) point, aligned across methods).
    np.savez(f"{OUT}/results/supp_figS3_5_batch_manyset_{dataset}_embeddings.npz",
             xgate=np.asarray(Xx, dtype=float), geneexpr=np.asarray(Xe, dtype=float),
             donor=np.asarray(don))
    rows.append(metrics_block(dataset, "xGATE", Xx, don))
    rows.append(metrics_block(dataset, "GeneExpression", Xe, don))
    df = pd.DataFrame(rows)
    outp = f"{OUT}/results/supp_figS3_5_batch_manyset_{dataset}.csv"
    df.to_csv(outp, index=False)
    pd.set_option("display.width", 200, "display.max_columns", 30)
    print(df.round(3).to_string(index=False))
    x = df[df.method == "xGATE"].iloc[0]; g = df[df.method == "GeneExpression"].iloc[0]
    print(f"\n[{dataset}] donor separability (lower = more batch-resilient):")
    print(f"  spread   xGATE {x.spread:.3g}  vs GE {g.spread:.3g}")
    print(f"  CH_donor std  xGATE {x.CH_donor_std:.2f} vs GE {g.CH_donor_std:.2f}")
    print(f"  sil_donor std xGATE {x.sil_donor_std:+.3f} vs GE {g.sil_donor_std:+.3f}")
    print(f"  donor logreg  xGATE {x.donor_logreg_acc:.3f} vs GE {g.donor_logreg_acc:.3f} (chance {x.chance})")
    print(f"[done] -> {outp}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["pancreas", "ts"], required=True)
    ap.add_argument("--n_sets", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    run(args.dataset, args.n_sets, args.seed)
