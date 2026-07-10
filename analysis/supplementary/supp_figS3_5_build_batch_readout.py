#!/usr/bin/env python
"""
R1.3 batch: READ-OUT variance decomposition on the BENCHMARK ground-truth pathways
(per user 2026-07-03).

Motivation. Suppl. Fig. 1c decomposes the variance of the xGATE pathway read-out into a
"state" factor and a "donor" factor, but there the state factor was active-vs-random-null
(the pkl only held positive pathways) and only xGATE was shown. To ground the state factor in
real biology and to compare against gene expression, we RESTRICT the tested pathways to the
SAME set used in the method benchmark (positives = active, negatives = inactive; ground truth
from xgate_paths import ROOT  # noqa: E402
from `pancreas_xgate_bench.csv` / `ts_fibroblast_xgate_sct.csv`). For every (donor, benchmark
pathway) we compute a scalar read-out by the canonical xGATE procedure on that DONOR's
co-expression graph, and a MATCHED gene-expression read-out on the same donor's expression,
then run the same two-factor ANOVA  z ~ C(state) + C(donor)  per method. The state:donor split
is the signal(state)-to-noise(donor) ratio; we compare it xGATE vs gene expression.

xGATE read-out  = canonical `embedding_recon`: pathway subgraph embedding -> VAE reconstruction
                  error vs size-matched random-subgraph null -> z (utilities/embeddings.py).
Gene-expr read-out = the SAME VAE machinery applied to the SVD expression embedding of the
                  pathway genes vs random size-matched gene sets (the Suppl. Note 3 embedding,
                  reduced to a per-(donor,pathway) scalar so it is apples-to-apples with xGATE).

The nulls here are INTERNAL to the score (exactly as in the published read-out and in the
Suppl. Fig. 1c decomposition); the STATE label is ground truth, not a null contrast.

Gene definitions match the benchmark exactly (KEGG via get_genes_in_pathway; verified to
reproduce the benchmark n_genes). Pancreas graph nodes are ENSEMBL, TS graph nodes are
hsa:entrez; expression is ENSEMBL in both, so hsa->ensembl is mapped where needed via mygene.

Usage:
  supp_figS3_5_build_batch_readout.py --dataset {pancreas,ts} --donor <ID>  [--B 100] [--seed 0]
  supp_figS3_5_build_batch_readout.py --dataset pancreas --donor HPAP022 --timing   # score 2 pathways, print timing
Writes results/supp_figS3_5_readout_<dataset>_<donor>.csv (one row per pathway: state, z_xgate, z_ge, ...).
"""
import argparse, sys, time, json
import numpy as np, pandas as pd
import random
import igraph as ig
import torch
from xGATE.utilities.embeddings import generate_embedding, embedding_recon
from xGATE.utilities.pathway_analysis import normalize
from xGATE.utilities.vae_model import VariationalAutoencoder, vae_loss_function, calculate_reconstruction_error
from xGATE.utilities import get_categorized_pathways, gather_pathways_between, get_genes_in_pathway

OUT = ROOT + ""
PANC_H5 = ROOT + "/data/batch_effect/pancreas_human.h5ad"
PANC_ADJ = ROOT + "/data/pancreas/donor_adj"
TS_H5 = "/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad"
TS_ADJ = ROOT + "/data/ts_fibroblast/donor_adj"
BENCH = {"pancreas": f"{OUT}/results/pancreas_xgate_bench.csv",
         "ts": f"{OUT}/results/ts_fibroblast_xgate_sct.csv"}
EMB_K = 20      # fixed dim for the saved pathway embeddings (scatter panels), matches supp_figS3_5_batch_manyset

_CATS = None
_MG = None


def kegg_hsa_genes(pname):
    """Benchmark gene source: KEGG hsa:entrez members of a pathway (dedup, matches benchmark n_genes)."""
    global _CATS
    if _CATS is None:
        _CATS = get_categorized_pathways()
    raw = get_genes_in_pathway(gather_pathways_between(pname, pname, _CATS))
    return sorted(set(raw))


def hsa_to_ensembl(hsa_names):
    global _MG
    if _MG is None:
        import mygene
        _MG = mygene.MyGeneInfo()
    ent = [n.replace("hsa:", "") for n in hsa_names]
    res = _MG.querymany(ent, scopes="entrezgene", fields="ensembl.gene", species="human", verbose=False)
    m = {}
    for r in res:
        if "ensembl" in r:
            g = r["ensembl"]
            ens = [e["gene"] for e in g] if isinstance(g, list) else [g["gene"]]
            m["hsa:" + str(r["query"])] = [str(e).split(".")[0] for e in ens]
    return m


def load_graph(path):
    adf = pd.read_csv(path, index_col=0)
    names = [str(g).split(".")[0] for g in adf.index]
    A = adf.values.astype(float)
    iu = np.triu_indices_from(A, k=1)
    w = A[iu]; nz = w != 0
    edges = list(zip(iu[0][nz].tolist(), iu[1][nz].tolist()))
    G = ig.Graph(n=len(names), edges=edges)
    G.vs["name"] = names
    G.es["weight"] = w[nz].tolist()
    return G, names


# ---------- gene-expression matched read-out (same VAE machinery as embedding_recon) ----------
def _svd_emb(Xsub, k):
    _, _, Vt = np.linalg.svd(Xsub, full_matrices=False)
    kk = min(k, Vt.shape[0])
    e = np.mean(Xsub @ Vt[:kk, :].T, axis=0)
    if e.size < k:
        e = np.concatenate([e, np.zeros(k - e.size)])
    return e[:k]


def xgate_path_emb(G, graph_genes, K):
    """Fixed-dim xGATE pathway embedding for the scatter (first K of generate_embedding, padded)."""
    gg = set(graph_genes)
    idx = [v.index for v in G.vs if v["name"] in gg]
    if len(idx) < 2:
        return np.full(K, np.nan)
    e = np.asarray(generate_embedding(G.subgraph(idx)), dtype=float)
    if e.size < K:
        e = np.concatenate([e, np.zeros(K - e.size)])
    return e[:K]


def geneexpr_recon(Xdonor, path_cols, B, k, rng):
    """Matched gene-expression read-out z: SVD expression embedding of the pathway genes vs
    B random size-matched gene sets, scored by the same VAE-reconstruction-error procedure as
    xGATE's embedding_recon (train VAE on one null draw, score pathway against a second draw)."""
    n = len(path_cols)
    if n < 2:
        return np.nan
    P = Xdonor.shape[1]
    emb_p = _svd_emb(Xdonor[:, path_cols], k)
    neg = [_svd_emb(Xdonor[:, rng.sample(range(P), n)], k) for _ in range(B)]
    rnd = [_svd_emb(Xdonor[:, rng.sample(range(P), n)], k) for _ in range(B)]
    df = pd.DataFrame(neg + rnd + [emb_p]).apply(normalize, axis=0)
    df = df.replace([np.inf, -np.inf], np.nan).fillna(0.0).values
    tr, sc, pe = df[:B], df[B:2 * B], df[-1:]
    dim = tr.shape[1]
    m = VariationalAutoencoder(dim, 16); opt = torch.optim.Adam(m.parameters(), lr=1e-3)
    X = torch.Tensor(tr)
    for _ in range(1000):
        m.zero_grad(); rx, mu, lv = m(X); loss = vae_loss_function(rx, X, mu, lv)
        if not torch.isfinite(loss):
            break
        loss.backward(); opt.step()
    re_p = np.array(calculate_reconstruction_error(pe.tolist(), m, vae_loss_function))
    re_n = np.array(calculate_reconstruction_error(sc.tolist(), m, vae_loss_function))
    return float((re_p[0] - re_n.mean()) / re_n.std()) if re_n.std() > 0 else np.nan


# ---------- per-donor driver ----------
def load_donor(dataset, donor):
    """Return (graph G, graph node-name set, Xexpr [cells x genes], expr gene->col map)."""
    import scanpy as sc, anndata as ad
    if dataset == "pancreas":
        G, names = load_graph(f"{PANC_ADJ}/adj_{donor}.csv")
        a = sc.read(PANC_H5)
        fd = a[(a.obs["cell_type"] == "type B pancreatic cell") &
               (a.obs["disease_state"] == "Control") & (a.obs["donor_id"] == donor)].copy()
        X = fd.X.toarray() if not isinstance(fd.X, np.ndarray) else fd.X
        g2i = {str(g).split(".")[0]: i for i, g in enumerate(fd.var_names)}
    else:
        G, names = load_graph(f"{TS_ADJ}/adj_{donor}.csv")
        a = ad.read_h5ad(TS_H5, backed="r")
        sel = (a.obs.cell_type == "fibroblast").values & (a.obs.donor_id == donor).values
        ii = np.where(sel)[0]
        if len(ii) > 1500:
            ii = np.sort(np.random.RandomState(42).choice(ii, 1500, replace=False))
        sub = a[ii].to_memory(); X = sub.X
        X = X.toarray() if not isinstance(X, np.ndarray) else X
        g2i = {str(g).split(".")[0]: i for i, g in enumerate(sub.var_names)}
    return G, set(names), X, g2i


def run(dataset, donor, B, seed, timing, emb_only=False):
    rng = random.Random(seed); np.random.seed(seed); torch.manual_seed(seed)
    bench = pd.read_csv(BENCH[dataset])
    bench["pathway"] = bench["pathway"].astype(str)
    print(f"[{dataset}/{donor}] benchmark pathways: {len(bench)} "
          f"({(bench.truth=='positive').sum()} pos / {(bench.truth=='negative').sum()} neg)", flush=True)

    # gene members per pathway in BOTH graph-space and expression-space
    hsa = {p: kegg_hsa_genes(p) for p in bench.pathway}
    all_hsa = sorted({g for v in hsa.values() for g in v})
    ens_map = hsa_to_ensembl(all_hsa)                      # hsa -> [ensembl]
    def to_ensembl(genes):
        out = []
        for g in genes:
            out.extend(ens_map.get(g, []))
        return sorted(set(out))

    G, node_names, X, g2i = load_donor(dataset, donor)
    print(f"[{dataset}/{donor}] graph |V|={G.vcount()} |E|={G.ecount()}; cells={X.shape[0]} genes={X.shape[1]}",
          flush=True)

    rows = bench.copy()
    if timing:
        rows = rows.iloc[[0, len(rows) - 1]].copy()          # one positive, one negative
        B = min(B, 20)
        print(f"[timing] scoring {len(rows)} pathways at B={B}", flush=True)

    out = []
    emb_x, emb_g, emb_state, emb_path = [], [], [], []
    for _, r in rows.iterrows():
        p = r["pathway"]; state = r["truth"]
        # graph-space genes
        if dataset == "pancreas":
            graph_genes = [g for g in to_ensembl(hsa[p]) if g in node_names]
            expr_genes = [g for g in to_ensembl(hsa[p]) if g in g2i]
        else:
            graph_genes = [g for g in hsa[p] if g in node_names]     # graph is hsa:entrez
            expr_genes = [g for g in to_ensembl(hsa[p]) if g in g2i]
        cols = [g2i[g] for g in expr_genes]
        n_graph, n_expr = len(graph_genes), len(cols)
        t0 = time.time()
        if emb_only:
            zx = zg = np.nan
            t1 = t2 = t0
        else:
            if n_graph >= 2:
                _p, zx = embedding_recon(G, None, graph_genes, 200, 200, B)
            else:
                zx = np.nan
            t1 = time.time()
            k = max(2, min(20, n_expr))
            zg = geneexpr_recon(X, cols, B, k, rng) if n_expr >= 2 else np.nan
            t2 = time.time()
        # fixed-dim embeddings for the scatter (same points as the read-out, labeled by state)
        emb_x.append(xgate_path_emb(G, graph_genes, EMB_K))
        emb_g.append(_svd_emb(X[:, cols], EMB_K) if n_expr >= 2 else np.full(EMB_K, np.nan))
        emb_state.append(state); emb_path.append(p)
        out.append(dict(dataset=dataset, donor=donor, pathway=p, state=state,
                        z_xgate=zx, z_ge=zg, n_graph=n_graph, n_expr=n_expr))
        print(f"  {p[:34]:34s} state={state:8s} n_g={n_graph:3d} n_e={n_expr:3d} "
              f"zX={zx:8.3f} zG={zg:8.3f}  [xg {t1-t0:5.1f}s ge {t2-t1:4.1f}s]", flush=True)

    df = pd.DataFrame(out)
    if not timing:
        np.savez(f"{OUT}/results/supp_figS3_5_readout_{dataset}_{donor}_emb.npz",
                 xgate=np.asarray(emb_x, float), geneexpr=np.asarray(emb_g, float),
                 state=np.asarray(emb_state), pathway=np.asarray(emb_path),
                 donor=np.asarray([donor] * len(emb_state)))
        if not emb_only:                     # emb_only keeps the existing z CSV untouched
            outp = f"{OUT}/results/supp_figS3_5_readout_{dataset}_{donor}.csv"
            df.to_csv(outp, index=False)
            print(f"[done] -> {outp}", flush=True)
        else:
            print(f"[emb_only done] -> supp_figS3_5_readout_{dataset}_{donor}_emb.npz", flush=True)
    else:
        print("[timing done] extrapolate: per-pathway xGATE seconds x n_pathways x n_donors", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=["pancreas", "ts"], required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--B", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--timing", action="store_true")
    ap.add_argument("--emb_only", action="store_true")   # export scatter embeddings, skip VAE scoring
    a = ap.parse_args()
    run(a.dataset, a.donor, a.B, a.seed, a.timing, a.emb_only)
