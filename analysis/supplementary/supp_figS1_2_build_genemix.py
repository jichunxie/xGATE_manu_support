#!/usr/bin/env python
"""
R2.1 gene-mix alpha-sweep (partial pathway activity).

Reviewer 2.1: does xGATE represent PARTIAL pathway activity as something
intermediate between fully-active and fully-inactive? We answer with a
controlled synthetic continuum on the canonical pancreas-ctrl co-expression
graph (the same graph that reproduces the published Insulin p=0.005).

Design (locked, plan T1):
  - Active pool   = genes of an ACTIVE pathway   (default: Insulin signaling pathway).
  - Inactive pool = genes of an INACTIVE pathway (default: Bacterial invasion of
    epithelial cells).
  - Restrict both pools to DETECTED nodes in the graph, size-match to a fixed N
    (default min(|active_det|,|inactive_det|), capped to --N if given).
  - For each alpha in {1.0,0.75,0.5,0.25,0.0}: synthetic set =
        round(alpha*N) active genes + (N - round(alpha*N)) inactive genes,
    drawn from the matched pools, so EVERY synthetic set has exactly N genes.
  - Score each synthetic set with the published xGATE engine (generate_embedding
    + VAE recon vs a size-matched UNIFORM random null, exactly as embedding_recon).
    Record recon error (pathway), null mean/std, p-value, z-score, and the
    embedding distance-to-random (centroid of the random cloud) -- the Fig2c
    "distance from the random graphs" quantity.
  - Repeat over several seeds (--seeds) for stability; report mean +/- SD.
  - Also persist, for ONE representative seed, the normalized embedding matrix
    (row 0 = synthetic-pathway embedding, rows 1.. = random-null embeddings) per
    alpha, so the figure script can MDS-project them exactly like Fig2c
    (fig2c_embedding_comparison.ipynb: MDS(n_components=3).fit_transform(data.T),
    row0=pathway triangle, rest=grey random cloud).

Outputs:
  results/supp_figS1_2_genemix_sweep.csv                 (one row per alpha x seed + agg)
  results/supp_figS1_2_genemix_embeddings_alpha<a>.csv   (MDS input, representative seed)

Reuses matched_null.py's NaN-safe recon machinery so results are apples-to-apples
with the published engine.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, random, os, json
import numpy as np
import pandas as pd

import torch
from xGATE.utilities import (
    create_network_from_adj_matrix, get_categorized_pathways,
    gather_pathways_between, get_genes_in_pathway,
)
from xGATE.utilities.embeddings import generate_embedding
from xGATE.utilities.pathway_analysis import normalize
from xGATE.utilities.vae_model import (
    VariationalAutoencoder, vae_loss_function, calculate_reconstruction_error,
)

def decode(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

def resolve_detected(name, cats, node_set):
    """KEGG pathway name -> detected hsa:<entrez> genes present in the graph."""
    ids = gather_pathways_between(name, name, cats)
    genes = get_genes_in_pathway(ids) if ids else []
    return [g for g in genes if g in node_set]

def recon_and_embedding(G, found_genes, B, rng, save_embeddings=False):
    """Mirror published embedding_recon (uniform null), but return p, z, recon
    error, AND (optionally) the normalized embedding matrix for MDS, plus the
    embedding distance from the random-null centroid.

    NaN-safe (matched_null.py pattern): sanitize non-finite embeddings, clip
    grads, retry VAE seeds. Returns dict or None-flagged on failure."""
    pathway_sub = G.subgraph([v.index for v in G.vs if v["name"] in set(found_genes)])
    k = pathway_sub.vcount()
    all_names = [v["name"] for v in G.vs]

    def draw_subgraphs(n):
        subs = []
        for _ in range(n):
            verts = set(rng.sample(all_names, k))
            subs.append(G.subgraph([v.index for v in G.vs if v["name"] in verts]))
        return subs

    neg = draw_subgraphs(B)        # train VAE (published: a 2nd uniform draw)
    rnd = draw_subgraphs(B)        # null
    emb_neg = [generate_embedding(s) for s in neg]
    emb_rnd = [generate_embedding(s) for s in rnd]
    emb_pw  = generate_embedding(pathway_sub)

    df_all = pd.concat([pd.DataFrame(emb_neg), pd.DataFrame(emb_rnd),
                        pd.DataFrame([emb_pw])], axis=0).reset_index(drop=True)
    df_all = df_all.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    dfn = df_all.apply(normalize, axis=0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    e_neg = dfn.iloc[:B].values
    e_rnd = dfn.iloc[B:2*B].values
    e_pw  = dfn.iloc[-1].values

    # distance-to-random: euclidean from pathway embedding to centroid of the
    # random null cloud, in normalized feature space (Fig2c geometry, pre-MDS).
    rnd_centroid = e_rnd.mean(axis=0)
    dist_to_random = float(np.linalg.norm(e_pw - rnd_centroid))

    dim = len(e_pw)
    re_pw = re_rnd = None
    for attempt in range(4):
        torch.manual_seed(1234 + attempt)
        model = VariationalAutoencoder(dim, 16)
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        X = torch.Tensor(e_neg.tolist())
        for _ in range(1000):
            model.zero_grad()
            rx, mu, lv = model(X)
            loss = vae_loss_function(rx, X, mu, lv)
            if not torch.isfinite(loss):
                break
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
        rp = np.array(calculate_reconstruction_error([e_pw.tolist()], model, vae_loss_function))
        rr = np.array(calculate_reconstruction_error(e_rnd.tolist(), model, vae_loss_function))
        if np.all(np.isfinite(rp)) and np.all(np.isfinite(rr)) and rr.std() > 0:
            re_pw, re_rnd = rp, rr
            break
    if re_pw is None:
        return dict(k=k, EP=pathway_sub.ecount(), p=np.nan, z=np.nan,
                    re_pathway=np.nan, re_null_mean=np.nan, re_null_std=np.nan,
                    dist_to_random=dist_to_random, emb=None)
    p = float(np.mean(np.append(re_rnd, re_pw) >= re_pw))
    z = float((re_pw[0] - re_rnd.mean()) / re_rnd.std())
    out = dict(k=k, EP=pathway_sub.ecount(), p=p, z=z,
               re_pathway=float(re_pw[0]), re_null_mean=float(re_rnd.mean()),
               re_null_std=float(re_rnd.std()), dist_to_random=dist_to_random,
               emb=None)
    if save_embeddings:
        # row 0 = synthetic pathway embedding; rows 1.. = random-null embeddings
        emb_mat = np.vstack([e_pw, e_rnd])
        out["emb"] = emb_mat
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adj", required=True)
    ap.add_argument("--active", default="Insulin signaling pathway")
    ap.add_argument("--inactive", default="Bacterial invasion of epithelial cells")
    ap.add_argument("--alphas", nargs="+", type=float,
                    default=[1.0, 0.75, 0.5, 0.25, 0.0])
    ap.add_argument("--N", type=int, default=0,
                    help="size to match both pools to (0 = min of detected pools)")
    ap.add_argument("--seeds", nargs="+", type=int, default=[12, 24, 36])
    ap.add_argument("--emb-seed", type=int, default=12,
                    help="which seed's per-alpha embeddings to persist for MDS")
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--out", required=True)
    ap.add_argument("--emb-prefix", required=True,
                    help="path prefix for per-alpha embedding CSVs")
    args = ap.parse_args()

    t0 = time.time()
    random.seed(args.seeds[0]); np.random.seed(args.seeds[0]); torch.manual_seed(args.seeds[0])

    print(f"[load] {args.adj}", flush=True)
    adj = pd.read_csv(args.adj, index_col=0) if args.adj.endswith(".csv") else pd.read_hdf(args.adj)
    adj.index = [decode(i) for i in adj.index]; adj.columns = [decode(c) for c in adj.columns]
    G = create_network_from_adj_matrix(adj); del adj
    node_set = set(v["name"] for v in G.vs)
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    cats = get_categorized_pathways()
    act_pool = resolve_detected(args.active, cats, node_set)
    ina_pool = resolve_detected(args.inactive, cats, node_set)
    print(f"[pools] active '{args.active}' detected={len(act_pool)}; "
          f"inactive '{args.inactive}' detected={len(ina_pool)}", flush=True)
    if len(act_pool) < 15 or len(ina_pool) < 15:
        print("[FATAL] a pool has <15 detected genes; cannot size-match meaningfully.", flush=True)
        sys.exit(2)

    N = min(len(act_pool), len(ina_pool))
    if args.N and args.N > 0:
        N = min(N, args.N)
    print(f"[size-match] N={N} genes per synthetic set", flush=True)

    rows = []
    emb_store = {}   # alpha -> embedding matrix (representative seed only)
    for seed in args.seeds:
        rng = random.Random(seed)
        # fix the size-matched pools per seed (subsample N from each, then mix)
        act_matched = rng.sample(act_pool, N)
        ina_matched = rng.sample(ina_pool, N)
        for a in args.alphas:
            n_act = int(round(a * N))
            n_ina = N - n_act
            # draw without replacement from each matched pool
            genes = rng.sample(act_matched, n_act) + rng.sample(ina_matched, n_ina)
            genes = list(dict.fromkeys(genes))   # de-dup (pools are disjoint normally)
            save = (seed == args.emb_seed)
            r = recon_and_embedding(G, genes, args.B, rng, save_embeddings=save)
            rows.append(dict(
                alpha=a, seed=seed, N=N, n_active_genes=n_act, n_inactive_genes=n_ina,
                n_detected=len(genes), k_in_graph=r["k"], EP=r["EP"],
                p_value=r["p"], z_score=r["z"], recon_error=r["re_pathway"],
                null_mean=r["re_null_mean"], null_std=r["re_null_std"],
                dist_to_random=r["dist_to_random"]))
            print(f"  seed={seed} a={a:.2f}  n_act={n_act:3d} n_ina={n_ina:3d}  "
                  f"p={r['p'] if r['p']==r['p'] else float('nan'):.4f} "
                  f"z={r['z'] if r['z']==r['z'] else float('nan'):7.2f} "
                  f"dist={r['dist_to_random']:.4f} t={time.time()-t0:.0f}s", flush=True)
            if save and r["emb"] is not None:
                emb_store[a] = r["emb"]

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)

    # aggregate over seeds for the report
    agg = (df.groupby("alpha")
             .agg(p_mean=("p_value", "mean"), p_sd=("p_value", "std"),
                  z_mean=("z_score", "mean"), z_sd=("z_score", "std"),
                  dist_mean=("dist_to_random", "mean"), dist_sd=("dist_to_random", "std"),
                  re_mean=("recon_error", "mean"))
             .reset_index())
    print("\n[aggregate over seeds]", flush=True)
    print(agg.to_string(index=False), flush=True)

    # persist representative-seed embeddings per alpha for MDS (Fig2c style)
    os.makedirs(os.path.dirname(args.emb_prefix), exist_ok=True)
    for a, mat in emb_store.items():
        tag = f"{a:.2f}".replace(".", "p")
        path = f"{args.emb_prefix}_alpha{tag}.csv"
        edf = pd.DataFrame(mat)
        edf.insert(0, "row_type", ["pathway"] + ["random"] * (mat.shape[0] - 1))
        edf.to_csv(path, index=False)
        print(f"  [emb] alpha={a:.2f} -> {path} ({mat.shape[0]}x{mat.shape[1]})", flush=True)

    # monotonicity check (verification gate)
    if agg["dist_mean"].notna().all():
        d = agg.sort_values("alpha")["dist_mean"].values
        # expect distance to random to be MONOTONE across alpha (active end far,
        # inactive end near, partial in between). Report Spearman with alpha.
        from scipy.stats import spearmanr
        rho, pv = spearmanr(agg["alpha"], agg["dist_mean"])
        print(f"\n[monotonicity] Spearman(alpha, dist_to_random) rho={rho:.3f} p={pv:.3g}", flush=True)
    print(f"\n[done] t={time.time()-t0:.0f}s -> {args.out}", flush=True)

if __name__ == "__main__":
    main()
