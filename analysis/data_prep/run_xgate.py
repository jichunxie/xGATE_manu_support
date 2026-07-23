#!/usr/bin/env python
"""
General xGATE runner: load a co-expression adjacency matrix (csv or pandas-h5),
build the igraph network, and score a list of KEGG pathway *names* using the
published xGATE engine (embedding_recon). Bypasses the analyze_pathways import
bug by calling embedding_recon directly. Deterministic seeds.

Usage:
  run_xgate.py --adj PATH [--key adj_matrix] --pathways "A" "B" ... \
               --out results.csv [--already-hsa] [--B 200] [--seed 12]

Node IDs in the adj are expected to be 'hsa:<entrez>' (xGATE convention). If the
matrix carries ensembl/symbol IDs instead, pass --convert {ensembl,symbol}.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, random, os
import numpy as np
import pandas as pd

import torch
from xGATE.utilities import (
    create_network_from_adj_matrix, get_categorized_pathways,
    gather_pathways_between, get_genes_in_pathway, embedding_recon, convert_gene_ids,
)

def decode(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

def load_adj(path, key):
    if path.endswith(".csv"):
        adj = pd.read_csv(path, index_col=0)
    else:
        adj = pd.read_hdf(path, key=key)
    adj.index   = [decode(i) for i in adj.index]
    adj.columns = [decode(c) for c in adj.columns]
    return adj

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adj", required=True)
    ap.add_argument("--key", default="adj_matrix")
    ap.add_argument("--pathways", nargs="+", required=True)
    ap.add_argument("--truth", nargs="*", default=None,
                    help="optional same-length list of active/inactive labels")
    ap.add_argument("--out", required=True)
    ap.add_argument("--convert", choices=["ensembl", "symbol"], default=None)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()

    os.environ["PYTHONHASHSEED"] = str(args.seed)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)

    t0 = time.time()
    print(f"[load] {args.adj}", flush=True)
    adj = load_adj(args.adj, args.key)
    print(f"  shape={adj.shape} ids={list(adj.index[:3])} "
          f"nnz_frac={(adj.values!=0).mean():.4g} t={time.time()-t0:.0f}s", flush=True)
    if args.convert:
        adj = convert_gene_ids(adj, args.convert)
        print(f"  converted ({args.convert}); ids now {list(adj.index[:3])}", flush=True)

    G = create_network_from_adj_matrix(adj)
    del adj
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    cats = get_categorized_pathways()
    node_names = set(v["name"] for v in G.vs)
    truth = args.truth if args.truth else [""]*len(args.pathways)

    rows = []
    for name, tr in zip(args.pathways, truth):
        ids = gather_pathways_between(name, name, cats)
        if not ids:
            print(f"  [WARN] no KEGG id for '{name}' (name mismatch)", flush=True)
            rows.append({"pathway": name, "truth": tr, "p_value": np.nan,
                         "z_score": np.nan, "n_genes": 0, "n_detected": 0}); continue
        genes = get_genes_in_pathway(ids)
        detected = [g for g in genes if g in node_names]
        p, z = embedding_recon(G, cats, genes, 200, 200, args.B)
        rows.append({"pathway": name, "truth": tr, "p_value": p, "z_score": z,
                     "n_genes": len(genes), "n_detected": len(detected)})
        print(f"  [{tr:8s}] {name:42s} p={p:.4f} z={z:8.2f} "
              f"det={len(detected)}/{len(genes)} t={time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n[done] t={time.time()-t0:.0f}s -> {args.out}", flush=True)
    print(df.to_string(index=False), flush=True)

if __name__ == "__main__":
    main()
