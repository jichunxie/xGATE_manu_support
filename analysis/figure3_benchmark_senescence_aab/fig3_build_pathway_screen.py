#!/usr/bin/env python
"""
R1.2: full-catalog xGATE-only p-value screen on a canonical graph.
Run xGATE on ALL KEGG hsa pathways passing a detected-gene filter; save raw p / z.
Used downstream for p-value-distribution + BH/BY diagnostics. xGATE-only (no
competing methods) -- this is a calibration check, NOT a method comparison; unlabeled
pathways are never treated as false positives.

Incremental CSV append + resume so a SLURM timeout doesn't lose progress.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, os, random
import numpy as np, pandas as pd
import torch
from xGATE.utilities import (create_network_from_adj_matrix, get_categorized_pathways,
                       gather_pathways_between, get_genes_in_pathway, embedding_recon)

def decode(x): return x.decode() if isinstance(x,(bytes,bytearray)) else str(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adj", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-detected", type=int, default=15)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()
    os.environ["PYTHONHASHSEED"] = str(args.seed)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    t0 = time.time()

    adj = pd.read_csv(args.adj, index_col=0) if args.adj.endswith(".csv") else pd.read_hdf(args.adj)
    adj.index=[decode(i) for i in adj.index]; adj.columns=[decode(c) for c in adj.columns]
    G = create_network_from_adj_matrix(adj); del adj
    node_set = set(v["name"] for v in G.vs)
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    cats = get_categorized_pathways()
    # cats maps pathway-name -> [hsa id]; iterate all unique names
    names = sorted(cats.keys())
    print(f"[catalog] {len(names)} KEGG pathway names", flush=True)

    done = set()
    if os.path.exists(args.out):
        done = set(pd.read_csv(args.out)["pathway"]); print(f"[resume] {len(done)} done", flush=True)
    else:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        pd.DataFrame(columns=["pathway","category","n_genes","n_detected","p_value","z_score"]).to_csv(args.out, index=False)

    n_run = 0
    for name in names:
        if name in done: continue
        try:
            ids = gather_pathways_between(name, name, cats)
            genes = get_genes_in_pathway(ids) if ids else []
            det = [g for g in genes if g in node_set]
            if len(det) < args.min_detected:
                row = dict(pathway=name, category=name.split(" - ")[0], n_genes=len(genes),
                           n_detected=len(det), p_value=np.nan, z_score=np.nan)
            else:
                p, z = embedding_recon(G, cats, genes, 200, 200, args.B)
                row = dict(pathway=name, category=name.split(" - ")[0], n_genes=len(genes),
                           n_detected=len(det), p_value=p, z_score=z)
        except Exception as e:
            row = dict(pathway=name, category="", n_genes=-1, n_detected=-1,
                       p_value=np.nan, z_score=np.nan)
            print(f"  [err] {name}: {type(e).__name__}", flush=True)
        pd.DataFrame([row]).to_csv(args.out, mode="a", header=False, index=False)
        n_run += 1
        if n_run % 20 == 0:
            print(f"  ...{n_run} run, last={name[:30]} t={time.time()-t0:.0f}s", flush=True)
    print(f"[done] ran {n_run} pathways t={time.time()-t0:.0f}s -> {args.out}", flush=True)

if __name__ == "__main__":
    main()
