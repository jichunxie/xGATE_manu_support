#!/usr/bin/env python
"""
R1.4 (T4): full edge-density & degree distributions -- naive-random does NOT match
the pathway, the global-degree-matched null DOES, yet active pathways still separate.

Reviewer 1.4: is uniform random node-sampling a fair null when pathway genes have
particular degree/density structure? We make the comparison explicit by dumping the
FULL distributions (not just medians) for representative active and inactive pathways:

  For each chosen pathway:
    - pathway induced subgraph: its edge density + its degree sequence.
    - NAIVE random null (B draws): uniform random gene sets of the same size ->
      induced edge density distribution + pooled degree distribution.
    - DEGREE-MATCHED null (B draws): random gene sets matched on each gene's GLOBAL
      degree (quantile bins) -> induced edge density distribution + pooled degree dist.

This shows: (a) naive-random density/degree distributions are shifted away from the
pathway (mismatch -> motivates the matched null); (b) the degree-matched null density
distribution overlaps the pathway's hubbiness (it genuinely matches the nuisance);
(c) yet active pathways still separate in xGATE (carried by supp_figS23_pancreas_matched_null.csv
p_matched, plotted in the figure).

Reuses matched_null.py's degree-matched sampler. Heavy (reads 414MB graph) -> sbatch.

Output: results/supp_figS22_distributions.npz  (per-pathway arrays for the figure script)
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, random, os
import numpy as np
import pandas as pd

sys.path.insert(0, ROOT + "/analysis/shared")
from xGATE.utilities import (
    create_network_from_adj_matrix, get_categorized_pathways,
    gather_pathways_between, get_genes_in_pathway,
)
from matched_null import make_degree_matched_sampler, uniform_sampler

def decode(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)

def induced(G, names):
    return G.subgraph([v.index for v in G.vs if v["name"] in set(names)])

def density(g):
    n = g.vcount()
    return (2.0 * g.ecount() / (n * (n - 1))) if n > 1 else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adj", required=True)
    ap.add_argument("--pathways", nargs="+", required=True)
    ap.add_argument("--truth", nargs="+", required=True)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--seed", type=int, default=12)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    random.seed(args.seed); np.random.seed(args.seed)
    t0 = time.time()
    adj = pd.read_csv(args.adj, index_col=0) if args.adj.endswith(".csv") else pd.read_hdf(args.adj)
    adj.index = [decode(i) for i in adj.index]; adj.columns = [decode(c) for c in adj.columns]
    G = create_network_from_adj_matrix(adj); del adj
    names = [v["name"] for v in G.vs]; node_set = set(names)
    gdeg = dict(zip(names, G.degree()))
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    matched = make_degree_matched_sampler(names, gdeg, n_bins=args.nbins)
    cats = get_categorized_pathways()
    rng = random.Random(args.seed)

    store = {}
    for name, tr in zip(args.pathways, args.truth):
        ids = gather_pathways_between(name, name, cats)
        genes = get_genes_in_pathway(ids) if ids else []
        found = [g for g in genes if g in node_set]
        if len(found) < 15:
            print(f"  [skip] {name}: {len(found)} detected (<15)", flush=True); continue
        psub = induced(G, found); k = psub.vcount()
        path_dens = density(psub)
        path_deg = np.array(psub.degree(), dtype=float)

        naive_dens, naive_deg = [], []
        match_dens, match_deg = [], []
        for _ in range(args.B):
            nv = set(uniform_sampler(names, k, rng))
            sg = induced(G, nv); naive_dens.append(density(sg)); naive_deg.extend(sg.degree())
            mv = set(matched(names, k, rng, target_genes=found))
            sm = induced(G, mv); match_dens.append(density(sm)); match_deg.extend(sm.degree())
        tag = name.replace(" ", "_").replace("/", "_")
        store[f"{tag}__truth"] = tr
        store[f"{tag}__path_dens"] = np.array([path_dens])
        store[f"{tag}__path_deg"] = path_deg
        store[f"{tag}__naive_dens"] = np.array(naive_dens)
        store[f"{tag}__naive_deg"] = np.array(naive_deg, dtype=float)
        store[f"{tag}__match_dens"] = np.array(match_dens)
        store[f"{tag}__match_deg"] = np.array(match_deg, dtype=float)
        print(f"  [{tr:8s}] {name:38s} k={k} pathDens={path_dens:.4f} "
              f"naiveDensMed={np.median(naive_dens):.4f} matchDensMed={np.median(match_dens):.4f} "
              f"t={time.time()-t0:.0f}s", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    # store pathway names list for the figure
    store["__pathways"] = np.array(args.pathways, dtype=object)
    store["__truth"] = np.array(args.truth, dtype=object)
    np.savez(args.out, **store)
    print(f"\n[done] t={time.time()-t0:.0f}s -> {args.out}", flush=True)

if __name__ == "__main__":
    main()
