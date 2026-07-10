#!/usr/bin/env python
"""Build a per-donor xGATE co-expression graph for control pancreatic beta cells,
for the corrected donor-batch analysis (supp_figS3_5_build_batch_manyset.py). Uses raw counts from
pancreas_human.h5ad .raw (SCTransform/log-norm is the gene-EXPRESSION baseline; the
graph is built from counts as xGATE normally does). Node names = ensembl (matching the
h5ad var_names), so the graph and the expression matrix share one gene ID space and the
same random gene sets can be applied to both.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, os, time
import numpy as np, pandas as pd, scipy.sparse as sp
from xGATE.utilities import (create_sifinet_object, quantile_thres2, cal_coexp,
                       create_network, filter_lowexp)

H5 = ROOT + "/data/batch_effect/pancreas_human.h5ad"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True)
    ap.add_argument("--out", required=True)
    # manual=True floors the edge threshold at the least_edge_prop quantile so small-donor
    # graphs stay non-degenerate (FDR-only threshold collapses to ~0 edges for HPAP022).
    ap.add_argument("--manual", action="store_true")
    ap.add_argument("--least_edge_prop", type=float, default=0.01)
    args = ap.parse_args()
    t0 = time.time()
    import scanpy as sc
    a = sc.read(H5)
    m = ((a.obs["cell_type"] == "type B pancreatic cell") &
         (a.obs["disease_state"] == "Control") & (a.obs["donor_id"] == args.donor)).values
    sub = a[m].copy()
    # raw counts (h5ad .raw), ensembl var_names to match the gene-expression matrix
    if sub.raw is not None:
        X = sub.raw.X; genes = [str(g).split(".")[0] for g in sub.raw.var_names]
    else:
        X = sub.X; genes = [str(g).split(".")[0] for g in sub.var_names]
    X = sp.csr_matrix(X)
    print(f"[subset] {args.donor}: {X.shape[0]} control beta cells x {X.shape[1]} genes", flush=True)
    df = pd.DataFrame(X.T.toarray(), index=genes)      # genes x cells
    df = df[~df.index.duplicated(keep="first")]
    df = df.loc[(df > 0).sum(axis=1) >= 0.05 * df.shape[1], :]
    df = df.loc[:, (df > 0).sum(axis=0) >= 0.05 * df.shape[0]]
    rm = df.mean(axis=1); df = df.loc[~((rm == 0) | (rm == 1)), :]
    print(f"[filter] genes x cells = {df.shape} t={time.time()-t0:.0f}s", flush=True)

    so = create_sifinet_object(df, rowfeature=True)
    so = quantile_thres2(so)
    so = cal_coexp(so, X=so.data_thres["dt"], X_full=so.data_thres["dt"])
    so = create_network(so, alpha=0.05, manual=args.manual, least_edge_prop=args.least_edge_prop)
    so = filter_lowexp(so, t1=10, t2=0.9, t3=0.9)
    adj = pd.DataFrame(np.where(np.abs(so.coexp - so.est_ms["mean"]) > so.thres, np.abs(so.coexp), 0.0))
    adj.index = df.index; adj.columns = df.index
    nV = adj.shape[0]; nE = int((adj.values != 0).sum() / 2)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adj.to_csv(args.out)
    print(f"[done] {args.donor} |V|={nV} |E|={nE} t={time.time()-t0:.0f}s -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
