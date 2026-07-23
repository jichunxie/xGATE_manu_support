#!/usr/bin/env python
"""
R1.4: global-degree-matched null sensitivity.

Reviewer concern: uniform random node sampling may not be the right null; will the
random subgraphs have the pathway's edge density / degree distribution? co-expression
has complex covariance.

Response logic (marginal vs structural matching):
  - DEFAULT null  = size-matched uniform random gene sets (published xGATE).
  - MATCHED null  = random gene sets matched on each gene's GLOBAL degree in the full
    co-expression graph. This is a legitimate nuisance control: it removes "these genes
    are just globally hubby" without conditioning on the within-pathway wiring we test.
    (We deliberately do NOT match the *induced* edge density / degree sequence, which
    would condition on the activity signal itself.)

For each benchmark pathway we report default vs matched p/z, plus graph diagnostics
(|V_P|, |E_P|, induced edge density, median global degree of pathway genes), and the
null distribution of induced edge density (to show edge density alone does not separate
active from inactive, whereas xGATE does).

Reuses the published generate_embedding + VAE so results are apples-to-apples.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, random, os
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

# ---------- samplers ----------
def uniform_sampler(all_names, k, rng, **kw):
    return rng.sample(all_names, k)

def make_degree_matched_sampler(all_names, global_deg, n_bins=20):
    """Return sampler that draws k genes matched to a target degree profile.
    Bins are quantile bins of global degree; each target gene is replaced by a
    random gene from the same (or nearest non-empty) bin, without replacement."""
    deg = np.array([global_deg[n] for n in all_names], dtype=float)
    # quantile bin edges (unique to avoid empty bins on ties)
    qs = np.unique(np.quantile(deg, np.linspace(0, 1, n_bins + 1)))
    bin_idx = np.clip(np.digitize(deg, qs[1:-1]), 0, len(qs) - 2)
    name_arr = np.array(all_names, dtype=object)
    bins = {}
    for b in np.unique(bin_idx):
        bins[int(b)] = name_arr[bin_idx == b].tolist()
    nbins_eff = len(bins)
    name_to_bin = {n: int(b) for n, b in zip(all_names, bin_idx)}

    def sampler(_all_names, k, rng, target_genes=None):
        # target degree profile = bins of the pathway genes
        tgt_bins = [name_to_bin[g] for g in target_genes]
        chosen, used = [], set()
        for tb in tgt_bins:
            # search outward from target bin for an available gene
            for off in range(nbins_eff):
                for b in ({tb + off, tb - off} if off else {tb}):
                    pool = bins.get(b)
                    if not pool:
                        continue
                    cand = [g for g in pool if g not in used]
                    if cand:
                        pick = rng.choice(cand)
                        chosen.append(pick); used.add(pick)
                        break
                else:
                    continue
                break
        # pad if short (shouldn't happen)
        while len(chosen) < k:
            pick = rng.choice(_all_names)
            if pick not in used:
                chosen.append(pick); used.add(pick)
        return chosen
    return sampler

# ---------- core recon (mirrors published embedding_recon, swappable sampler) ----------
def recon(G, found_genes, sampler, B, rng, target_genes=None):
    pathway_sub = G.subgraph([v.index for v in G.vs if v["name"] in set(found_genes)])
    k = pathway_sub.vcount()
    all_names = [v["name"] for v in G.vs]

    def draw_subgraphs(n):
        subs = []
        for _ in range(n):
            verts = set(sampler(all_names, k, rng, target_genes=target_genes))
            subs.append(G.subgraph([v.index for v in G.vs if v["name"] in verts]))
        return subs

    neg = draw_subgraphs(B)        # train VAE
    rnd = draw_subgraphs(B)        # null
    emb_neg = [generate_embedding(s) for s in neg]
    emb_rnd = [generate_embedding(s) for s in rnd]
    emb_pw  = generate_embedding(pathway_sub)

    df_all = pd.concat([pd.DataFrame(emb_neg), pd.DataFrame(emb_rnd),
                        pd.DataFrame([emb_pw])], axis=0).reset_index(drop=True)
    # sanitize non-finite embedding entries (eccentricity on disconnected subgraphs etc.)
    df_all = df_all.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    dfn = df_all.apply(normalize, axis=0).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    e_neg = dfn.iloc[:B].values.tolist()
    e_rnd = dfn.iloc[B:2*B].values.tolist()
    e_pw  = dfn.iloc[-1].values.tolist()

    dim = len(e_pw)
    # train VAE with gradient clipping; retry a few seeds if RE comes back non-finite
    re_pw = re_rnd = None
    for attempt in range(4):
        torch.manual_seed(1234 + attempt)
        model = VariationalAutoencoder(dim, 16)
        opt = torch.optim.Adam(model.parameters(), lr=1e-3)
        X = torch.Tensor(e_neg)
        for _ in range(1000):
            model.zero_grad()
            rx, mu, lv = model(X)
            loss = vae_loss_function(rx, X, mu, lv)
            if not torch.isfinite(loss):
                break
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
        rp  = np.array(calculate_reconstruction_error([e_pw], model, vae_loss_function))
        rr  = np.array(calculate_reconstruction_error(e_rnd, model, vae_loss_function))
        if np.all(np.isfinite(rp)) and np.all(np.isfinite(rr)) and rr.std() > 0:
            re_pw, re_rnd = rp, rr
            break
    if re_pw is None:   # all attempts non-finite -> flag, do not emit spurious 0
        return dict(k=k, EP=pathway_sub.ecount(),
                    densP=(2.0*pathway_sub.ecount()/(k*(k-1)) if k > 1 else 0.0),
                    p=np.nan, z=np.nan, null_dens_med=np.nan, null_dens_p95=np.nan)
    p = float(np.mean(np.append(re_rnd, re_pw) >= re_pw))
    z = float((re_pw[0] - re_rnd.mean()) / re_rnd.std())

    # induced-density diagnostics
    dens = lambda g: (2.0 * g.ecount() / (g.vcount() * (g.vcount() - 1))) if g.vcount() > 1 else 0.0
    null_dens = np.array([dens(s) for s in rnd])
    return dict(k=k, EP=pathway_sub.ecount(), densP=dens(pathway_sub),
                p=p, z=z, null_dens_med=float(np.median(null_dens)),
                null_dens_p95=float(np.quantile(null_dens, 0.95)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adj", required=True)
    ap.add_argument("--pathways", nargs="+", required=True)
    ap.add_argument("--truth", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--nbins", type=int, default=20)
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()

    os.environ["PYTHONHASHSEED"] = str(args.seed)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    t0 = time.time()

    adj = pd.read_csv(args.adj, index_col=0) if args.adj.endswith(".csv") else pd.read_hdf(args.adj)
    adj.index = [decode(i) for i in adj.index]; adj.columns = [decode(c) for c in adj.columns]
    G = create_network_from_adj_matrix(adj); del adj
    print(f"[network] |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

    names = [v["name"] for v in G.vs]
    gdeg = dict(zip(names, G.degree()))
    node_set = set(names)
    matched = make_degree_matched_sampler(names, gdeg, n_bins=args.nbins)
    cats = get_categorized_pathways()

    rng_u = random.Random(args.seed)
    rng_m = random.Random(args.seed)
    rows = []
    for name, tr in zip(args.pathways, args.truth):
        ids = gather_pathways_between(name, name, cats)
        genes = get_genes_in_pathway(ids) if ids else []
        found = [g for g in genes if g in node_set]
        if len(found) < 5:
            print(f"  [skip] {name}: only {len(found)} genes detected", flush=True); continue
        med_gdeg = float(np.median([gdeg[g] for g in found]))
        d = recon(G, found, uniform_sampler, args.B, rng_u)
        m = recon(G, found, matched,        args.B, rng_m, target_genes=found)
        call_d = d["p"] < 0.05; call_m = m["p"] < 0.05
        rows.append(dict(pathway=name, truth=tr, n_detected=len(found),
            EP=d["EP"], densP=round(d["densP"], 5), med_global_deg=med_gdeg,
            p_default=round(d["p"],4), z_default=round(d["z"],2),
            p_matched=round(m["p"],4), z_matched=round(m["z"],2),
            default_null_densmed=round(d["null_dens_med"],5),
            matched_null_densmed=round(m["null_dens_med"],5),
            call_default=call_d, call_matched=call_m, call_stable=(call_d==call_m)))
        print(f"  [{tr:8s}] {name:38s} p_def={d['p']:.4f}(z{d['z']:6.1f}) "
              f"p_mat={m['p']:.4f}(z{m['z']:6.1f}) densP={d['densP']:.4f} "
              f"stable={call_d==call_m} t={time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n[done] t={time.time()-t0:.0f}s -> {args.out}", flush=True)
    print(df.to_string(index=False), flush=True)
    if len(df):
        print(f"\n[summary] calls stable under degree-matching: "
              f"{df.call_stable.sum()}/{len(df)}", flush=True)

if __name__ == "__main__":
    main()
