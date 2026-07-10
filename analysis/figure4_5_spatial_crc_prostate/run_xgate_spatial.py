#!/usr/bin/env python
"""
Spatial xGATE runner: score KEGG pathway activity per spatial subdomain
(Fig 4 CRC Visium HD; Fig 5 prostate Xenium Prime).

The R spatial pipeline (analysis/figure4_5_spatial_crc_prostate/analysis/{crc,prostate})
identifies domains with IRIS and splits each into K-means subclusters, giving one
subdomain per (domain, cluster). This script produces the xGATE pathway-activity CSV
that the plotting scripts (05/06_preprocess_pathway_data.R) consume but which was
previously generated off-repo. For each subdomain it:

  1. builds an xGATE gene co-expression graph from the subdomain's SCT-normalized
     expression matrix, using the same SifiNet quantile-association construction xGATE
     uses elsewhere (mirrors analysis/data_prep/panc_donor_graph.py), then
  2. scores each KEGG pathway with the published xGATE engine (embedding_recon),
     yielding an empirical p-value per pathway.

INPUT  --expr-dir DIR
  One CSV per subdomain, genes x cells: first column = gene symbols (index), header =
  cell IDs, values = the subdomain's SCT-transformed counts. Filenames encode the
  subdomain label as D<domain>_<cluster> (e.g. D1_2.csv = domain 1, cluster 2), matching
  the `subdomain` labels the R pipeline writes (see extract_domain_cluster in R/utils.R).
  Export one matrix per subdomain from the subdomain-annotated spatial object the R
  pipeline produces (e.g. sp_input_with_subdomains.rds): subset cells to each subdomain
  and write its genes x cells SCT matrix.

  --pathways / --pathways-file
  KEGG pathway *names* to score (same names used in Fig 3 / the pathway master list).
  Pass inline with --pathways "A" "B" ... or one-per-line via --pathways-file.

OUTPUT --out FILE
  Long CSV with columns: pathway, domain, cluster, subdomain, p-value, z_score,
  n_genes, n_detected, Type. The four columns pathway/domain/cluster/p-value are the
  schema required by 05/06_preprocess_pathway_data.R. Point the pipeline's PATHWAY_FILE
  (CRC: colon_pathway_analysis_results_xgate_subdomain.csv;
   prostate: prostate_pathway_analysis_results.csv) at this output.

Usage:
  run_xgate_spatial.py --expr-dir DIR --pathways-file pathways.txt \
      --out results/colon_pathway_analysis_results_xgate_subdomain.csv \
      [--convert symbol] [--B 200] [--seed 12] [--manual]

Node IDs default to gene symbols (10x spatial features); --convert symbol maps them to
the 'hsa:<entrez>' space xGATE scores in (default on). Deterministic seeds.
"""
import sys, os, re, time, glob, random, argparse
import numpy as np, pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shared"))
from xgate_paths import ROOT  # noqa: E402

import torch  # noqa: E402
from xGATE.utilities import (  # noqa: E402
    create_sifinet_object, quantile_thres2, cal_coexp, create_network, filter_lowexp,
    create_network_from_adj_matrix, convert_gene_ids, get_categorized_pathways,
    gather_pathways_between, get_genes_in_pathway, embedding_recon,
)

# Label like "D1_2" / "domain1_cluster2" / "1_2": first int = domain, second = cluster.
LABEL_RE = re.compile(r"(\d+)\D+(\d+)")


def parse_subdomain(fname):
    stem = os.path.splitext(os.path.basename(fname))[0]
    m = LABEL_RE.search(stem)
    if not m:
        raise ValueError(f"cannot parse D<domain>_<cluster> from filename '{fname}'")
    return int(m.group(1)), int(m.group(2)), stem


def load_expr(path):
    """Load a genes x cells expression matrix (first column = gene symbols)."""
    df = pd.read_csv(path, index_col=0)
    df.index = [str(g).split(".")[0] for g in df.index]     # strip version suffixes
    df = df[~df.index.duplicated(keep="first")]
    return df


def build_adj(expr, manual, least_edge_prop):
    """Build the xGATE co-expression adjacency for one subdomain's expression matrix.

    Mirrors analysis/data_prep/panc_donor_graph.py: the same 5% cell / 5% gene detection
    filter and SifiNet quantile-association network construction xGATE uses elsewhere.
    """
    df = expr.loc[(expr > 0).sum(axis=1) >= 0.05 * expr.shape[1], :]
    df = df.loc[:, (df > 0).sum(axis=0) >= 0.05 * df.shape[0]]
    rm = df.mean(axis=1)
    df = df.loc[~((rm == 0) | (rm == 1)), :]
    if df.shape[0] < 2 or df.shape[1] < 2:
        return None
    so = create_sifinet_object(df, rowfeature=True)
    so = quantile_thres2(so)
    so = cal_coexp(so, X=so.data_thres["dt"], X_full=so.data_thres["dt"])
    so = create_network(so, alpha=0.05, manual=manual, least_edge_prop=least_edge_prop)
    so = filter_lowexp(so, t1=10, t2=0.9, t3=0.9)
    adj = pd.DataFrame(np.where(np.abs(so.coexp - so.est_ms["mean"]) > so.thres,
                                np.abs(so.coexp), 0.0))
    adj.index = df.index
    adj.columns = df.index
    return adj


def read_pathways(args):
    if args.pathways_file:
        with open(args.pathways_file) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    return args.pathways


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr-dir", required=True,
                    help="directory of per-subdomain genes x cells SCT expression CSVs")
    ap.add_argument("--pathways", nargs="*", default=None, help="KEGG pathway names")
    ap.add_argument("--pathways-file", default=None, help="one KEGG pathway name per line")
    ap.add_argument("--out", required=True)
    ap.add_argument("--convert", choices=["ensembl", "symbol"], default="symbol",
                    help="gene-id space of the expression matrices (default: symbol)")
    ap.add_argument("--B", type=int, default=200, help="null-distribution size")
    ap.add_argument("--seed", type=int, default=12)
    ap.add_argument("--manual", action="store_true",
                    help="floor the edge threshold for small subdomains (see panc_donor_graph)")
    ap.add_argument("--least_edge_prop", type=float, default=0.01)
    args = ap.parse_args()

    pathways = read_pathways(args)
    if not pathways:
        ap.error("supply pathways via --pathways or --pathways-file")

    os.environ["PYTHONHASHSEED"] = str(args.seed)
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)

    files = sorted(glob.glob(os.path.join(args.expr_dir, "*.csv")))
    if not files:
        ap.error(f"no *.csv expression matrices found in {args.expr_dir}")
    cats = get_categorized_pathways()
    t0 = time.time()
    rows = []

    for f in files:
        domain, cluster, label = parse_subdomain(f)
        expr = load_expr(f)
        print(f"[{label}] {expr.shape[0]} genes x {expr.shape[1]} cells", flush=True)
        adj = build_adj(expr, args.manual, args.least_edge_prop)
        if adj is None or adj.shape[0] < 2:
            print(f"  [WARN] {label}: too few genes/cells after filtering; skipped", flush=True)
            continue
        if args.convert:
            adj = convert_gene_ids(adj, args.convert)
        G = create_network_from_adj_matrix(adj)
        node_names = set(v["name"] for v in G.vs)
        print(f"  |V|={G.vcount()} |E|={G.ecount()} t={time.time()-t0:.0f}s", flush=True)

        for name in pathways:
            ids = gather_pathways_between(name, name, cats)
            if not ids:
                print(f"  [WARN] no KEGG id for '{name}' (name mismatch)", flush=True)
                rows.append({"pathway": name, "domain": domain, "cluster": cluster,
                             "subdomain": label, "p-value": np.nan, "z_score": np.nan,
                             "n_genes": 0, "n_detected": 0, "Type": "xGATE"})
                continue
            genes = get_genes_in_pathway(ids)
            detected = [g for g in genes if g in node_names]
            p, z = embedding_recon(G, cats, genes, 200, 200, args.B)
            rows.append({"pathway": name, "domain": domain, "cluster": cluster,
                         "subdomain": label, "p-value": p, "z_score": z,
                         "n_genes": len(genes), "n_detected": len(detected),
                         "Type": "xGATE"})
            print(f"    {name:42s} p={p:.4f} z={z:8.2f} det={len(detected)}/{len(genes)}",
                  flush=True)

    df = pd.DataFrame(rows, columns=["pathway", "domain", "cluster", "subdomain",
                                     "p-value", "z_score", "n_genes", "n_detected", "Type"])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n[done] {len(df)} rows across {df['subdomain'].nunique()} subdomains "
          f"t={time.time()-t0:.0f}s -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
