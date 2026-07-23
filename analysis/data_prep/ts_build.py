#!/usr/bin/env python
"""
Build an xGATE co-expression graph from a Tabula Sapiens compartment h5ad for a given
cell type, either pooled (benchmark) or per-donor (batch effect). Records graph-build
timing in the same columns as computational_benchmark_summary.csv.

TS h5ad (cellxgene format): X is normalized; raw counts in .raw.X (preferred for xGATE).
Genes named by ensembl in var_names (var.feature_id) -> convert to hsa:entrez.
"""
from xgate_paths import ROOT  # noqa: E402
import argparse, sys, time, os, json, resource, socket, datetime
import numpy as np, pandas as pd, scipy.sparse as sp
from xGATE.utilities import (create_sifinet_object, quantile_thres2, cal_coexp,
                       create_network, filter_lowexp, convert_gene_ids)

TS = "/path/to/group/Data/human_atlas/Tabula Sapiens"

def cpu_model():
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compartment", required=True)
    ap.add_argument("--celltype", required=True)
    ap.add_argument("--donor", default=None)       # None -> pooled
    ap.add_argument("--subsample", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="TS")        # dataset_label for complexity csv
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()
    np.random.seed(args.seed)
    t0 = time.time()
    import anndata as ad
    a = ad.read_h5ad(f"{TS}/{args.compartment}.h5ad", backed="r")
    m = (a.obs.cell_type == args.celltype).values
    if args.donor:
        m &= (a.obs.donor_id == args.donor).values
    idx = np.where(m)[0]
    if args.subsample and len(idx) > args.subsample:
        idx = np.sort(np.random.choice(idx, args.subsample, replace=False))
    print(f"[subset] {args.celltype} donor={args.donor or 'POOLED'}: {len(idx)} cells", flush=True)
    sub = a[idx].to_memory()
    # choose raw counts
    if sub.raw is not None:
        X = sub.raw.X; genes = list(sub.raw.var_names)
    elif "counts" in (sub.layers.keys() if sub.layers else []):
        X = sub.layers["counts"]; genes = list(sub.var_names)
    else:
        X = sub.X; genes = list(sub.var_names)
    X = sp.csr_matrix(X)            # cells x genes
    genes = [str(g).split(".")[0] for g in genes]
    print(f"[counts] {X.shape} (cells x genes); tissue mix: "
          f"{sub.obs.tissue_in_publication.value_counts().head(6).to_dict()}", flush=True)

    df = pd.DataFrame(X.T.toarray(), index=genes)   # genes x cells
    df = df[~df.index.duplicated(keep="first")]
    df = df.loc[(df > 0).sum(axis=1) >= 0.05 * df.shape[1], :]
    df = df.loc[:, (df > 0).sum(axis=0) >= 0.05 * df.shape[0]]
    rm = df.mean(axis=1); df = df.loc[~((rm == 0) | (rm == 1)), :]
    print(f"[filter] genes x cells = {df.shape} t={time.time()-t0:.0f}s", flush=True)
    n_cells = df.shape[1]

    so = create_sifinet_object(df, rowfeature=True)
    so = quantile_thres2(so)
    so = cal_coexp(so, X=so.data_thres['dt'], X_full=so.data_thres['dt'])
    so = create_network(so, alpha=0.05, manual=False, least_edge_prop=0.01)
    so = filter_lowexp(so, t1=10, t2=0.9, t3=0.9)
    adj = pd.DataFrame(np.where(np.abs(so.coexp - so.est_ms['mean']) > so.thres, np.abs(so.coexp), 0.0))
    adj.index = df.index; adj.columns = df.index
    nV = adj.shape[0]; nE = int((adj.values != 0).sum() / 2)
    wall = time.time() - t0
    print(f"[adj] |V|={nV} |E|={nE} t={wall:.0f}s", flush=True)
    adj = convert_gene_ids(adj, "ensembl")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adj.to_csv(args.out)

    # timing row in summary-csv format
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
    cols=["file","dataset","n_genes","n_vertices","n_edges","file_size_mb","wall_time_s",
          "cpu_user_s","cpu_sys_s","cpu_total_s","peak_rss_mb","delta_rss_mb","n_cpus_available",
          "cpu_model","hostname","timestamp","master_file","n_cells","master_file_size_mb","cluster","dataset_label"]
    ru = resource.getrusage(resource.RUSAGE_SELF)
    row={c:"" for c in cols}
    row.update(dict(file=os.path.basename(args.out), dataset=args.celltype.replace(" ","_"),
        n_genes=nV, n_vertices=nV, n_edges=nE,
        file_size_mb=round(os.path.getsize(args.out)/1e6,2), wall_time_s=round(wall,1),
        cpu_user_s=round(ru.ru_utime,1), cpu_sys_s=round(ru.ru_stime,1),
        cpu_total_s=round(ru.ru_utime+ru.ru_stime,1), peak_rss_mb=round(peak_mb,1),
        n_cells=n_cells, hostname=socket.gethostname(),
        cpu_model=cpu_model(), n_cpus_available=os.cpu_count(),
        timestamp=datetime.datetime.now().isoformat(timespec="seconds"),
        cluster=(args.donor or "pooled"), dataset_label=args.label))
    tcsv = args.out.replace(".csv", "_timing.csv")
    pd.DataFrame([row])[cols].to_csv(tcsv, index=False)
    print(f"[done] saved {adj.shape} -> {args.out}; timing -> {tcsv}", flush=True)

if __name__ == "__main__":
    main()
