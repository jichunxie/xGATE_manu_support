#!/usr/bin/env python
"""
R1.1 FUCCI benchmark (GSE146773): build a co-expression graph from FACS-sorted
proliferating U2OS FUCCI cells (SMART-seq2, deep). Ground truth: cell-cycle /
DNA-replication / repair pathways should be active (coordinated) in proliferating
cells; immune / unrelated pathways should be inactive.

Counts CSV is cells x genes (ensembl columns); we transpose to genes x cells.
"""
from xgate_paths import ROOT  # noqa: E402
import sys, time, os, resource, socket, datetime
import numpy as np, pandas as pd
from xGATE.utilities import (create_sifinet_object, quantile_thres2, cal_coexp,
                       create_network, filter_lowexp, convert_gene_ids)

D = ROOT + "/data/fucci_u2os"
OUT = ROOT + "/data/fucci_u2os/adj_matrix_fucci_u2os.csv"
SUMMARY = ROOT + "/results/computational_benchmark_summary_additions.csv"

def cpu_model():
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return ""

def main():
    t0 = time.time()
    counts = pd.read_csv(f"{D}/GSE146773_Counts.csv.gz", index_col=0)  # cells x genes
    print(f"[load] counts {counts.shape} (cells x genes) t={time.time()-t0:.0f}s", flush=True)
    df = counts.T   # genes x cells
    df.index = [str(g).split(".")[0] for g in df.index]
    df = df[~df.index.duplicated(keep="first")]
    # gene QC then cell QC
    df = df.loc[(df > 0).sum(axis=1) >= 0.05 * df.shape[1], :]
    df = df.loc[:, (df > 0).sum(axis=0) >= 0.05 * df.shape[0]]
    rm = df.mean(axis=1)
    df = df.loc[~((rm == 0) | (rm == 1)), :]
    print(f"[filter] genes x cells = {df.shape} t={time.time()-t0:.0f}s", flush=True)

    so = create_sifinet_object(df, rowfeature=True)
    so = quantile_thres2(so)
    so = cal_coexp(so, X=so.data_thres['dt'], X_full=so.data_thres['dt'])
    so = create_network(so, alpha=0.05, manual=False, least_edge_prop=0.01)
    so = filter_lowexp(so, t1=10, t2=0.9, t3=0.9)
    print(f"[graph] thres={so.thres:.4g} t={time.time()-t0:.0f}s", flush=True)

    adj = pd.DataFrame(np.where(np.abs(so.coexp - so.est_ms['mean']) > so.thres,
                                np.abs(so.coexp), 0.0))
    adj.index = df.index; adj.columns = df.index
    nV = adj.shape[0]; nE = int((adj.values != 0).sum() / 2); n_cells = df.shape[1]
    wall = time.time() - t0
    print(f"[adj] |V|={nV} |E|={nE} t={wall:.0f}s", flush=True)
    adj = convert_gene_ids(adj, "ensembl")
    adj.to_csv(OUT)
    print(f"[done] saved {adj.shape} -> {OUT} t={time.time()-t0:.0f}s", flush=True)

    # full timing + hardware row, same columns as computational_benchmark_summary.csv
    ru = resource.getrusage(resource.RUSAGE_SELF)
    cols = ["file", "dataset", "n_genes", "n_vertices", "n_edges", "file_size_mb", "wall_time_s",
            "cpu_user_s", "cpu_sys_s", "cpu_total_s", "peak_rss_mb", "delta_rss_mb",
            "n_cpus_available", "cpu_model", "hostname", "timestamp", "master_file", "n_cells",
            "master_file_size_mb", "cluster", "dataset_label"]
    row = {c: "" for c in cols}
    row.update(dict(file=os.path.basename(OUT), dataset="fucci_u2os", n_genes=nV, n_vertices=nV,
        n_edges=nE, file_size_mb=round(os.path.getsize(OUT) / 1e6, 2), wall_time_s=round(wall, 1),
        cpu_user_s=round(ru.ru_utime, 1), cpu_sys_s=round(ru.ru_stime, 1),
        cpu_total_s=round(ru.ru_utime + ru.ru_stime, 1), peak_rss_mb=round(ru.ru_maxrss / 1024.0, 1),
        n_cpus_available=os.cpu_count(), cpu_model=cpu_model(), hostname=socket.gethostname(),
        timestamp=datetime.datetime.now().isoformat(timespec="seconds"), n_cells=n_cells,
        cluster="pooled", dataset_label="FUCCI U2OS"))
    pd.DataFrame([row])[cols].to_csv(SUMMARY, index=False)
    print(f"[timing] hardware row -> {SUMMARY}", flush=True)

if __name__ == "__main__":
    main()
