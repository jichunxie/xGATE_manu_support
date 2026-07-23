#!/usr/bin/env python
"""Extract one pancreas donor's control beta-cell RAW counts (genes ensembl x cells) from
pancreas_human.h5ad .raw, for SCTransform (per-donor batch baseline, R1.3 1st dataset)."""
from xgate_paths import ROOT  # noqa: E402
import sys, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
donor=sys.argv[1]; out=sys.argv[2]
a=ad.read_h5ad(ROOT + "/data/pancreas/pancreas_human.h5ad")
m=(a.obs['cell_type']=='type B pancreatic cell')&(a.obs['disease_state']=='Control')&(a.obs['donor_id']==donor)
fd=a[m.values].copy()
raw=fd.raw.to_adata() if fd.raw is not None else fd
X=raw.X; X=X.toarray() if sp.issparse(X) else np.asarray(X)
genes=[str(g).split('.')[0] for g in raw.var_names]
df=pd.DataFrame(X.T,index=genes); df=df[~df.index.duplicated(keep='first')]
df=df.loc[df.sum(axis=1)>0,:]
df.to_csv(out)
print(f"{donor}: {df.shape} genes x cells -> {out}", flush=True)
