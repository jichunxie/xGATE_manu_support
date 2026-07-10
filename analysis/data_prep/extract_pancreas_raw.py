#!/usr/bin/env python
"""Extract RAW pancreas ctrl beta-cell counts (genes ensembl x cells) from the .raw layer
of pancreas_human.h5ad, for SCTransform input (its X is already normalized)."""
from xgate_paths import ROOT  # noqa: E402
import sys, numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
a=ad.read_h5ad(ROOT + "/data/pancreas/pancreas_human.h5ad")
m=(a.obs.get('cell_type')=='type B pancreatic cell')
if 'disease_state' in a.obs: m=m&(a.obs['disease_state']=='Control')
fd=a[m.values].copy()
raw=fd.raw.to_adata() if fd.raw is not None else fd
X=raw.X; X=X.toarray() if sp.issparse(X) else np.asarray(X)
genes=[str(g).split('.')[0] for g in raw.var_names]
df=pd.DataFrame(X.T,index=genes); df=df[~df.index.duplicated(keep='first')]
df=df.loc[df.sum(axis=1)>0,:]
df.to_csv(ROOT + "/data/pancreas/pancreas_raw_ensembl.csv")
print(f"pancreas RAW: {df.shape} genes x cells; integer? {np.allclose(df.values, np.round(df.values))}")
