#!/usr/bin/env python
"""Extract one TS donor's fibroblast RAW counts (genes ensembl x cells) from Stromal.h5ad,
for SCTransform input (per-donor batch graphs, R1.3 2nd dataset). Light gene filter only;
SCTransform + graph build do the rest."""
import sys, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
TS="/path/to/group/Data/human_atlas/Tabula Sapiens/Stromal.h5ad"
donor=sys.argv[1]; out=sys.argv[2]
a=ad.read_h5ad(TS, backed="r")
m=(a.obs.cell_type=="fibroblast").values & (a.obs.donor_id==donor).values
idx=np.where(m)[0]
print(f"[{donor}] {len(idx)} fibroblast cells", flush=True)
sub=a[idx].to_memory()
if sub.raw is not None: X=sub.raw.X; genes=list(sub.raw.var_names)
else: X=sub.X; genes=list(sub.var_names)
X=sp.csr_matrix(X)
genes=[str(g).split(".")[0] for g in genes]
df=pd.DataFrame(X.T.toarray(), index=genes)        # genes x cells
df=df[~df.index.duplicated(keep="first")]
df=df.loc[df.sum(axis=1)>0, :]
df.to_csv(out)
print(f"[{donor}] raw {df.shape} genes x cells; integer? {np.allclose(df.values, np.round(df.values))} -> {out}", flush=True)
