#!/usr/bin/env python
"""Merge TS xGATE: 10 positives (valid, scored on the SCT graph in the old run) +
10 re-scored negatives (correct JSON/competing set) -> ts_fibroblast_xgate_sct.csv."""
from xgate_paths import ROOT  # noqa: E402
import pandas as pd
RES=ROOT + "/results"
old=pd.read_csv(f"{RES}/ts_fibroblast_xgate_sct_OLDposneg.csv")
pos=old[old.truth=="positive"].copy()                 # 10 positives, valid on SCT graph
neg=pd.read_csv(f"{RES}/ts_neg_sct.csv")              # 10 correct negatives
assert len(pos)==10, f"expected 10 positives, got {len(pos)}"
assert len(neg)==10, f"expected 10 negatives, got {len(neg)}"
out=pd.concat([pos,neg],ignore_index=True)
out.to_csv(f"{RES}/ts_fibroblast_xgate_sct.csv",index=False)
print(f"merged -> ts_fibroblast_xgate_sct.csv ({len(pos)} pos + {len(neg)} neg)")
print(out[["pathway","truth","p_value"]].to_string(index=False))
