#!/usr/bin/env python
"""Cache the full human KEGG pathway panel as {pathway_name: [hsa:entrez, ...]}.

Used by supp_figS3_5_build_batch_manyset.py so the donor-separability comparison samples REAL pathways
(no random gene sets, no nulls). Fetch once; both datasets reuse the cache.
"""
from xgate_paths import ROOT  # noqa: E402
import json, sys, time
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser

OUT = ROOT + "/data/kegg_pathways_hsa.json"


def main():
    lines = REST.kegg_list("pathway", "hsa").read().splitlines()
    panel = {}
    for i, line in enumerate(lines):
        entry, desc = line.split("\t")
        pid = entry.replace("path:", "")
        name = desc.split(" - ")[0].strip()
        try:
            kgml = REST.kegg_get(pid, "kgml").read()
            genes = sorted({g for el in KGML_parser.read(kgml).genes for g in el.name.split()
                            if g.startswith("hsa:")})
        except Exception as e:
            print(f"  skip {pid} {name}: {e}", flush=True)
            continue
        if genes:
            panel[name] = genes
        if i % 25 == 0:
            print(f"[{i}/{len(lines)}] {name}: {len(genes)} genes", flush=True)
        time.sleep(0.1)
    json.dump(panel, open(OUT, "w"))
    sizes = sorted(len(v) for v in panel.values())
    print(f"[done] {len(panel)} pathways -> {OUT}; size min/median/max = "
          f"{sizes[0]}/{sizes[len(sizes)//2]}/{sizes[-1]}")


if __name__ == "__main__":
    main()
