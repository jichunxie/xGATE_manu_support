# data/

Input data are **not committed** to this repository because raw single-cell and
spatial matrices are large. This directory is a placeholder; download the public
datasets below and place them here (or point
`configs/paths.yaml:data_root` at wherever you store them).

See `docs/data_availability.md` for manuscript accessions and source portals.

## Expected layout

The biology notebooks under `analysis/manu_support_meta/` resolve their inputs from
this data root (paths were rewired from private cluster locations to a documented
placeholder root; see `configs/README.md`). Expected files:

```
data/
├── pancreas_human.rds                          # GSE148073 pancreatic beta-cell counts (Seurat)
├── pancreas_ctrl_final.csv                      # exported control beta-cell expression, genes x cells
├── pancreas_t1d_final.csv                       # exported T1D beta-cell expression, genes x cells
├── adj_matrix_pancreas_ctrl_final.csv           # co-expression adjacency (healthy)
├── adj_matrix_pancreas_t1d_final.csv            # co-expression adjacency (T1D)
├── adj_matrix_pancreas_aab_cluster1_final.csv   # AAb+ cluster 1 adjacency
├── adj_matrix_pancreas_aab_cluster2_final.csv   # AAb+ cluster 2 adjacency
├── pathway_results_aab_cluster1.csv             # AAb pathway results (cluster 1)
├── pathway_results_aab_cluster2.csv             # AAb pathway results (cluster 2)
├── hepatocyte_human.rds                         # liver hepatocyte counts (Seurat)
├── pathway_genes_master_list.json               # KEGG-derived pathway gene sets (root list)
├── pancreas_genesets_ensembl.json               # pancreas gene sets (Ensembl IDs)
├── liver_genesets_ensembl.json                  # liver gene sets (Ensembl IDs)
└── (benchmark inputs: FUCCI GSE146773, Tabula Sapiens fibroblasts, hallmark, kegg)
```

Regenerated intermediates used by the revision benchmark (e.g.
`*_competing_percall.csv`, `*_xgate_*.csv`) are produced by the pipeline into
`results/` — they do **not** need to be downloaded.

## Where to obtain the public datasets

| Dataset | Accession / source | Status |
|---------|--------------------|--------|
| Liver hepatocytes | GEO **GSM6058681** | public |
| Pancreatic β-cells (healthy / T1D / AAb+) | GEO **GSE148073** | public |
| ETO senescence time-course | GEO **GSE226225** | public |
| FUCCI U2OS cell-cycle | GEO **GSE146773** | public |
| Tabula Sapiens fibroblasts | Tabula Sapiens consortium — <https://tabula-sapiens.sf.czbiohub.org> | public |
| CRC spatial (Fig 4, Visium HD) | 10x Genomics — [human CRC Visium HD](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc) | public |
| CRC scRNA reference (Fig 4 IRIS deconvolution) | Guimarães et al. 2024, Nat. Commun. 15:5694 (CZI CELLxGENE) — [doi:10.1038/s41467-024-49916-4](https://doi.org/10.1038/s41467-024-49916-4) | public |
| Prostate spatial (Fig 5, Xenium Prime) | 10x Genomics — [human prostate Xenium Prime](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-prostate) | public |
| Prostate Cell Atlas scRNA reference (Fig 5 SingleR/IRIS) | Tuong et al. 2021, Cell Reports 37(12):110132 — [doi:10.1016/j.celrep.2021.110132](https://doi.org/10.1016/j.celrep.2021.110132) | public |
| MSigDB Hallmark gene sets | MSigDB (`h.all`) | public |
| KEGG pathways (hsa) | KEGG | public |

See `docs/data_availability.md` for the same list plus the ethics statement.

## Not included / why

- **Raw matrices** (`*.rds`, `*.h5ad`, `*.mtx`, `*.pkl`, embedding CSVs): excluded
  by size. All are publicly downloadable from the accessions above.
- **Regenerated outputs**: produced by the scripts; not shipped.
