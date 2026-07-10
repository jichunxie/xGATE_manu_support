# R environment

The R analyses (competing-method benchmarks and the spatial pipelines) use CRAN +
Bioconductor + Seurat. A full R lockfile is not provided; instead the exact conda
snapshot `envs/xgate_r.yml` pins every R and Bioconductor package that was used.

## Recommended: conda snapshot

Choose any environment name; the name below is an example, not an author cluster
environment.

```bash
conda env create -f envs/xgate_r.yml -n xgate-manu-r
conda activate xgate-manu-r
```

Key pinned versions (see `envs/xgate_r.yml` for the complete list):

| Package | Version |
|---------|---------|
| R (r-base) | 4.4.3 |
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| AUCell | 1.28.0 |
| clusterProfiler (ORA) | 4.14.0 |
| fgsea (GESECA) | 1.32.2 |
| igraph (R) | 2.1.4 |

## Alternative: install into an existing R

If you prefer your own R (≥ 4.1) rather than the conda snapshot, install the
packages the scripts need:

```r
# CRAN
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "viridis", "jsonlite"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "AUCell", "clusterProfiler", "fgsea",
  "SingleR", "SingleCellExperiment", "SummarizedExperiment",
  "scuttle", "BiocParallel", "ComplexHeatmap"
))
```

The spatial pipelines additionally use the helper functions bundled in
`analysis/figure4_5_spatial_crc_prostate/R/`; the dependency
installer there is `install_dependencies.R`.

## Notes

- `scGSEA` (via `gficf`) and `PAGODA` (`pagoda2`) can be difficult to build on
  some systems (Rook / C++ toolchain). They are only needed to regenerate the
  scGSEA and PAGODA benchmark columns; the assembled metric CSVs already contain
  those columns.
- Versions here are what was used for the manuscript; newer patch releases are
  generally fine. Pin exactly only if you need bit-for-bit reproduction.
