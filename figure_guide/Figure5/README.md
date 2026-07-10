# Figure 5 - prostate spatial (Xenium Prime)

Spatial transcriptomics of a human prostate adenocarcinoma sample (single-cell
spatial): domain / subdomain detection, cell-type enrichment, and pathway maps.

**Data:** 10x Genomics Xenium Prime FFPE human prostate
(<https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-prostate>).
Run in a local R/Seurat environment; see `docs/r_environment.md`.
Author-confirmed: **Figure 5 = prostate**.

## Manuscript items

| Panel | Backend scripts | Required inputs | Expected outputs | Status |
|-------|-----------------|-----------------|------------------|--------|
| Fig 5 | `analysis/figure4_5_spatial_crc_prostate/analysis/prostate/scripts/00_...10_*.R` | prostate raw (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, `spatial_location.csv`) | spatial panels for `Figure5.jpeg` | verified; external data required |
| `Supp_prostate_pathways` | same prostate R pipeline | same prostate raw inputs | `Supp_prostate_pathways.jpeg` panels | likely; external data required |

Shared helper functions: `analysis/figure4_5_spatial_crc_prostate/R/`.

## Run

```bash
cd analysis/figure4_5_spatial_crc_prostate/analysis/prostate/scripts
Rscript 00_config_prostate.R
Rscript 01_process_raw_data.R
# continue numbered scripts through 10 in order
```

## Makefile target

- None. This is a numbered R pipeline, not a `analysis/shared/` Makefile target.

## Notes

- No `run.sh` wrapper is provided because local raw-data paths and R/Seurat
  state must be configured before running.
- Final `Figure5.jpeg` assembly from spatial panels is a manual layout step
  outside version control.
