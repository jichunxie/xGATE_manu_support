# Figure 4 - CRC spatial (Visium HD)

Spatial transcriptomics of a human colorectal cancer Visium HD sample: domain
detection, sub-clustering, immune-cold/-hot pathway correlation heatmap, and ROI
pathway-activity maps.

**Data:** 10x Genomics human CRC Visium HD
(<https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc>).
Run in a local R/Seurat environment; see `docs/r_environment.md`.

## Manuscript items

| Panel | Backend scripts | Required inputs | Expected outputs |
|-------|-----------------|-----------------|------------------|
| Fig 4 (a-f) | `analysis/figure4_5_spatial_crc_prostate/analysis/crc/scripts/00_...12_*.R` | CRC raw (`countList_spatial_crc.RDS`, `sc_input_crc.RDS`) | spatial panels for `Figure4.jpeg` |
| `Supp_CRC_pathways` | same CRC R pipeline | same CRC raw inputs | `Supp_CRC_pathways.jpeg` panels |

Shared helper functions: `analysis/figure4_5_spatial_crc_prostate/R/`.

## Run

```bash
cd analysis/figure4_5_spatial_crc_prostate/analysis/crc/scripts
Rscript 00_config.R
Rscript 01_identify_spatial_domains.R
# continue numbered scripts through 12 in order
```

## Makefile target

- None. This is a numbered R pipeline, not a `analysis/shared/` Makefile target.

## Notes

- No `run.sh` wrapper is provided because local raw-data paths and R/Seurat
  state must be configured before running.
- Final `Figure4.jpeg` assembly from spatial panels is a manual layout step
  outside version control.
