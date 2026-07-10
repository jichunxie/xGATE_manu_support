# Data availability

This document lists the datasets used in the xGATE manuscript and their public
sources.

All datasets are **publicly available**; accessions below are taken from the
manuscript's Data Availability statement. No dataset is committed to this
repository (size); download from the sources below into `data/` (see
`data/README.md`).

## Datasets

| Dataset | Used in | Accession / source | Status |
|---------|---------|--------------------|--------|
| Liver hepatocytes | Fig 3 liver panels; `Supp_liver_*` | GEO **GSM6058681** | public |
| Pancreatic β-cells (healthy / T1D / AAb+) | Fig 2; Fig 3 pancreas panels | GEO **GSE148073** | public |
| ETO senescence time-course | Fig 3e; `Supp_senescence*` | GEO **GSE226225** | public |
| FUCCI U2OS cell-cycle | Fig 3 benchmark; `Supp_fucci_benchmark` | GEO **GSE146773** | public |
| Tabula Sapiens fibroblasts | Fig 3 benchmark; batch/donor (`Supp_batch_ts`) | Tabula Sapiens consortium — <https://tabula-sapiens.sf.czbiohub.org> | public |
| CRC spatial (Visium HD) | Fig 4; `Supp_CRC_pathways` | 10x Genomics — <https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc> | public |
| CRC scRNA reference (IRIS domain deconvolution, Fig 4) | Fig 4 (IRIS reference) | Guimarães et al. 2024, *Nat. Commun.* 15:5694, via CZI CELLxGENE — <https://doi.org/10.1038/s41467-024-49916-4> | public |
| Prostate adenocarcinoma spatial (**Xenium Prime**) | Fig 5; `Supp_prostate_pathways` | 10x Genomics — <https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-prostate> | public |
| Prostate Cell Atlas scRNA reference (SingleR annotation + IRIS, Fig 5) | Fig 5 (annotation/IRIS reference) | Tuong et al. 2021, *Cell Reports* 37(12):110132 — <https://doi.org/10.1016/j.celrep.2021.110132> | public |
| MSigDB Hallmark gene sets | hallmark benchmark variant | MSigDB `h.all` | public |
| KEGG pathways (hsa) | pathway gene sets | KEGG | public |

> Note: the prostate sample is a **Xenium Prime** dataset (single-cell spatial),
> not Visium HD; the CRC sample is Visium HD.

## Ethics

The study analyzed only publicly available, de-identified human data; its
secondary use was reviewed by the Duke University Health System IRB and granted
exempt status (see manuscript).

## Benchmark truth labels

Active/inactive pathway labels per dataset are fixed by external biology (not by
a disease-vs-normal contrast). The label lists are given in the supplementary
tables and Methods of the manuscript.

## What is and is not in this repository

- **Included:** all analysis code (Python, R, notebooks), environment
  specifications, and documentation.
- **Not included:** raw and large processed matrices (`*.rds`, `*.h5ad`,
  `*.mtx`, `*.pkl`, embedding CSVs) and all regenerated outputs. Obtain the raw
  inputs from the sources above and place them under `data/` (see
  `data/README.md`); regenerate outputs with the scripts.

## Deposition blockers

**None** — every dataset is publicly available with an accession above. No
deposition is pending.
