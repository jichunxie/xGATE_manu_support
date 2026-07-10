# Figure guide - by-figure index

This directory is the **figure-level frontend/navigation layer** for the
manuscript-support repository. It tells a reader what to run for each figure or
table while keeping the shared analysis backend intact.

```
analysis/shared/
  Shared backend scripts for revised analyses. Do not physically split these by
  figure; same-directory imports, xgate_paths.py, xgate_py.sh, and Makefile
  targets depend on this layout.

figure_guide/
  Figure-level frontend/navigation layer. Each figure/table folder explains
  panels covered, status, backend scripts/notebooks, required inputs, expected
  outputs, commands, and unresolved issues.
```

Intended use: choose a manuscript item, open its folder, read `README.md`, then
run the listed Makefile target, wrapper, notebook, or R pipeline when inputs are
available. The guide does not duplicate scientific logic.

For the single flat table (all figures + supplementary in one view) see
[`docs/paper_code_map.md`](../docs/paper_code_map.md). Reproduction commands and
environment setup are in [`docs/reproducibility.md`](../docs/reproducibility.md).

## Main figures & table

| Item | Guide | What it shows | Producer (summary) |
|------|-------|---------------|--------------------|
| Figure 1 | [Figure1/](Figure1/README.md) | Pipeline schematic | drawn (no code) |
| Figure 2 | [Figure2/](Figure2/README.md) | Pancreas: co-expression graphs, marginal expression, embeddings | biology notebooks plus a parameterized Fig 2b re-code |
| Figure 3 | [Figure3/](Figure3/README.md) | Benchmark (a,b), FDR/Jaccard (c,d), senescence (e) | `analysis/shared/` + senescence notebook |
| Figure 4 | [Figure4/](Figure4/README.md) | CRC spatial (Visium HD) | CRC R pipeline |
| Figure 5 | [Figure5/](Figure5/README.md) | Prostate spatial (Xenium Prime) | prostate R pipeline |
| Table 1 | [Table1/](Table1/README.md) | Computational complexity + edge density | `analysis/table1_complexity/table1_runtime.py` |
| Supplementary | [Supplementary/](Supplementary/README.md) | All `Supp_*` panels | mixed (`analysis/shared/` + biology) |

## Why an index instead of moving files

`analysis/shared/` is a **shared analysis library**: 61 scripts import the common
`xgate_paths` module, the benchmark chain feeds several figures at once, and the
`Makefile`/launcher assume that layout. Physically splitting it per figure would
break those imports and the reproducible pipeline. This index gives the by-figure
organization without that risk (author-approved "Option A"; see
[`docs/restructure_plan.md`](../docs/restructure_plan.md)).
