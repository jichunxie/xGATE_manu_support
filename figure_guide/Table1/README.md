# Table 1 - computational complexity + edge density

Aggregated runtime and peak memory of co-expression-graph construction across 162
cell subpopulations from seven studies, with the reviewer-requested edge-density
column. Shows runtime/memory scale with candidate gene pairs (about |V|^2) and cell
count, not edge density.

Run from an activated local Python environment created from `environment.yml`,
`requirements.txt`, or the exact snapshot in `envs/`.

## Manuscript item

| Item | Backend script | Required input | Expected outputs | Command |
|------|----------------|----------------|------------------|---------|
| Table 1 | `analysis/table1_complexity/table1_runtime.py` | `computational_benchmark_summary.csv` plus `computational_benchmark_summary_additions.csv` (**not shipped**) | `results/table1_runtime_table.csv`, `figures/table1_complexity_analysis.{png,pdf}` | `make complexity` |

Outputs are generated and git-ignored by policy.

## Input (not shipped)

The two `computational_benchmark_summary*.csv` files hold the measured runtime and
peak-memory numbers for graph construction across the 162 subpopulations. They are
**not committed** to the repository. Regenerate them by running the graph-construction
timing pipeline (`analysis/data_prep/`) across your subpopulations, or request them
from the authors, then place `computational_benchmark_summary.csv` (and
`computational_benchmark_summary_additions.csv`) at the repository root before running
`make complexity`. Without them, `make complexity` will exit with a file-not-found error.

## Run

```bash
make complexity
```

or:

```bash
bash figure_guide/Table1/run.sh
```

## Makefile target

- `make complexity`

## Wrapper

- `figure_guide/Table1/run.sh` changes to the repository root and runs
  `make complexity`.

Reference output: CPU vs |V|² Spearman rho=+0.954 (p=1.9e-85); CPU vs edge
density rho=+0.083 (p=0.29); peak memory vs |V|² rho=+0.917 (p=7.1e-66).
Also produces `Supp_complexity_analysis`.
