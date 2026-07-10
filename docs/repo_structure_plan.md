# Repository structure

This repository is organized for manuscript reproducibility. Browse by figure/table in
`figure_guide/`; the backend code lives in `analysis/`, grouped by purpose.

## Layout

| Path | Contents |
|---|---|
| `analysis/shared/` | Path resolver (`xgate_paths.py`), plot style, launcher (`xgate_py.sh`), and shared helper modules imported across analyses. |
| `analysis/data_prep/` | Dataset extraction, normalization, and graph-building scripts (`extract_*`, `build_*`, `run_xgate.py`, `ts_*`, `panc_donor_graph.py`, `fucci_build.py`, `*_export.R`). |
| `analysis/competing_methods/` | Competing-method benchmark: `competing_methods.R`, `run_all_competing_metrics.R` (ORA/AUCell/scGSEA/GESECA/PAGODA), and `prep_competing_*.py` input prep. |
| `analysis/figure2_pancreas/` | Figure 2 pancreas analyses (co-expression graphs, MDS embeddings, AAB) + Figure 2b re-code inputs. |
| `analysis/figure3_benchmark_senescence_aab/` | Figure 3 benchmark metrics, FDR/BH, Jaccard-vs-p-value, AAB notebooks, and the senescence time-course notebook in its `senescence/` subfolder (Fig 3e). |
| `analysis/figure4_5_spatial_crc_prostate/` | Spatial R pipelines: CRC (Figure 4) and prostate (Figure 5). Kept as one unit because both share the helper functions in its `R/` folder via relative `source("R/…")`. |
| `analysis/table1_complexity/` | Table 1 computational-complexity / runtime analysis (`table1_runtime.py`; 162 clusters / 7 studies). |
| `analysis/supplementary/` | Supplementary panels: batch-effect (`supp_figS3_5_batch_figures.py`, canonical), matched-null p-values, partial-activity, threshold sensitivity, plus `misc/supp_figS20_pathway_comptime.ipynb` (S20 time-complexity, Jupyter). |
| `analysis/manu_support_meta/` | Auxiliary files from the original manuscript-support tree (`requirements.txt`, `verify_xgate.py`). |
| `figure_guide/` | Figure/table-level entry points and readable analysis descriptions. |
| `docs/` | Reproducibility, paper→code map, data availability, R environment, dependency plan, script name map. |
| `configs/`, `data/`, `tests/`, `envs/` | Config templates, data-download instructions, sanity checks, conda snapshots. |
| `figures/`, `results/` | Generated locally by scripts; **git-ignored** (only `.gitkeep` tracked). Not committed. |
| `computational_benchmark_summary*.csv` | Table 1 timing/memory measurements read by `table1_runtime.py`. **Not shipped** — regenerate from the graph-construction timing pipeline (`analysis/data_prep/`) or request from the authors; place at the repo root before `make complexity`. |

## Script names

Figure/table scripts are named by the manuscript item they produce:

- `figXY_*` — main figure X, panel Y (e.g. `fig3c_jaccard_pvalue.py`).
- `supp_figSN_*` — supplementary figure SN (e.g. `supp_figS10_threshold_sensitivity.py`).
- `table1_*` — Table 1.
- `figX_build_*` / `supp_figSN_build_*` — upstream builder feeding that figure (no panel itself).

Shared library modules (`analysis/shared/`), dataset prep (`analysis/data_prep/`), and
competing-method scripts (`analysis/competing_methods/`) keep descriptive names. Users navigate
via **Makefile targets**, **`figure_guide/`**, and **`docs/script_name_map.md`**, which map each
script to its manuscript analysis. (Some scripts still emit intermediate output files under
their old `r*` stems; those output names are unchanged in this pass — see the note in
`docs/script_name_map.md`.)

## Path resolution

`analysis/shared/xgate_py.sh` puts `analysis/shared` on `PYTHONPATH` and exports
`XGATE_ROOT`; the **xGATE method** is imported from the installed package
(`pip install git+https://github.com/jichunxie/xGATE.git@v1.0`), not the path. `analysis/shared/xgate_paths.py` resolves the repo root and the
`results/` and `figures/` output roots (config- or env-driven; no absolute paths hard-coded).
