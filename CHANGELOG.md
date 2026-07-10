# Changelog

All notable changes to this manuscript-reproducibility repository.

## [1.0.0] - 2026-07-10

Initial public release accompanying the xGATE manuscript.

### Changed
- **Script naming.** Figure/table scripts renamed to match the manuscript item they
  produce: `figXY_*` (main figure X, panel Y), `supp_figSN_*` (supplementary figure SN),
  `table1_*` (Table 1); upstream builders use `figX_build_*` / `supp_figSN_build_*`.
  Shared library, dataset-prep, and competing-method scripts keep descriptive names.
- **Output filenames aligned.** Figure/CSV outputs now use stems that match their
  producing script (e.g. `fig3c_jaccard_pvalue.pdf`, `table1_runtime_table.csv`,
  `supp_figS3_batch_metrics`, `fig3_benchmark_metrics_bh.csv`, `supp_figS6_9_benchmark_<dataset>.*`).
  One input, `supp_figS23_<dataset>_matched_null.csv`, is produced by an off-repo step — rename any
  local `r14_<dataset>_matched_null.csv` to match.
- **Repository layout.** Competing-method scripts moved from `analysis/data_prep/` to a new
  `analysis/competing_methods/`. `data_prep/` now holds only dataset extraction /
  normalization / graph building.
- **Naming.** The legacy internal name `CPACT` / `CellSubNet` fully retired in favor of
  **xGATE** across code and docs, including data filenames and column values.
- License clarified as **MIT** (`LICENSE.md`), consistent with the manuscript.

### Removed
- Retired `analysis/supplementary/make_all_figures.py` (superseded by the individual
  per-figure scripts; no code callers).
- Deleted a dead duplicate of `r14_pvalue_compare.py` that also contained hard-coded
  private paths.
- **Data policy: no data CSVs shipped.** Removed the Table 1 timing inputs
  (`computational_benchmark_summary*.csv`) and the generated `table1_runtime_table.csv`.
  Table 1 now requires the user to supply/regenerate those measurements — see
  `figure_guide/Table1/README.md`.

### Fixed
- Repaired stale `Makefile` targets and documentation references left by earlier renames
  (`benchmark_extra_plots.py`, `assemble_benchmark.py`, `figure3_senescence/`, `liver/`).

### Dependencies
- Pinned `requirements.txt` / `environment.yml` to the verified `xgate` conda env
  (Python 3.10.20; numpy 1.26.4, pandas 2.3.3, scipy 1.15.2, scikit-learn 1.7.2,
  networkx 3.4.2, matplotlib 3.10.9, torch 2.12.0, biopython 1.87, …). Note: PyYAML and
  JupyterLab are not in that env snapshot and are marked optional.

### Notes
- If reproducing from data generated under the old `CPACT` names, rename the matching
  local data files/columns to `xGATE` — see `docs/script_name_map.md`.
