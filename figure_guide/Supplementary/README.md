# Supplementary figures - by producer

All `Supp_*` panels grouped by producing code. Full per-panel table (inputs,
outputs, status) is in [`docs/paper_code_map.md`](../../docs/paper_code_map.md).
Python panels run through `make <target>` after activating your local Python
environment; biology/spatial panels use local Jupyter or R environments.

## Status summary

- Verified runnable with Makefile targets: benchmark extras, null model panels,
  complexity panels, partial-activity panels.
- External/generated inputs required: batch/donor panels, biology notebooks,
  CRC/prostate spatial panels.
- Source directory confirmed: `Supp_Treg`, `Supp_pca_immune_fraction` are produced by the
  spatial pipeline under `analysis/figure4_5_spatial_crc_prostate/`.
- To be confirmed: exact `Supp_Senescence_reclustering` producing cell/script.

## Benchmark panels
`analysis/figure3_benchmark_senescence_aab/fig3a_assemble_benchmark.py` + `supp_figS10_threshold_sensitivity.py`
produce `Supp_{liver,pancreas,fucci,ts}_benchmark` and
`Supp_threshold_sensitivity`.
Command: `make benchmark extra`.

## FDR / null-dependence panels
`analysis/shared/r12_*` produces `Supp_jaccard_pvalue_full`, `Supp_pvalue_diagnostics`,
`Supp_bh_by_simulation`, `Supp_bh_by_asymptotic`. Command: `make fdr`.

## Batch / donor panels
`analysis/supplementary/supp_figS3_5_batch_figures.py` (the **canonical** batch-figure script) produces
`Supp_batch_metrics`, `Supp_batch_pancreas`, and `Supp_batch_ts`. Command: `make batch`.
Required input: upstream embeddings in `results/`.
Older batch variants have been retired to `archive/old_numbering_or_obsolete/batch_superseded/`.

## Null-model panels
`analysis/supplementary/supp_figS23_matched_null_pvalues.py` (+ `matched_null.py`, `supp_figS22_density_degree.py`)
produces `Supp_matched_null_pvalues` / `Supp_two_nulls`,
`Supp_density_degree_1/2`.
Command: `make nulls`.

## Complexity panels
`analysis/table1_complexity/table1_runtime.py` produces `Supp_complexity_analysis` (S19).
Command: `make complexity`. (See [Table1](../Table1/README.md).)
`analysis/supplementary/misc/supp_figS20_pathway_comptime.ipynb` produces
`Supp_time_complexity` (S20; per-pathway comp-time vs pathway size). Local Jupyter;
reads `misc/Results/pathway_comp_time_results.csv`.

## Partial-activity panels
`analysis/supplementary/supp_figS1_2_partial_activity.py` (+ `oxphos`) produces
`Supp_partial_activity(.hub)`.
Command: `make partial`.

## Biology panels
- `Supp_senescence`, `Supp_Senescence_reclustering` — `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb` (reclustering panel: exact cell **TBD**).
- `Supp_pancreas_pathway_activity` — `analysis/figure3_benchmark_senescence_aab/fig3g_pancreas_ctrl_t1d.ipynb`.

## Spatial panels
- `Supp_CRC_pathways` — CRC pipeline (see [Figure4](../Figure4/README.md)).
- `Supp_prostate_pathways` — prostate pipeline (see [Figure5](../Figure5/README.md)).
- `Supp_Treg`, `Supp_pca_immune_fraction` — spatial pipeline under
  `analysis/figure4_5_spatial_crc_prostate/` (R). Source directory confirmed;
  exact script-to-panel call may be verified later if needed.

## Competing-method benchmarks
`analysis/competing_methods/competing_methods.R` + `run_all_competing_metrics.R` run the
unified 4-dataset ORA / AUCell / scGSEA / GESECA / PAGODA comparison; normalized inputs are
prepared by `analysis/competing_methods/prep_competing_*.py`.

## Wrappers

- No supplementary `run.sh` wrapper is added. Use the Makefile targets above so
  shared backend scripts stay in `analysis/shared/`.
