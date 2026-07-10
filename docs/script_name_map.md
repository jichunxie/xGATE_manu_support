# Script name map (script → manuscript meaning)

Analysis scripts are named by the manuscript item they produce, so the repository is browsable
by figure/table. The `make` target that runs each entry-point script is shown where one exists.

**Naming convention**

- `figXY_*` — produces **main figure X, panel Y** (e.g. `fig3c_…` → Figure 3c).
- `supp_figSN_*` — produces **supplementary figure SN** (e.g. `supp_figS10_…` → Supp. Fig. S10).
- `table1_*` — produces **Table 1**.
- `figX_build_*` / `supp_figSN_build_*` — **upstream builder** that generates inputs for that
  figure but does not itself emit a panel.
- `analysis/shared/` — shared library modules; names are stable (imported by many scripts).
- `analysis/data_prep/`, `analysis/competing_methods/` — dataset preparation and competing-method
  benchmarking; descriptive names, not figure-scoped.

> Output-filename note: output PDFs/CSVs use stems that match their producing script
> (e.g. `fig3c_jaccard_pvalue.py` → `figures/fig3c_jaccard_pvalue.pdf`; `table1_runtime.py` →
> `results/table1_runtime_table.csv`; the Fig 3 benchmark table is `results/fig3_benchmark_metrics_bh.csv`;
> the per-dataset S6–S9 check/cross tables are `figures/supp_figS6_9_benchmark_<dataset>.*`). These
> figure outputs are still manually assembled/renamed to the manuscript `Figure*.jpeg` / `Supp_*.pdf`
> outside version control. One input is intentionally left under its original stem:
> `results/supp_figS23_<dataset>_matched_null.csv` is read by `supp_figS23_matched_null_pvalues.py`
> but produced by an **off-repo** degree-matched-null step — if you generated it under the old
> `r14_<dataset>_matched_null.csv` name, rename your local files to match.

## Figure / table entry points (produce panels)

| Script | Manuscript item | `make` target |
|---|---|---|
| `analysis/figure2_pancreas/fig2a_coexpr_graphs.ipynb` | Fig 2a co-expression graphs | — |
| `analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` (rendered: `fig2b_Mean_Exp.html`) | Fig 2b marginal expression (R/Seurat) | — |
| `analysis/figure2_pancreas/fig2c_embedding_comparison.ipynb` | Fig 2c MDS embedding space | — |
| `analysis/figure3_benchmark_senescence_aab/fig3a_assemble_benchmark.py` | Fig 3a,b + Supp benchmark PDFs (S6–S9) | `make benchmark` |
| `analysis/figure3_benchmark_senescence_aab/fig3b_metrics_panel.py` | Fig 3a,b benchmark metrics panel | `make fig3` |
| `analysis/figure3_benchmark_senescence_aab/fig3c_jaccard_pvalue.py` | Fig 3c + Supp_jaccard_pvalue_full (S12) | `make fdr` |
| `analysis/figure3_benchmark_senescence_aab/fig3d_fdr_4datasets.py` | Fig 3d + Supp_bh_by_simulation (S13) | `make fdr` |
| `analysis/figure3_benchmark_senescence_aab/fig3d_fdr_realized.py` | Fig 3d realized FDR (BH / BY) | — |
| `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb` | Fig 3e + Supp_senescence (S15) | — |
| `analysis/figure3_benchmark_senescence_aab/fig3g_pancreas_ctrl_t1d.ipynb` | Fig 3g + Supp_pancreas_pathway_activity (S17) | — |
| `analysis/figure3_benchmark_senescence_aab/fig3gh_aab_figure.ipynb` | Fig 3g,h AAB cluster comparison | — |
| `analysis/table1_complexity/table1_runtime.py` | Table 1 + Supp_complexity_analysis (S19); 162 clusters / 7 studies | `make complexity` |
| `analysis/supplementary/misc/supp_figS20_pathway_comptime.ipynb` | Supp_time_complexity (S20; per-pathway computation time vs pathway size) | — (Jupyter) |
| `analysis/supplementary/supp_figS1_2_partial_activity.py` | Supp_partial_activity S1 (+ `oxphos` = hub variant S2) | `make partial` |
| `analysis/supplementary/supp_figS3_5_batch_figures.py` | Supp_batch_metrics / _pancreas / _ts (S3–S5) — **canonical batch script** | `make batch` |
| `analysis/supplementary/supp_figS10_threshold_sensitivity.py` | Supp_threshold_sensitivity (S10; threshold + extra benchmark plots) | `make extra` |
| `analysis/figure3_benchmark_senescence_aab/supp_figS11_pvalue_diagnostics.py` / `supp_figS11_pvalue_plot.py` | Supp_pvalue_diagnostics (S11) | — |
| `analysis/figure3_benchmark_senescence_aab/supp_figS14_bh_asymptotic.py` | Supp_bh_by_asymptotic (S14) | — |
| `analysis/supplementary/supp_figS22_density_degree.py` | Supp_density_degree_1 / _2 (S22) | — |
| `analysis/supplementary/supp_figS23_matched_null_pvalues.py` | Supp_matched_null_pvalues (S23; default vs degree-matched null) | `make nulls` |

## Upstream builders (`figX_build_*` / `supp_figSN_build_*`; produce inputs consumed above)

| Script | Role |
|---|---|
| `analysis/figure3_benchmark_senescence_aab/fig3_build_fucci.py` | FUCCI U2OS graph/figure build (feeds Fig 3 benchmark) |
| `analysis/figure3_benchmark_senescence_aab/fig3_build_pathway_screen.py` | pathway-activity screen (feeds Fig 3) |
| `analysis/figure3_benchmark_senescence_aab/fig3gh_build_aab_generate.ipynb` | AAB cluster data generation (feeds Fig 3g,h) |
| `analysis/figure3_benchmark_senescence_aab/fig3gh_build_aab_xgate.ipynb` | run xGATE on AAB clusters (feeds Fig 3g,h) |
| `analysis/figure3_benchmark_senescence_aab/supp_figS21_build_randomnull.py` | random-null p-values (feeds Supp_two_nulls S21) |
| `analysis/supplementary/supp_figS3_5_build_batch_donor.py` / `_build_batch_manyset.py` / `_build_batch_readout.py` / `_build_readout_anova.py` / `_build_donor_classifier.py` | batch-effect input builds / donor readouts (feed S3–S5) |
| `analysis/supplementary/supp_figS22_build_distributions.py` | null degree/density distributions build (feeds S22) |
| `analysis/supplementary/supp_figS1_2_build_genemix.py` | gene-mix runs (feed S1,S2) |

## Dataset preparation (`analysis/data_prep/`)

| Script | Role |
|---|---|
| `run_xgate.py`, `ts_build.py`, `ts_extract_embeddings.py`, `ts_merge_posneg.py`, `panc_donor_graph.py`, `fucci_build.py` | dataset graph / embedding builds |
| `extract_*`, `build_from_csv.py`, `build_kegg_panel.py`, `prep_liver_pancreas_ensembl.py`, `lognorm_export.R`, `sctransform_export.R` | data extraction / normalization for benchmark datasets |

## Competing methods (`analysis/competing_methods/`)

| Script | Role |
|---|---|
| `competing_methods.R`, `run_all_competing_metrics.R` | ORA / AUCell / scGSEA runs + unified 4-dataset per-call metrics |
| `prep_competing_fucci.py` / `prep_competing_full.py` / `prep_competing_ts.py` | prep normalized inputs for the competing-method comparison |

## Shared helpers (`analysis/shared/`; imported by many scripts — names are stable)

| Module | Role |
|---|---|
| `xgate_paths.py` | config-driven path resolver (imported repo-wide) |
| `plot_style.py` | shared matplotlib styling |
| `bench_lib.py` | benchmark metric helpers |
| `build_assets.py` | rebuilds Fig 3 / Supp benchmark composite assets after BH q<0.05 calls (imported by `fig3d_fdr_realized.py`) |
| `matched_null.py` | degree-matched null generation |
| `xgate_py.sh` | launcher (sets `PYTHONPATH` to `analysis/shared`; xGATE imported from the installed package) |

> Naming note: the method is **xGATE**. The earlier internal names `CPACT` / `CellSubNet`
> have been fully retired from the code and docs. Some **input/output data filenames and
> data-column values** were renamed to match (e.g. `*_xGATE_pancreas_embeddings.csv`,
> `..._xgate_subdomain.csv`, the `Type="xGATE"` column, and `AAB…xGATE.csv`). If you are
> reproducing from data generated under the old names, rename those local files/columns from
> `CPACT`/`cpact` to `xGATE`/`xgate` so the scripts find them. Affected artifacts:
> `*_CPACT_pancreas_embeddings.csv` → `*_xGATE_pancreas_embeddings.csv`;
> `colon_pathway_analysis_results_cpact_subdomain.csv` → `..._xgate_subdomain.csv`;
> `AAB…CPACT.csv` → `AAB…xGATE.csv`; and the `Type="CPACT"` column value → `Type="xGATE"`.
