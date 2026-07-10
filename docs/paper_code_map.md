# Paper → code map

Maps each manuscript item to the script/notebook that produces it, its required
input, and its expected output. Status legend:

- **verified** — producing script traced and confirmed in this repo.
- **likely** — mapping is strongly supported but not independently re-run here.
- **re-coded** — original plotting code was lost; a transparent replacement
  script exists and uses manuscript-matched inputs.
- **needs confirmation** — provenance uncertain; author to confirm.

Paths are relative to the repository root. Revision panels are written to
`figures/` with the internal names shown; the manuscript `Figure*.jpeg`
and `Supp_*.pdf` files are assembled/renamed from these panels (a manual step
outside version control).

## Main figures

| Manuscript item | Script / notebook | Required input | Expected output | Status | Notes |
|---|---|---|---|---|---|
| **Fig 1** pipeline schematic | — (drawn, not code-generated) | — | `Figure1.jpeg` | verified | Illustration; no producing script. |
| **Fig 2a** co-expression graphs (active/inactive pathways) | `analysis/figure2_pancreas/fig2a_coexpr_graphs.ipynb` | `data/adj_matrix_pancreas_ctrl_final.csv`, pathway gene sets | graph panels → `Figure2.jpeg` | likely | Uses healthy-pancreas co-expression adjacency. |
| **Fig 2b** marginal gene expression boxplots | `analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` (render to HTML locally) | `data/pancreas_human.rds` (GSE148073 Seurat object), KEGG pathway genes (via `KEGGREST`/`org.Hs.eg.db`) | boxplot JPGs → `Figure2.jpeg` | verified | R/Seurat notebook: selects control `type B pancreatic cell`s, `SCTransform(vst.flavor="v2")`, one boxplot per pathway (same pathways as Fig 2a, e.g. `AMPK signaling pathway`, `Bacterial invasion of epithelial cells`). Supersedes the earlier lost Python recode. |
| **Fig 2c** MDS embedding space | `analysis/figure2_pancreas/fig2c_embedding_comparison.ipynb` | `Pancreas/Embeddings/*_xGATE_pancreas_embeddings.csv` (+ FeatherGraph/Graph2Vec/NetLSD) | MDS scatter → `Figure2.jpeg` | likely | Embedding CSVs are large; regenerable, git-ignored. |
| **Fig 3a,b** precision–recall / specificity / MCC / AUCPR across 4 datasets | `analysis/figure3_benchmark_senescence_aab/fig3a_assemble_benchmark.py` → `analysis/figure3_benchmark_senescence_aab/fig3b_metrics_panel.py` | `results/fig3_benchmark_metrics_bh.csv` (from per-method `*_percall.csv`) | `figures/fig3b_metrics_panel.pdf` → `Figure3.jpeg` | verified | BH-corrected xGATE calls; scGSEA-AUCPR NA-guard fix baked in. `make fig3`. |
| **Fig 3c** null p-value vs gene-set Jaccard overlap | `analysis/figure3_benchmark_senescence_aab/fig3c_jaccard_pvalue.py` | per-dataset inactive-pathway p-values | `figures/fig3c_jaccard_pvalue.pdf` → `Figure3.jpeg` | verified | `make fdr`. |
| **Fig 3d** realized FDR (BH / BY) | `analysis/figure3_benchmark_senescence_aab/fig3d_fdr_4datasets.py` (`fig3d_fdr_realized.py`) | `results/fig3_benchmark_metrics_bh.csv` + labels | `figures/fig3d_fdr_4datasets.pdf` → `Figure3.jpeg` | verified | `make fdr`. |
| **Fig 3e** senescence pathway trends vs marker genes | `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb` | `Senescence Study/Results/final_analysis.csv` | trend line plot → `Figure3.jpeg` | verified | ETO 10-day time-course. |
| **Fig 4** CRC spatial (domains, subclusters, correlation heatmap, ROIs) | `analysis/figure4_5_spatial_crc_prostate/analysis/crc/scripts/00_…12_*.R`; xGATE subdomain scores via `analysis/figure4_5_spatial_crc_prostate/run_xgate_spatial.py` | CRC Visium HD counts (`countList_spatial_crc.RDS`, `sc_input_crc.RDS`; built from the 10x raw per `Python_crc.txt` / bin2cell) + CRC scRNA IRIS reference (Guimarães 2024, CELLxGENE) | spatial panels → `Figure4.jpeg` | verified | Numbered R pipeline; shared funcs in `.../R/`. Pathway CSV `colon_pathway_analysis_results_xgate_subdomain.csv` produced by `run_xgate_spatial.py`. |
| **Fig 5** prostate spatial (Xenium Prime) | `analysis/figure4_5_spatial_crc_prostate/analysis/prostate/scripts/00_…10_*.R`; xGATE subdomain scores via `analysis/figure4_5_spatial_crc_prostate/run_xgate_spatial.py` | prostate raw (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, `spatial_location.csv`; 10x Xenium Prime) + Prostate Cell Atlas scRNA reference (Tuong 2021, SingleR/IRIS) | spatial panels → `Figure5.jpeg` | verified | **Author-confirmed Fig 5 = prostate spatial.** `01_process_raw_data.R` builds `sp_input` from the 10x raw. Pathway CSV `prostate_pathway_analysis_results.csv` produced by `run_xgate_spatial.py`. |

## Table

| Manuscript item | Script | Required input | Expected output | Status | Notes |
|---|---|---|---|---|---|
| **Table 1** computational complexity + edge density | `analysis/table1_complexity/table1_runtime.py` | `computational_benchmark_summary.csv`, `computational_benchmark_summary_additions.csv` (**not shipped** — user-supplied timing measurements; regenerate or request) | `results/table1_runtime_table.csv`, `figures/table1_complexity_analysis.*` | verified | `make complexity` (needs the input at repo root). Input covers 162 clusters / 7 studies; edge-density column added per reviewer. |

## Supplementary figures

| Supp figure | Script | Input | Output (internal name) | Status |
|---|---|---|---|---|
| `Supp_partial_activity(.hub)` | `analysis/supplementary/supp_figS1_2_partial_activity.py` (+ `oxphos` arg) | gene-mix runs | `supp_figS1_2_partial_activity*.pdf` | verified |
| `Supp_batch_metrics` | `analysis/supplementary/supp_figS3_5_batch_figures.py` | `results/supp_figS3_5_batch_*` | `supp_figS3_batch_metrics.pdf` | verified |
| `Supp_batch_pancreas` | `analysis/supplementary/supp_figS3_5_batch_figures.py` | `supp_figS3_5_batch_manyset_pancreas*` | `supp_figS4_batch_pancreas.pdf` | verified |
| `Supp_batch_ts` | `analysis/supplementary/supp_figS3_5_batch_figures.py` | `supp_figS3_5_batch_manyset_ts*` | `supp_figS5_batch_ts.pdf` | verified |
| `Supp_liver_benchmark_pathway_activity` | `analysis/figure3_benchmark_senescence_aab/fig3a_assemble_benchmark.py` / `analysis/supplementary/supp_figS10_threshold_sensitivity.py` | liver per-call CSVs | `supp_figS6_9_benchmark_Liver.pdf` | verified |
| `Supp_pancreas_benchmark_pathway_activity` | same | pancreas per-call CSVs | `supp_figS6_9_benchmark_Pancreas.pdf` | verified |
| `Supp_fucci_benchmark` | same | FUCCI per-call CSVs | `supp_figS6_9_benchmark_FUCCI_U2OS.pdf` | verified |
| `Supp_ts_benchmark` | same | TS per-call CSVs | `supp_figS6_9_benchmark_TS_Fibroblast.pdf` | verified |
| `Supp_threshold_sensitivity` | `analysis/supplementary/supp_figS10_threshold_sensitivity.py` | benchmark metrics | `supp_figS10_threshold_sensitivity.pdf` | verified |
| `Supp_pvalue_diagnostics` | `analysis/figure3_benchmark_senescence_aab/supp_figS11_pvalue_diagnostics.py` / `supp_figS11_pvalue_plot.py` | `supp_figS11_pancreas_ctrl_focused.csv` | `supp_figS11_pvalue_diagnostics.pdf` | likely |
| `Supp_jaccard_pvalue_full` | `analysis/figure3_benchmark_senescence_aab/fig3c_jaccard_pvalue.py` | all-pairs p-values | `fig3c_jaccard_pvalue.pdf` (full) | verified |
| `Supp_bh_by_simulation` | `analysis/figure3_benchmark_senescence_aab/fig3d_fdr_4datasets.py` / `supp_figS14_bh_asymptotic.py` | simulated FDR draws | `fig3d_fdr_*` | likely |
| `Supp_bh_by_asymptotic` | `analysis/figure3_benchmark_senescence_aab/supp_figS14_bh_asymptotic.py` | benchmark labels | `fig3d_fdr_*` | likely |
| `Supp_senescence` | `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb` | senescence results | senescence supp panels | likely |
| `Supp_Senescence_reclustering` | senescence notebook (reclustering section) | senescence counts | reclustering panel | needs confirmation |
| `Supp_pancreas_pathway_activity` | `analysis/figure3_benchmark_senescence_aab/fig3g_pancreas_ctrl_t1d.ipynb` | pancreas counts + gene sets | pancreas pathway panel | likely |
| `Supp_Treg` | `analysis/figure4_5_spatial_crc_prostate/` (R spatial pipeline) | spatial data | `Supp_Treg` panel | mapped — author-confirmed source directory; exact script-to-panel call should be verified if needed |
| `Supp_complexity_analysis` | `analysis/table1_complexity/table1_runtime.py` | `computational_benchmark_summary.csv`, `computational_benchmark_summary_additions.csv` | `table1_complexity_analysis.pdf` | verified |
| `Supp_time_complexity` | `analysis/supplementary/misc/supp_figS20_pathway_comptime.ipynb` | `analysis/supplementary/misc/Results/pathway_comp_time_results.csv` | per-pathway comp-time vs pathway-size panel | likely (Jupyter) |
| `Supp_two_nulls` / `Supp_matched_null_pvalues` | `analysis/supplementary/supp_figS23_matched_null_pvalues.py` (+ `matched_null.py`) | pancreas-ctrl default vs degree-matched null | `Supp_matched_null_pvalues.pdf` | verified (`make nulls`) |
| `Supp_density_degree_1` / `_2` | `analysis/supplementary/supp_figS22_density_degree.py` | null distributions | `supp_figS22_density_degree_*.pdf` | likely |
| `Supp_CRC_pathways` | CRC spatial R pipeline (Fig 4 scripts) | CRC raw | `Supp_CRC_pathways.jpeg` | likely |
| `Supp_prostate_pathways` | prostate spatial R pipeline (Fig 5 scripts) | prostate raw | `Supp_prostate_pathways.jpeg` | likely |
| `Supp_pca_immune_fraction` | `analysis/figure4_5_spatial_crc_prostate/` (R spatial pipeline) | spatial data | `Supp_pca_immune_fraction` panel | mapped — author-confirmed source directory; exact script-to-panel call should be verified if needed |

## Benchmarking (competing methods)

| Manuscript item | Script | Input | Output | Status |
|---|---|---|---|---|
| GESECA benchmark | `analysis/competing_methods/competing_methods.R` (`geseca` via `fgsea`) | pathway gene sets + counts | `results/*_competing_percall.csv` (GESECA rows) | verified (unified pipeline; standalone `GESECA pathway analysis.Rmd` superseded, not in export) |
| PAGODA benchmark | `analysis/competing_methods/competing_methods.R` (`pagoda2`) | counts | `results/*_competing_percall.csv` (PAGODA rows) | verified (unified pipeline; standalone `Pagoda pathway analysis.Rmd` superseded, not in export) |
| ORA / AUCell / scGSEA (pancreas) | `analysis/competing_methods/competing_methods.R` + `run_all_competing_metrics.R` | pancreas counts + gene sets | `results/*_competing_percall.csv` | verified (unified pipeline; standalone `Competing_methods_ORA_AUCell_scGSEA.Rmd` superseded, not in export) |
| ORA / AUCell / scGSEA (liver) | `analysis/competing_methods/competing_methods.R` + `run_all_competing_metrics.R` (unified 4-dataset pipeline; see row below) | liver counts + gene sets | per-method pathway CSVs | verified (standalone liver Rmd superseded by the unified pipeline; not in export) |
| 4-dataset competing methods (revision) | `analysis/competing_methods/competing_methods.R`, `run_all_competing_metrics.R` | per-dataset counts | `results/*_competing_percall.csv` | verified |

## Open mapping questions (author confirmation)

Resolved:
- **Fig 5 = prostate spatial (Xenium Prime)** — confirmed.
- **Canonical Fig 2 pancreas source** — confirmed as the *updated* pancreas dataset
  version (public accession GSE148073; the committed adjacency matrices under
  `data/` derive from it), not the older public-repo notebooks.
- Liver competing methods now run through the unified 4-dataset pipeline
  (`analysis/competing_methods/competing_methods.R` + `run_all_competing_metrics.R`); the standalone
  liver Rmd/notebook is superseded and not shipped in this export.
- `supp_figS10_threshold_sensitivity.py` renamed to `analysis/supplementary/supp_figS10_threshold_sensitivity.py`
  (more descriptive); `Makefile` `extra` target and all docs updated to match.

Resolved (Phase 2.5):
- **`Supp_Treg`**, **`Supp_pca_immune_fraction`** — author-confirmed both originate from
  `analysis/figure4_5_spatial_crc_prostate/` (R spatial pipeline). Exact
  script-to-panel call may be verified later if needed; not a blocker.
- **Supp. Tables 1–3 and 7** — confirmed manually curated `.xlsx` (no producing script
  required); verified, not a blocker.

Resolved (Phase 3):
- **Fig 4/5 xGATE-on-spatial step** — the per-subdomain xGATE pathway CSVs
  (`colon_pathway_analysis_results_xgate_subdomain.csv`, `prostate_pathway_analysis_results.csv`)
  are now produced by `analysis/figure4_5_spatial_crc_prostate/run_xgate_spatial.py`
  (builds a co-expression graph from each subdomain's SCT matrix, then scores KEGG pathways
  with the xGATE engine). Previously generated off-repo.
- **CRC raw → RDS preprocessing** — documented: the Visium HD bin-to-cell step follows the
  external `bin2cell` tutorial linked in
  `analysis/figure4_5_spatial_crc_prostate/analysis/crc/scripts/Python_crc.txt`.
- **Spatial single-cell references** — added to `docs/data_availability.md` and
  `data/README.md`: CRC IRIS reference = Guimarães et al. 2024 (CELLxGENE); prostate
  SingleR/IRIS reference = Prostate Cell Atlas, Tuong et al. 2021.

Resolved (Phase 5):
- **Supp. Tables 4–6** — author-confirmed manually curated (built by hand); no producing
  script required. Verified, not a blocker.

Still open:
1. **`Supp_Senescence_reclustering`** — exact producing cell/script to be confirmed (do not guess).
   Author-acknowledged; to be resolved in a later pass.
