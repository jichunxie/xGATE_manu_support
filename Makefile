# Reproduce the xGATE manuscript analyses (Fig 3, Table 1, and supplementary panels).
#   make <target>   (run from the repository root)
# Python targets run under the `xgate` conda env via the path-setting wrapper
# analysis/shared/xgate_py.sh. See docs/reproducibility.md for setup + data layout.
# Backend script names map to manuscript analyses in docs/script_name_map.md.

PY  := ./analysis/shared/xgate_py.sh
F3  := analysis/figure3_benchmark_senescence_aab
SUP := analysis/supplementary
T1  := analysis/table1_complexity

.PHONY: all benchmark fig3 fdr nulls complexity partial extra batch help

help:
	@echo "targets: benchmark fig3 fdr nulls complexity partial extra batch all"
	@echo "  (python via $(PY); R/spatial/biology pipelines: see docs/reproducibility.md)"

# scGSEA-AUCPR-fixed benchmark metrics + Fig3 panel b + supp_figS6_9_benchmark_* PDFs
benchmark:
	$(PY) $(F3)/fig3a_assemble_benchmark.py

fig3: benchmark
	$(PY) $(F3)/fig3b_metrics_panel.py

# empirical-FDR / BH panels + jaccard vs p-value (Fig 3c,d)
fdr:
	$(PY) $(F3)/fig3d_fdr_4datasets.py
	$(PY) $(F3)/fig3c_jaccard_pvalue.py

# default vs global-degree-matched null p-value comparison (Supp_matched_null_pvalues)
nulls:
	$(PY) $(SUP)/supp_figS23_matched_null_pvalues.py

# Table 1 runtime / complexity panel
complexity:
	$(PY) $(T1)/table1_runtime.py

# partial-activity (gene-mix) panels
partial:
	$(PY) $(SUP)/supp_figS1_2_partial_activity.py
	$(PY) $(SUP)/supp_figS1_2_partial_activity.py oxphos

# threshold sensitivity + extra benchmark plots
extra:
	$(PY) $(SUP)/supp_figS10_threshold_sensitivity.py

# batch-effect figures (heavier; upstream embeddings must exist in results/)
batch:
	$(PY) $(SUP)/supp_figS3_5_batch_figures.py

all: benchmark fig3 fdr nulls complexity partial extra
	@echo "analyses regenerated. R/spatial/biology pipelines: see docs/reproducibility.md"
