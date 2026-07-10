# Figure 3 - benchmark, FDR control, senescence time-course

Method benchmark against five competitors (a,b), null-dependence / FDR control
(c,d), and a senescence time-course validation (e).

Panels a-d are Python scripts; panel e is a Jupyter notebook. Panels are written to `figures/` and assembled
into `Figure3.jpeg` on a laptop (manual layout, outside version control).

## Manuscript items

| Panel | Backend script/notebook | Required inputs | Expected output | Command | Status |
|-------|-------------------------|-----------------|-----------------|---------|--------|
| 3a,b PR / specificity / MCC / AUCPR | `analysis/figure3_benchmark_senescence_aab/fig3a_assemble_benchmark.py` → `fig3b_metrics_panel.py` | `results/fig3_benchmark_metrics_bh.csv` (from per-method `*_percall.csv`) | `figures/fig3b_metrics_panel.pdf` | `make fig3` | verified |
| 3c null p vs Jaccard overlap | `analysis/figure3_benchmark_senescence_aab/fig3c_jaccard_pvalue.py` | per-dataset inactive-pathway p-values | `figures/fig3c_jaccard_pvalue.pdf` | `make fdr` | verified |
| 3d realized FDR (BH/BY) | `analysis/figure3_benchmark_senescence_aab/fig3d_fdr_4datasets.py` (`fig3d_fdr_realized.py`) | `fig3_benchmark_metrics_bh.csv` + labels | `figures/fig3d_fdr_4datasets.pdf` | `make fdr` | verified |
| 3e senescence trends vs markers | `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb` | `Senescence Study/Results/final_analysis.csv` | trend line plot for `Figure3.jpeg` | notebook | verified; external data required |

## Run

```bash
# a-d, after activating your local Python environment
make fig3 fdr
# e, in your local Jupyter environment
jupyter lab   # open senescence_visualization.ipynb
```

## Makefile targets

- `make benchmark`: builds benchmark metrics and supplementary benchmark panels.
- `make fig3`: builds Figure 3a,b after `benchmark`.
- `make fdr`: builds Figure 3c,d.

## Notes

- `analysis/shared/` is the shared backend. Do not split these scripts by
  figure; several outputs feed main and supplementary panels.
- Final `Figure3.jpeg` assembly from PDF panels is a manual layout step outside
  version control.

## Datasets

Four benchmark datasets: liver hepatocytes (GSM6058681), pancreatic β-cells
(GSE148073), FUCCI U2OS (GSE146773), Tabula Sapiens fibroblasts. scGSEA-AUCPR
NA-guard fix baked into `fig3a_assemble_benchmark.py`; xGATE calls are BH `q<0.05`.
Reference: scGSEA AUCPR = Liver 0.868 / Pancreas 0.736 / FUCCI 0.893 / TS 0.977.
