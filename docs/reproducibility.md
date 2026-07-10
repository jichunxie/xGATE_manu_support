# Reproducibility guide

Command-oriented guide to reproducing the xGATE manuscript analyses. Read
`docs/data_availability.md` first — large inputs must be downloaded from the
public sources listed there.

> **Two repositories.** The **method/package** is the separate
> [`xGATE`](https://github.com/jichunxie/xGATE) repo. This
> **`xGATE_manu_support`** repo only *reproduces the manuscript results* and
> depends on that package (see §1). Figure-generation and benchmark pipelines
> live here; reusable method code lives in `xGATE`.

The Python figure/table scripts resolve their paths from `configs/paths.yaml`
(via `analysis/shared/xgate_paths.py`) or fall back to the repository root, so no
absolute paths are hard-coded. Biology notebooks still use an editable placeholder
root `/path/to/xGATE` (see the note at the top of each notebook and
`configs/README.md`).

---

## 1. Environment setup

Create local environments from the provided dependency files. The original
authors used HPC environments during development; public users should create
their own Python and R environments from these files.

```bash
conda env create -f environment.yml
# or, for exact manuscript snapshots:
conda env create -f envs/xgate.yml -n xgate-manu-python
# R/spatial notes:
# see docs/r_environment.md and envs/xgate_r.yml
```

Python requirements are also listed in `requirements.txt`. R package versions
and installation notes are in `docs/r_environment.md`.
The `Bio` Python module used by KEGG pathway helpers is provided by Biopython,
listed as `biopython>=1.83` in both Python dependency files.

Python `.py` analyses launch through the wrapper, which puts the xGATE package on
`PYTHONPATH` and fixes library ordering when needed. If your Python environment
is already active, run Makefile targets directly. Advanced users can set
`XGATE_CONDA_PREFIX` to their own environment prefix before using the wrapper.

```bash
./analysis/shared/xgate_py.sh analysis/shared/<script>.py
```

**Providing the xGATE method package.** Preferred — install the pinned release
from the method repo:

```bash
pip install "git+https://github.com/jichunxie/xGATE.git@v1.0"
```

This is already listed in `requirements.txt` / `environment.yml`, so a normal
environment build installs it. The code imports it as `from xGATE.utilities import
...`. For an archival-exact build, pin the release commit hash instead of `@v1.0`.

---

## 2. Configure paths and data

```bash
cp configs/paths.example.yaml    configs/paths.yaml
cp configs/datasets.example.yaml configs/datasets.yaml
# edit both for your machine
```

Place downloaded inputs under `data/` following `data/README.md`.

---

## 3. Sanity checks (no large data required)

```bash
# imports, config load, notebook path hygiene, input-schema helper
python tests/test_sanity.py
# optional wrapper form, after activating your Python environment:
./analysis/shared/xgate_py.sh tests/test_sanity.py
```

Expected in the documented full Python environment: all checks print `PASS` and
the script exits 0. In a minimal environment, the xGATE import check may print
`SKIP` for missing third-party dependencies such as `Bio`; this means the
repository hygiene checks passed, but the full Python dependency set has not
been installed.

`python tests/test_sanity.py` is the no-data check; it needs no downloaded inputs.

Table 1 is the smallest analysis to reproduce, but its inputs
(`computational_benchmark_summary.csv` plus `computational_benchmark_summary_additions.csv`
— measured runtime/memory across 162 subpopulations / 7 studies) are **not shipped**.
Regenerate them from the graph-construction timing pipeline (`analysis/data_prep/`) or
request them from the authors, place them at the repository root, then run:

```bash
make complexity
# equivalent figure-guide wrapper:
bash figure_guide/Table1/run.sh
```

Expected outputs: `results/table1_runtime_table.csv` and
`figures/table1_complexity_analysis.{png,pdf}` (generated, git-ignored). Without the
input CSVs at the repo root, `make complexity` exits with a file-not-found error.

Reference Table 1 values (once the input is present):

- CPU vs `|V|^2`: rho `+0.954`, p `1.9e-85`
- CPU vs edge density: rho `+0.083`, p `0.29`
- Peak memory vs `|V|^2`: rho `+0.917`, p `7.1e-66`

---

## 4. Reproduce main figures

Run from the repository root. Panels land in `figures/`.

```bash
make benchmark     # fig3_benchmark_metrics_bh.csv + supp_figS6_9_benchmark_* + Fig 3 panel b input
make fig3          # Fig 3a,b metrics panel      -> figures/fig3b_metrics_panel.pdf
make fdr           # Fig 3c Jaccard + Fig 3d FDR  -> fig3c_jaccard_pvalue.pdf, fig3d_fdr_4datasets.pdf
make complexity    # Table 1 / complexity panel   -> table1_complexity_analysis.pdf
make all           # everything above + extras
```

| Figure | Command | Output panel |
|--------|---------|--------------|
| Fig 3a,b | `make fig3` | `figures/fig3b_metrics_panel.pdf` |
| Fig 3c | `make fdr` | `figures/fig3c_jaccard_pvalue.pdf` |
| Fig 3d | `make fdr` | `figures/fig3d_fdr_4datasets.pdf` |
| Table 1 | `make complexity` | `results/table1_runtime_table.csv` |

**Fig 3e (senescence)**, **Fig 2 (pancreas)**, **Fig 4 (CRC)**, **Fig 5 (prostate)**
come from the biology notebooks / R pipelines (see §6).

---

## 5. Reproduce key supplementary analyses

```bash
make nulls         # Supp_matched_null_pvalues (default vs degree-matched null)
make partial       # Supp_partial_activity (+ oxphos variant)
make extra         # Supp_threshold_sensitivity + extra benchmark plots
make batch         # Supp_batch_metrics / _pancreas / _ts   (needs upstream embeddings in results/)
```

`analysis/supplementary/supp_figS3_5_batch_figures.py` is the canonical batch-figure script (`make batch`).
Older batch variants were retired to `archive/old_numbering_or_obsolete/batch_superseded/`.

See `docs/paper_code_map.md` for the full Supp-figure → script table.

### Data preparation / upstream builds

The benchmark and graph inputs consumed by the targets above are produced by public
data-prep / graph-build scripts in `analysis/data_prep/` (`extract_*` and `*_export.R`,
`build_from_csv.py`, `build_kegg_panel.py`, `run_xgate.py`, `ts_build.py`,
`panc_donor_graph.py`, `fucci_build.py`) and the competing-method scripts in
`analysis/competing_methods/` (`prep_competing_*.py`, `competing_methods.R`,
`run_all_competing_metrics.R`). These form the upstream pipeline; run them to regenerate the
CSV/embedding inputs before the figure targets. Each script's manuscript meaning is listed in
`docs/script_name_map.md`.

---

## 6. Biology figures (notebooks + spatial R pipelines)

**Notebooks** (Jupyter in your local Python notebook environment):

```bash
jupyter lab
# Fig 2  : analysis/figure2_pancreas/*.ipynb   (pancreas = GSE148073)
# Fig 3e : analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb
# Liver benchmark competing methods: analysis/competing_methods/competing_methods.R + run_all_competing_metrics.R
```

Notebooks resolve inputs from the placeholder data root `/path/to/xGATE/data`;
point that at your `data/` directory (or edit the first cell).

> **Fig 2 notes.** The canonical pancreas source for the revised Figure 2 is the
> *updated* pancreas dataset (public accession **GSE148073**); the committed
> adjacency matrices derive from it. **Figure 2b** (marginal-expression boxplot)
> is produced by the R/Seurat notebook
> `analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` (rendered copy
> `fig2b_Mean_Exp.html`). It loads the GSE148073 Seurat object
> (`data/pancreas_human.rds`), selects control `type B pancreatic cell`s,
> normalizes with `SCTransform(vst.flavor="v2")`, pulls pathway genes from KEGG via
> `KEGGREST`/`org.Hs.eg.db`, and draws one boxplot per pathway — the same pathways
> as Figure 2a (`AMPK signaling pathway`, `Bacterial invasion of epithelial cells`).
> (An earlier standalone Python re-code was not retained; this notebook supersedes it.)

Figure 2b command:

```bash
Rscript -e 'rmarkdown::render("analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd")'
```

Local validation used the same command shape with `--input-dir
<local-pancreas-data-dir>` and wrote
`figures/figure2b_marginal_expression_recode.pdf` plus
`figures/figure2b_marginal_expression_recode.png`. These generated
outputs are git-ignored by policy; regenerate them after placing the public
GSE148073-derived matrix locally.

**Spatial pipelines** (R/Seurat environment), run the numbered scripts in order:

```bash
cd analysis/figure4_5_spatial_crc_prostate/analysis/crc/scripts
Rscript 00_config.R && Rscript 01_identify_spatial_domains.R && ...   # Fig 4
# prostate/: 00_config_prostate.R … 10_visualize_celltype_enrichment.R  # Fig 5
```

The xGATE pathway-activity CSV each pipeline reads as `PATHWAY_FILE`
(CRC: `colon_pathway_analysis_results_xgate_subdomain.csv`;
prostate: `prostate_pathway_analysis_results.csv`) is produced by the Python
scorer `analysis/figure4_5_spatial_crc_prostate/run_xgate_spatial.py`, which builds
a co-expression graph from each subdomain's SCT matrix and scores KEGG pathways with
the xGATE engine:

```bash
./analysis/shared/xgate_py.sh analysis/figure4_5_spatial_crc_prostate/run_xgate_spatial.py \
  --expr-dir <dir of per-subdomain genes-by-cells SCT CSVs> \
  --pathways-file <KEGG pathway names, one per line> \
  --out data/raw/colon_pathway_analysis_results_xgate_subdomain.csv
```

CRC Visium HD bin-to-cell preprocessing (raw 10x → `sc_input_crc.RDS` /
`countList_spatial_crc.RDS`) follows the bin2cell tutorial linked in
`.../crc/scripts/Python_crc.txt`; the prostate `01_process_raw_data.R` builds its
`sp_input` directly from the 10x raw files.

---

## 7. Local machine vs HPC

- The Python figure/analysis targets are CPU-only and run on a laptop given the
  inputs.
- Building co-expression graphs from raw counts and the batch/donor analyses are
  memory-heavy (tens of GB) — use a workstation or cluster node. The original
  HPC SLURM submission scripts are preserved (sanitized) under
  `archive/needs_review/hpc_jobs/` for reference; the `make` targets do not
  depend on a scheduler.

---

## 8. Limitations

- Large single-cell and spatial inputs are not committed; download them from the
  public sources listed in `docs/data_availability.md`.
- Composite `Figure*.jpeg` / stitched `Supp_*.pdf` assembly from the panels is a
  manual layout step, outside version control.
- Expected exact numbers: scGSEA AUCPR = Liver 0.868 / Pancreas 0.736 / FUCCI
  0.893 / TS 0.977 (benchmark), reproduced by `make benchmark`.
