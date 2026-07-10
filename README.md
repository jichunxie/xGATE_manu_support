# xGATE — manuscript code & reproducibility

Analysis code, figure/table generation, and reproducibility documentation
accompanying the xGATE manuscript. It covers single-cell pathway-activity
analysis, the four-dataset method benchmark, false-discovery-rate and null
analyses, batch/donor robustness, computational-complexity profiling, a
senescence time-course, and spatial transcriptomics of human cancers.

> **Paper:** _xGATE: Using gene co-expression graph topological fingerprints to enhance pathway activity scoring._

---

## Two repositories

The project is split into two GitHub repositories:

| Repository | Role | Contents |
|------------|------|----------|
| [**`xGATE`**](https://github.com/jichunxie/xGATE) | **Method / package** | The reusable xGATE method implementation, package code, core API, install instructions, and lightweight usage examples. Not intended to reproduce every manuscript figure. |
| **`xGATE_manu_support`** (this repo) | **Manuscript reproducibility** | Code to reproduce the paper's figures, supplementary figures, tables, benchmarks, and manuscript-specific analyses, plus the mapping from each manuscript item to its script/inputs/outputs. |

**This repo depends on the `xGATE` package rather than re-developing it.** Reusable
method code belongs in `xGATE`; manuscript-specific figure/benchmark/analysis
scripts belong here. The method is installed as an **external pinned dependency**
(`pip install git+https://github.com/jichunxie/xGATE.git@v1.0`) and imported as
`from xGATE.utilities import ...` — see
[How this repo depends on xGATE](#how-this-repo-depends-on-xgate).

---

## Method summary (brief)

The full method is documented in the [`xGATE`](https://github.com/jichunxie/xGATE)
repository. In short:

xGATE takes (i) pathways of interest from a knowledge base (e.g. KEGG) and (ii) a
single-cell RNA-seq dataset for a cell population, builds a gene co-expression
graph per pathway, and uses a graph-based deep-learning model to capture the
graph's topological structure. It outputs a p-value and effect size for the
*activity* of each pathway — distinguishing coordinated (active) from
non-coordinated (inactive) pathways, which marginal mean-expression methods
cannot separate. Pathway calls are Benjamini–Hochberg FDR-corrected per dataset.

> **Naming note.** This repository uses the manuscript name **xGATE** throughout.
> Any historical internal names in older intermediate files have been renamed or
> documented only where needed for compatibility.

---

## Repository structure

```text
xGATE_manu_support/
├── README.md
├── LICENSE.md                    # MIT License
├── Makefile                      # figure/table reproduction targets
├── environment.yml               # Python environment
├── requirements.txt              # Python dependencies (+ how to install the xGATE method)
├── CITATION.cff                  # how to cite this repository
├── CHANGELOG.md
├── configs/                      # path and dataset config templates
├── data/                         # input layout and public download instructions
├── docs/                         # reproducibility, paper-code map, data availability
├── envs/                         # optional environment snapshots
├── figure_guide/                 # figure/table-level entry points
├── analysis/                     # manuscript analysis code
│   ├── shared/                   #   path resolver, launcher, shared helpers
│   ├── data_prep/                #   dataset extraction, normalization, graph construction
│   ├── competing_methods/        #   ORA/AUCell/scGSEA/GESECA/PAGODA benchmark + prep
│   ├── figure2_pancreas/         #   Fig 2 pancreas co-expression graphs, donor sensitivity, and embeddings
│   ├── figure3_benchmark_senescence_aab/   # Fig 3 benchmark/FDR/AAB (+ senescence/ subfolder = Fig 3e)
│   ├── figure4_5_spatial_crc_prostate/     # Fig 4 CRC and Fig 5 prostate spatial analyses
│   ├── table1_complexity/        #   Table 1 complexity/runtime analysis
│   └── supplementary/            #   supplementary batch/null/partial-activity analyses
└── tests/                        # lightweight sanity checks
```

Generated outputs are written locally to `figures/` and `results/`; these
directories are ignored by Git and are not part of the tracked repository.

Start with [`figure_guide/`](figure_guide/README.md) to reproduce a specific figure or
table. The [`analysis/`](analysis/) directory holds the corresponding backend scripts.
The **xGATE method** package (canonical repo: <https://github.com/jichunxie/xGATE>) is
installed as a pinned external dependency (see below). Backend script names map to
manuscript analyses in [`docs/script_name_map.md`](docs/script_name_map.md).

---

## Installation

```bash
# 1. clone, then create local environments from the provided dependency files
conda env create -f environment.yml                 # Python analyses + notebooks
# For exact manuscript snapshots, see envs/ and docs/r_environment.md.

# 2. configure paths + datasets for your machine
cp configs/paths.example.yaml    configs/paths.yaml
cp configs/datasets.example.yaml configs/datasets.yaml   # then edit both
```

Full Python/R notes: [`docs/reproducibility.md`](docs/reproducibility.md) and
[`docs/r_environment.md`](docs/r_environment.md). The original authors used HPC
environments during development; public users should create their own local
environment from the listed dependency files.

### How this repo depends on xGATE

The manuscript analyses need the **xGATE method package** (a separate repo),
installed as a **pinned external dependency**:

```bash
pip install "git+https://github.com/jichunxie/xGATE.git@v1.0"
```

This is included in `requirements.txt` / `environment.yml`, so a normal environment
build installs it automatically. Scripts and notebooks then import it as
`from xGATE.utilities import ...` / `import xGATE.utilities`; the launcher
`analysis/shared/xgate_py.sh` only adds `analysis/shared` to `PYTHONPATH` (the
method itself comes from the installed package). For an archival-exact build, pin
the release commit hash instead of the `@v1.0` tag. See
[`docs/xgate_dependency_plan.md`](docs/xgate_dependency_plan.md).

The Python figure/table scripts read their paths from `configs/paths.yaml` (via
`analysis/shared/xgate_paths.py`) or fall back to the repo root — no absolute
paths are hard-coded. Biology notebooks still use an editable placeholder root
(see the note at the top of each notebook).

---

## Data availability

Input data are **not** committed because the raw single-cell and spatial matrices
are large. Download the public datasets and place them under `data/` — see
[`data/README.md`](data/README.md) and
[`docs/data_availability.md`](docs/data_availability.md). All dataset accessions
currently listed there are public.

---

## Reproduce the main figures

From the repository root, after activating your local Python environment and
placing input data under `data/`:

```bash
make fig3          # Fig 3a,b  benchmark precision-recall / MCC / AUCPR
make fdr           # Fig 3c,d  Jaccard vs null p-value + realized FDR (BH/BY)
make complexity    # Table 1   computational complexity + edge density
make all           # all Python figure/table targets + extras
```

Biology figures (Jupyter / R environments; Figure 2b is an R/Seurat notebook):
- **Fig 2** — `analysis/figure2_pancreas/*.ipynb`; Fig 2b:
  `analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` (render to HTML locally)
- **Fig 3e** — `analysis/figure3_benchmark_senescence_aab/senescence/fig3e_senescence_visualization.ipynb`
- **Fig 4 (CRC)** / **Fig 5 (prostate)** — numbered R pipelines under
  `analysis/figure4_5_spatial_crc_prostate/analysis/{crc,prostate}/scripts/`

Full command list and figure↔script table:
[`docs/reproducibility.md`](docs/reproducibility.md),
[`docs/paper_code_map.md`](docs/paper_code_map.md).

## Reproduce supplementary analyses

```bash
make nulls         # Supp_matched_null_pvalues
make partial       # Supp_partial_activity (+ oxphos)
make extra         # Supp_threshold_sensitivity + extra benchmark plots
make batch         # Supp_batch_* (needs upstream embeddings)
```

---

## Expected outputs

- Metric CSVs and panels under `results/` and `figures/`.
- Reference values: scGSEA AUCPR = Liver 0.868 / Pancreas 0.736 / FUCCI 0.893 /
  TS 0.977 (from `make benchmark`).
- Composite `Figure*.jpeg` / `Supp_*.pdf` are assembled from these panels by hand
  (layout step, outside version control).

## Sanity checks

```bash
python tests/test_sanity.py     # imports, config parse, privacy scan, input-schema doc
bash figure_guide/Table1/run.sh # lightweight Table 1 smoke test
```

---

## Citation

> _Using gene co-expression graph topological fingerprints to enhance pathway
> activity scoring._

Citation will be added upon publication / preprint release.

## Contact

- For questions about the manuscript support code, please contact **Yuxia Xie**
  (<yuxia.xie@duke.edu>).
- PI: **Jichun Xie** (<jichun.xie@duke.edu>).
- Method / package: <https://github.com/jichunxie/xGATE>.

## Release status

This repository accompanies the xGATE manuscript, which is not yet published; the
citation DOI will be finalized on publication. The analyses depend on the standalone
[xGATE](https://github.com/jichunxie/xGATE) package, pinned to release `v1.0` — see
[`docs/xgate_dependency_plan.md`](docs/xgate_dependency_plan.md).

## License

This manuscript-support repository is released under the **MIT License**. Please
see [`LICENSE.md`](LICENSE.md) for permitted use.
