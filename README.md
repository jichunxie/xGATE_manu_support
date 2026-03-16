# xGATE Manuscript Support

This repository contains the analysis scripts, Jupyter notebooks, and R Markdown documents that accompany the xGATE manuscript. It covers single-cell pathway activity analysis, spatial transcriptomics of human cancers, senescence studies, and benchmarking against competing methods.

The core method is implemented in the **xGATE** Python package:
> https://github.com/jichunxie/xGATE

---

## Table of Contents

1. [Repository Structure](#repository-structure)
2. [Prerequisites](#prerequisites)
3. [Python Setup & xGATE Installation](#python-setup--xgate-installation)
4. [R Setup](#r-setup)
5. [Notebooks & Scripts Overview](#notebooks--scripts-overview)
6. [Troubleshooting](#troubleshooting)
7. [Contact](#contact)

---

## Repository Structure

```
xGATE_manu_support/
├── requirements.txt                     # Python dependencies (includes xGATE)
├── verify_xgate.py                      # Quick install verification script
├── Pancreas/                            # Pancreas single-cell xGATE analyses
│   ├── Pancreas_xGATE_CTRL_T1D.ipynb
│   ├── AAB_xGATE.ipynb
│   ├── AMPK_vs_Bacterial_Invasion_Coexp_Networks.ipynb
│   ├── Embedding_comparison.ipynb
│   └── Competing_methods_ORA_AUCell_scGSEA.Rmd
├── Liver/                               # Liver / hepatocyte pathway analyses
│   ├── xGATE_Hepatocyte_Pathway_Analysis.ipynb
│   └── Competing_methods_ORA_AUCell_scGSEA.Rmd
├── Human_Cancer_Spatial_Analysis/       # Spatial transcriptomics (CRC, prostate)
│   ├── install_dependencies.R
│   └── R/                              # Shared R analysis functions
├── Benchmarking/                        # GESECA & Pagoda pathway benchmarks (R)
├── Senescence Study/                    # Senescence visualization notebook
└── Supp/                               # Supplementary figure notebooks
```

---

## Prerequisites

| Requirement | Version |
|---|---|
| Python | ≥ 3.8 |
| pip | ≥ 21 |
| Git | any recent version |
| R | ≥ 4.1 (for R notebooks only) |

On **macOS**, make sure Xcode command-line tools are installed before proceeding:

```bash
xcode-select --install
```

---

## Python Setup & xGATE Installation

All Python notebooks in this repository depend on **xGATE**. Follow the steps below to create an isolated environment and install every required package, including xGATE itself.

### Step 1 — Create and activate a virtual environment

```bash
# from the repository root
python3 -m venv .venv
source .venv/bin/activate   # macOS / Linux
# .venv\Scripts\activate    # Windows
```

### Step 2 — Install all dependencies (including xGATE)

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

`requirements.txt` installs xGATE directly from its GitHub repository:

```
git+https://github.com/jichunxie/xGATE.git
```

No separate installation step is needed.

> **Alternative — install xGATE only:**
> ```bash
> pip install git+https://github.com/jichunxie/xGATE.git
> ```

### Step 3 — Verify the installation

Run the bundled verification script from the repository root:

```bash
python verify_xgate.py
```

Expected output (one of the following):

```
xGATE imported as 'xgate' module
```

or

```
xGATE imported as 'xGATE' module
```

If you see an error instead, see [Troubleshooting](#troubleshooting).

### Step 4 — Launch Jupyter Notebook

```bash
jupyter notebook
```

Then open any notebook from the `Pancreas/`, `Liver/`, `Senescence Study/`, or `Supp/` folders.

> **Important — import name:** xGATE registers itself under the lower-case module name `xgate`. Use `import xgate` in your code. If that fails, try `import xGATE` (older builds). The `verify_xgate.py` script checks both automatically.

---

## R Setup

The R components (R Markdown documents and helper scripts) use CRAN, Bioconductor, and Seurat packages.

Install all R dependencies with:

```r
# run once from the repository root in an R session
source("Human_Cancer_Spatial_Analysis/install_dependencies.R")
```

This script installs:

- **CRAN:** `ggplot2`, `dplyr`, `Seurat`, `patchwork`, `ComplexHeatmap`, `viridis`, and more
- **Bioconductor:** `SingleR`, `SingleCellExperiment`, `SummarizedExperiment`, `scuttle`, `BiocParallel`

After installation, open any `.Rmd` file (e.g., in `Benchmarking/` or `Pancreas/`) in RStudio and knit or run interactively.

---

## Notebooks & Scripts Overview

### Pancreas (`Pancreas/`)

| File | Description |
|---|---|
| `Pancreas_xGATE_CTRL_T1D.ipynb` | xGATE pathway activity analysis comparing healthy controls vs. T1D |
| `AAB_xGATE.ipynb` | AAB-level xGATE analysis |
| `AAB_Figure.ipynb` | Figure generation for AAB results |
| `Figure_2d_AAB_generate.ipynb` | Generates Figure 2d |
| `AMPK_vs_Bacterial_Invasion_Coexp_Networks.ipynb` | Co-expression network comparison of AMPK and bacterial invasion pathways |
| `Embedding_comparison.ipynb` | Comparison of graph embedding methods (CPACT, FeatherGraph, Graph2Vec, NetLSD) |
| `Competing_methods_ORA_AUCell_scGSEA.Rmd` | Benchmark of ORA, AUCell, and scGSEA on pancreas data |
| `Contour plots for Pancreas and Hepatocytes.Rmd` | Contour visualizations |

### Liver (`Liver/`)

| File | Description |
|---|---|
| `xGATE_Hepatocyte_Pathway_Analysis.ipynb` | xGATE pathway analysis on hepatocyte data |
| `Competing_methods_ORA_AUCell_scGSEA.Rmd` | Benchmark of competing methods on liver data |
| `Contour plots for Pancreas and Hepatocytes.Rmd` | Shared contour plot generation |

### Human Cancer Spatial Analysis (`Human_Cancer_Spatial_Analysis/`)

Spatial transcriptomics analyses of colorectal cancer (CRC) and prostate cancer.  
R functions are organized under `R/` (data ingestion, processing, clustering, visualization, cell annotation, correlation heatmaps).

### Benchmarking (`Benchmarking/`)

| File | Description |
|---|---|
| `GESECA pathway analysis.Rmd` | GESECA pathway scoring benchmark |
| `Pagoda pathway analysis.Rmd` | Pagoda2 pathway scoring benchmark |

### Senescence Study (`Senescence Study/`)

`senescence_visualization.ipynb` — Visualization of pathway activity in a senescence dataset.

### Supplementary (`Supp/`)

`supp_visualization.ipynb` — Supplementary figure generation.

---

## Troubleshooting

### `xGATE import failed`

1. Confirm the virtual environment is active: `which python` should point to `.venv/bin/python`.
2. Confirm xGATE installed:
   ```bash
   pip show xgate
   ```
   If not found, reinstall:
   ```bash
   pip install --force-reinstall git+https://github.com/jichunxie/xGATE.git
   ```
3. If `pip install` fails with a Git error, make sure Git is installed and you have network access:
   ```bash
   git --version
   curl -I https://github.com
   ```
4. On macOS, if you get compiler errors during installation, ensure Xcode tools are installed:
   ```bash
   xcode-select --install
   ```

### Jupyter kernel does not see the virtual environment

Register the environment as a Jupyter kernel:

```bash
pip install ipykernel
python -m ipykernel install --user --name xgate-env --display-name "Python (xgate-env)"
```

Then select **Python (xgate-env)** from the kernel menu in JupyterLab.

### R package installation failures

- For Bioconductor version mismatches, update BiocManager first:
  ```r
  install.packages("BiocManager")
  BiocManager::install(version = "3.18")
  ```
- Seurat installation issues on macOS may require setting the correct compiler; see the [Seurat install guide](https://satijalab.org/seurat/articles/install.html).

---

## Contact

- **xGATE package issues:** https://github.com/jichunxie/xGATE
- **Manuscript analysis issues:** open an issue in this repository
