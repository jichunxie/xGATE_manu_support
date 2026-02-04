# xGATE_manu_support

Repository overview
-------------------
This repository contains analysis scripts and Jupyter/R notebooks for single-cell data analysis and spatial cancer analyses (folders include `Human_Cancer_Spatial_Analysis`, `Pancreas`, `Liver`, etc.). The project contains both R and Python components — installation instructions are provided for each and a top-level installer is available.

Installation (recommended)
-------------------------
R (recommended central installer):

1. From the repository root, run the top-level R installer which will install core R packages, source the R analysis install script, and install `xGATE` from GitHub:

```r
# from the repository root in R
source("install_dependencies.R")
```

2. The top-level installer calls the existing R install script at [Human_Cancer_Spatial_Analysis/install_dependencies.R](Human_Cancer_Spatial_Analysis/install_dependencies.R) to install R dependencies.

Note about `xGATE`:

`xGATE` used by the notebooks in this repository is a Python package available at https://github.com/jichunxie/xGATE. It is installed via `pip` from the repository (see `requirements.txt`).

Python install (xGATE):

```bash
# from the repository root
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

This installs `xGATE` from the GitHub repository (via the `git+https://...` entry in `requirements.txt`).

Manual pip install (if you prefer):

```bash
pip install git+https://github.com/jichunxie/xGATE.git
```

Verify Python `xGATE` installation:

Create a short Python script or run in Python REPL:

```python
try:
    import xgate as x
    print("Imported xgate as 'xgate' module")
except Exception:
    try:
        import xGATE as x
        print("Imported xGATE as 'xGATE' module")
    except Exception:
        raise SystemExit("xGATE import failed — check pip install logs and requirements.txt")
```

Python (notebooks)
------------------
Notebooks in the repo require a Python environment. A minimal `requirements.txt` is provided; create a virtual environment and install with:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Edit `requirements.txt` to add domain-specific packages used in particular notebooks.

Notes
-----
- The install scripts assume an internet connection and appropriate permissions to install system packages.
- On macOS, ensure Xcode command-line tools are installed: `xcode-select --install`.

Contact
-------
For issues with `xGATE` itself see: https://github.com/jichunxie/xGATE
