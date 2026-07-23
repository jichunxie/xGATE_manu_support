# xGATE dependency

How `xGATE_manu_support` (this reproducibility repo) depends on the separate
`xGATE` method package. As of the `v1.0` release the method is a **pinned external
dependency**; the earlier vendored snapshot has been removed.

## Repository roles

- **`xGATE`** — canonical method/package repository. Reusable implementation,
  core API, install instructions, lightweight examples.
- **`xGATE_manu_support`** — manuscript reproducibility repository. Figure/table/
  benchmark/analysis scripts and the paper→code mapping. Depends on `xGATE`.

## How the dependency is installed

The method is installed pinned to the manuscript release tag:

```bash
pip install "git+https://github.com/jichunxie/xGATE.git@v1.0"
```

This line is included in `requirements.txt` and `environment.yml`, so a normal
environment build installs it automatically. Both the scripts and the biology
notebooks import it uniformly as:

```python
from xGATE.utilities import ...   # scripts
import xGATE.utilities            # notebooks
```

`analysis/shared/xgate_py.sh` only adds `analysis/shared` to `PYTHONPATH`; the
method itself comes from the installed package.

For an **archival-exact** build, replace `@v1.0` with the release commit hash so the
dependency cannot drift if the tag is ever moved.

## Background: the import reconciliation (done in the xGATE repo at v1.0)

Earlier, the method's `setup.py` declared `name="xGATE"` but `find_packages()`
shipped a top-level **`utilities`** package (no `xGATE/__init__.py` wrapping it), so
after `pip install`:

- `import utilities`         → worked
- `import xGATE.utilities`   → **failed** (used by the biology notebooks)

This split forced a vendored snapshot with a two-entry `PYTHONPATH` hack. It was
fixed in the `xGATE` repo for `v1.0` by restructuring the package so the installable
layout is `xGATE/` containing `xGATE/utilities/`. A clean `pip install` now provides
`import xGATE.utilities`, so this repo standardizes on `from xGATE.utilities import
...` everywhere and the vendored snapshot is gone.
