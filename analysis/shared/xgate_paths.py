"""Central path resolution for the xGATE manuscript-support analyses.

Removes hard-coded absolute paths from the analysis scripts. Resolution order
for the repository root:

  1. ``$XGATE_ROOT``            — exported by ``xgate_py.sh`` (preferred).
  2. ``configs/paths.yaml``     — ``repo_root:`` value, if the file exists and
                                  PyYAML is installed.
  3. auto-detected             — two directories up from this file
                                  (``analysis/shared/`` -> repo root).

Optional per-root overrides from ``configs/paths.yaml`` (``data_root``,
``results_root``, ``figures_root``) are honored when present.

Usage in a script::

    from xgate_paths import ROOT, RESULTS, FIGURES, DATA
    df = pd.read_csv(RESULTS + "/fig3_benchmark_metrics_bh.csv")
"""
from __future__ import annotations

import os
from pathlib import Path


def _auto_root() -> Path:
    # this file lives at <repo>/analysis/shared/xgate_paths.py
    return Path(__file__).resolve().parents[2]


def _from_config(repo_root: Path) -> dict:
    cfg_path = repo_root / "configs" / "paths.yaml"
    if not cfg_path.is_file():
        return {}
    try:
        import yaml  # type: ignore
    except ImportError:
        return {}
    try:
        return yaml.safe_load(cfg_path.read_text()) or {}
    except Exception:  # noqa: BLE001
        return {}


def _resolve() -> dict:
    env_root = os.environ.get("XGATE_ROOT")
    root = Path(env_root).resolve() if env_root else _auto_root()

    cfg = _from_config(root)
    if not env_root and cfg.get("repo_root") and cfg["repo_root"] != "/path/to/xGATE":
        root = Path(cfg["repo_root"]).resolve()
        cfg = _from_config(root)  # re-read relative to configured root

    def pick(key: str, default: Path) -> str:
        val = cfg.get(key)
        if val and val != f"/path/to/xGATE/{default.name}" and "/path/to/" not in str(val):
            return str(val)
        return str(default)

    return {
        "ROOT": str(root),
        "DATA": pick("data_root", root / "data"),
        "RESULTS": pick("results_root", root / "results"),
        "FIGURES": pick("figures_root", root / "figures"),
    }


_p = _resolve()
ROOT: str = _p["ROOT"]
DATA: str = _p["DATA"]
RESULTS: str = _p["RESULTS"]
FIGURES: str = _p["FIGURES"]

__all__ = ["ROOT", "DATA", "RESULTS", "FIGURES"]
