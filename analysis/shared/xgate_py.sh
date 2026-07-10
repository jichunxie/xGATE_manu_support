#!/usr/bin/env bash
# Launch Python for the xGATE revision analyses with the xGATE package on PYTHONPATH.
#
# Resolves the repository root from this script's own location, so no absolute
# user/cluster paths are baked in. Point XGATE_CONDA_PREFIX at your conda env
# (the one built from ../../envs/xgate.yml); it defaults to the active env.
#
#   XGATE_CONDA_PREFIX=/path/to/conda/envs/xgate ./xgate_py.sh some_script.py

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
export XGATE_ROOT="$REPO_ROOT"

# conda env: override with XGATE_CONDA_PREFIX, else use the active env
PREFIX="${XGATE_CONDA_PREFIX:-${CONDA_PREFIX:-}}"
if [[ -z "$PREFIX" ]]; then
  echo "xgate_py.sh: set XGATE_CONDA_PREFIX to your 'xgate' conda env (see envs/xgate.yml)" >&2
  exit 1
fi

# put conda libstdc++ ahead of system libs (avoids CXXABI mismatch), and the
# shared analysis helpers (<repo>/analysis/shared) on the path. The xGATE method
# package is installed into the conda env (pip install git+...@v1.0; see
# requirements.txt / environment.yml), so scripts and notebooks import it as
# `from xGATE.utilities import ...` with no PYTHONPATH entry.
export LD_LIBRARY_PATH="$PREFIX/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$REPO_ROOT/analysis/shared:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

exec "$PREFIX/bin/python" "$@"
