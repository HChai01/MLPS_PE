#!/bin/sh
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
PYTHON_BIN=${PYTHON:-python3}
VENV_DIR="$ROOT_DIR/.mlps_pe_runtime"
REQUIREMENTS_FILE="$ROOT_DIR/python/toolkit/requirements-toolkit.txt"

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "Python is required to build MLPS_PE; set PYTHON=/path/to/python" >&2
    exit 1
fi

"$PYTHON_BIN" -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install --upgrade pip
"$VENV_DIR/bin/python" -m pip install -r "$REQUIREMENTS_FILE"

cat <<EOF
Built MLPS_PE for $(uname -s)-$(uname -m):
  $ROOT_DIR/mlps_pe

The executable uses the bundled runtime at:
  $VENV_DIR
EOF
