#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
PP_BIN="$ROOT/vendor/pp/bin/pp"
OUT=${1:-"$ROOT/build"}
PYTHON_BIN=$(command -v python3)

if [ ! -x "$PP_BIN" ]; then
  echo "pp is not built; run scripts/bootstrap_pp.sh first" >&2
  exit 1
fi

mkdir -p "$OUT"
"$PP_BIN" \
  --grant "fs:$ROOT/src:ro" \
  --grant "fs:$ROOT/data:ro" \
  --grant "fs:$ROOT/studies:ro" \
  --grant "fs:$OUT:wo" \
  --grant process \
  --schedule parallel:4 \
  --reconcile "$OUT" \
  "$ROOT/ppath.pp" -- "$ROOT" "$PYTHON_BIN"

