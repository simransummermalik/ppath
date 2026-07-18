#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
export OPAMROOT="$ROOT/.opam-root"

if [ -d /opt/homebrew/opt/zlib/lib/pkgconfig ]; then
  PP_ZLIB_PKGCONFIG=/opt/homebrew/opt/zlib/lib/pkgconfig
  if [ -n "${PKG_CONFIG_PATH:-}" ]; then
    export PKG_CONFIG_PATH="$PP_ZLIB_PKGCONFIG:$PKG_CONFIG_PATH"
  else
    export PKG_CONFIG_PATH="$PP_ZLIB_PKGCONFIG"
  fi
fi

if [ ! -d "$ROOT/vendor/pp/.git" ]; then
  mkdir -p "$ROOT/vendor"
  git clone https://github.com/pranavra0/pp.git "$ROOT/vendor/pp"
fi

if [ ! -d "$ROOT/_opam" ]; then
  opam init --bare --disable-sandboxing --no-setup -y
  opam switch create "$ROOT" 5.3.0 -y
fi

opam install --switch="$ROOT" "$ROOT/vendor/pp" --deps-only -y
opam exec --switch="$ROOT" -- dune build --root="$ROOT/vendor/pp"
echo "built $ROOT/vendor/pp/bin/pp"
