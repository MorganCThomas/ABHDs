#!/usr/bin/env bash
# Vendor the pre-built Mol* viewer (pinned version) so protein reports can show 3D structures offline.
# Source: jsDelivr mirror of the npm package `molstar` (MIT licence, https://molstar.org). ~5.3 MB JS + 73 KB CSS into data/vendor/ (git-ignored).
# Without it the reports load the same version from jsDelivr (internet needed).
set -euo pipefail
V=5.13.1
D="$(git rev-parse --show-toplevel)/data/vendor/molstar-$V"
mkdir -p "$D"
for f in molstar.js molstar.css; do curl -fsSL -o "$D/$f" "https://cdn.jsdelivr.net/npm/molstar@$V/build/viewer/$f"; done
ls -la "$D"
