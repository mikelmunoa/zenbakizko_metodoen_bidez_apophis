#!/usr/bin/env bash
# Backwards-compatible wrapper. Cross-platform logika `download_kernels.jl`-en dago.
# Sistema guztietarako (Linux / macOS / Windows) zuzenean ere erabil daiteke:
#     julia data/download_kernels.jl

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec julia "$SCRIPT_DIR/download_kernels.jl"
