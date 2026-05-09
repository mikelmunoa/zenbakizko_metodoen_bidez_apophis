#!/usr/bin/env bash
# Download SPICE kernels required by GravitationSimulation.jl
# Run from the repository root: bash data/download_kernels.sh

set -euo pipefail

KERNELS_DIR="$(dirname "$0")/kernels"
mkdir -p "$KERNELS_DIR"

NAIF_BASE="https://naif.jpl.nasa.gov/pub/naif/generic_kernels"

echo "Downloading SPICE kernels to $KERNELS_DIR ..."

# Leap seconds kernel (~50 KB)
wget -N -q --show-progress -P "$KERNELS_DIR" \
    "$NAIF_BASE/lsk/naif0012.tls"

# Planetary constants kernel (~120 KB)
wget -N -q --show-progress -P "$KERNELS_DIR" \
    "$NAIF_BASE/pck/pck00011.tpc"

# DE440 planetary ephemeris (~114 MB — planets 1900-2100)
wget -N -q --show-progress -P "$KERNELS_DIR" \
    "$NAIF_BASE/spk/planets/de440s.bsp"

# Apophis kernel (obtain from JPL Horizons or NASA CNEOS)
# https://ssd.jpl.nasa.gov/ftp/ssd/small_body/spk/asteroid/
echo ""
echo "NOTE: Apophis kernel (2099942_*.bsp) must be downloaded manually from:"
echo "  https://ssd.jpl.nasa.gov/ftp/ssd/small_body/spk/asteroid/"
echo "Place it in: $KERNELS_DIR/"

echo ""
echo "Done. Kernels downloaded to: $KERNELS_DIR"
