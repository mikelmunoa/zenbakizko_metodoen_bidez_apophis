#!/usr/bin/env julia
# Download SPICE kernels required by GravitationSimulation.jl.
# Cross-platform (Linux / macOS / Windows): uses Julia's Downloads stdlib.
#
# Usage from the repository root:
#   julia data/download_kernels.jl

using Downloads

const KERNELS_DIR = joinpath(@__DIR__, "kernels")
const NAIF_BASE   = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"

const FILES = [
    # (URL, gutxi gorabeherako tamaina deskribapenerako)
    ("$NAIF_BASE/lsk/naif0012.tls",        "leap seconds (~50 KB)"),
    ("$NAIF_BASE/pck/pck00011.tpc",        "planetary constants (~120 KB)"),
    ("$NAIF_BASE/spk/planets/de440s.bsp",  "DE440s planetary ephemeris (~114 MB)"),
]

mkpath(KERNELS_DIR)
println("Downloading SPICE kernels to ", KERNELS_DIR, " ...")

for (url, label) in FILES
    name = basename(url)
    dest = joinpath(KERNELS_DIR, name)
    if isfile(dest)
        println("  ✓ $name dagoeneko bertan dago — saltatzen.")
        continue
    end
    println("  → $name  ($label)")
    Downloads.download(url, dest)
end

println()
println("NOTE: Apophis kernel (2099942_*.bsp) must be downloaded manually from:")
println("  https://ssd.jpl.nasa.gov/ftp/ssd/small_body/spk/asteroid/")
println("Place it in: ", KERNELS_DIR)
println()
println("Done. Kernels downloaded to: ", KERNELS_DIR)
