# GravitationSimulation.jl

[![CI](https://github.com/mikelmunoa/zenbakizko_metodoen_bidez_apophis/actions/workflows/CI.yml/badge.svg)](https://github.com/mikelmunoa/zenbakizko_metodoen_bidez_apophis/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mikelmunoa.github.io/zenbakizko_metodoen_bidez_apophis/dev)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

High-precision N-body gravitational simulation of near-Earth asteroids in Julia, with post-Newtonian corrections and zonal harmonic perturbations. Developed as a Bachelor's thesis (GrAL) at the University of the Basque Country (UPV/EHU).

## Features

- **IRKGL16** — 16th-order implicit Gauss-Legendre symplectic integrator
- **Post-Newtonian GR** — Damour-Deruelle (1985) correction; ~6 km effect on Apophis by 2029 flyby
- **Zonal harmonics** — Earth J₂/J₃/J₄ and Sun J₂ perturbations in body-pole reference frames
- **Configurable physics** via `f_master!`: enable/disable GR and harmonics independently (ASSIST-style)
- **LittleEphemeris.jl** — SPICE-based Chebyshev ephemeris interpolation with zero-allocation evaluation
- **Apophis 2029** — close approach trajectory analysis (minimum distance ~38,000 km)

## Quick Start

```julia
using Pkg
Pkg.add(url="https://github.com/mikelmunoa/zenbakizko_metodoen_bidez_apophis")

using GravitationSimulation, StaticArrays
```

### Damour-Deruelle perturbation

```julia
# Apophis state (km, km/s, barycentric ICRF)
u_apophis = SVector(x, y, z, vx, vy, vz)
sun_state = SVector(xs, ys, zs, vxs, vys, vzs)
mu_sun    = 1.32712440018e11  # km³/s²

dd = dd_perturbation(u_apophis, sun_state, mu_sun)
# → SVector{3}: acceleration in km/s²
```

### Master integrator with configurable physics

```julia
physics = (gr=true, j2_sun=false, j2_earth=true, j3_earth=true, j4_earth=false)
p = (bodies=bodies, mus=mus, idx=idx, physics=physics)

sol = RK4(u0, t0, T, h, p, f_master!)
```

## Installation & Data Setup

The package requires SPICE kernels (~500 MB) not included in the repository.
Download them with:

```bash
bash data/download_kernels.sh
```

Or manually from [NAIF/JPL](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/):
- `de440.bsp` — planetary ephemeris
- `naif0012.tls` — leap seconds kernel
- `2099942_2025_2031.bsp` — Apophis kernel

## Mathematical Background

The equations of motion are:

$$\ddot{\mathbf{r}} = \sum_{j \neq i} \frac{\mu_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3} + \mathbf{a}_\text{PN} + \mathbf{a}_{J_2} + \mathbf{a}_{J_3} + \mathbf{a}_{J_4}$$

The post-Newtonian term follows Damour & Deruelle (1985):

$$\mathbf{a}_\text{PN} = \frac{GM_\odot}{c^2 r^3}\left[\left(\frac{4GM_\odot}{r} - v^2\right)\mathbf{r} + 4(\mathbf{r}\cdot\mathbf{v})\mathbf{v}\right]$$

See [`docs/src/physics.md`](docs/src/physics.md) for the full formulation.

## Citation

If you use this software in your research, please cite using [`CITATION.cff`](CITATION.cff):

```bibtex
@software{munoa2025gravitationsimulation,
  author  = {Muñoa Illarregi, Mikel},
  title   = {{GravitationSimulation.jl}: High-precision N-body simulation for near-Earth asteroids},
  year    = {2025},
  url     = {https://github.com/mikelmunoa/zenbakizko_metodoen_bidez_apophis},
  license = {MIT}
}
```

## References

- Damour, T. & Deruelle, N. (1985). *Ann. Inst. Henri Poincaré*, 43(1), 107–132.
- Holman, M. J. et al. (2023). ASSIST. *PSJ*, 4(4). [doi:10.3847/PSJ/acc9a9](https://doi.org/10.3847/PSJ/acc9a9)
- Folkner, W. M. et al. (2021). DE440/DE441. JPL IOM 392R-21-005.

## License

MIT — see [LICENSE](LICENSE).
