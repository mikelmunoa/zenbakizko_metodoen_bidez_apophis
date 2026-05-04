# GravitationSimulation.jl

High-precision N-body gravitational simulation of near-Earth asteroids in Julia.

## Overview

This package implements the gravitational equations of motion for asteroid trajectory
propagation, including:

- Newtonian N-body gravity using SPICE-based planetary ephemerides
- Post-Newtonian (Damour-Deruelle 1985) general relativistic corrections
- Zonal harmonic perturbations: Earth J₂/J₃/J₄, Sun J₂
- Configurable physics via [`f_master!`](@ref)

The primary application is the trajectory analysis of asteroid (99942) Apophis
during its April 2029 close approach to Earth (~38,000 km minimum distance).

## Package structure

```
src/
├── GravitationSimulation.jl   # Module definition and exports
├── equations.jl               # Physics functions (f_master!, dd_perturbation, ...)
└── integrator.jl              # RK4 integrator
```

## Quick example

```julia
using GravitationSimulation, StaticArrays

# Compute D-D perturbation at 1 AU
AU = 1.495978707e8
u   = SVector(AU, 0.0, 0.0, 0.0, 30.0, 0.0)
sun = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

dd = dd_perturbation(u, sun, 1.32712440018e11)
```

## Contents

```@contents
Pages = ["physics.md", "api.md", "tutorials/apophis_2029.md"]
Depth = 2
```
