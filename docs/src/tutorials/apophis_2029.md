# Tutorial: Apophis 2029 flyby

This tutorial shows how to propagate Apophis from early 2029 through its
April 13 close approach using `f_master!` with full physics.

## Setup

```julia
using GravitationSimulation
using LittleEphemeris
using StaticArrays
using LinearAlgebra

# Load SPICE kernels (adjust paths as needed)
furnsh("data/kernels/naif0012.tls")
furnsh("data/kernels/de440.bsp")
furnsh("data/kernels/2099942_2025_2031.bsp")
```

## Initial conditions from ephemeris

```julia
DAY = 86400.0
t0  = utc2et("2029-01-01T00:00:00")  # seconds past J2000
T   = 120 * DAY                       # 120 days

# Get Apophis state from SPICE
u0_arr, _ = spkssb(2099942, t0, "J2000")
u0 = SVector{6, Float64}(u0_arr...)
```

## Configure physics

```julia
# Full physics: GR + Earth harmonics (no Sun J2 — very small)
physics = (gr=true, j2_sun=false, j2_earth=true, j3_earth=true, j4_earth=true)

# Body indices and gravitational parameters
# (see notebooks/ for full setup)
p = build_params(t0, physics)
```

## Integrate

```julia
h   = 60.0     # 1-minute step (seconds)
sol = RK4(u0, t0, t0 + T, h, p, f_master!)

# Find minimum Earth distance
earth_dist = [norm(sol[i][1:3] .- Earth(t0 + (i-1)*h)[1:3]) for i in eachindex(sol)]
min_dist, idx = findmin(earth_dist)
t_closest = t0 + (idx - 1) * h

println("Minimum distance: $(round(min_dist, digits=0)) km")
println("Closest approach: $(et2utc(t_closest, "ISOC", 0))")
```

## Effect of post-Newtonian correction

To quantify the GR effect, compare trajectories with and without:

```julia
p_nogr = (p..., physics=(gr=false, j2_sun=false, j2_earth=true, j3_earth=true, j4_earth=true))

sol_gr   = RK4(u0, t0, t0 + T, h, p,      f_master!)
sol_nogr = RK4(u0, t0, t0 + T, h, p_nogr, f_master!)

# Position difference at closest approach
Δr = norm(sol_gr[idx][1:3] .- sol_nogr[idx][1:3])
println("GR position shift at flyby: $(round(Δr, digits=1)) km")
# Expected: ~6 km
```
