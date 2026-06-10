module GravitationSimulation

using LinearAlgebra, SPICE, LittleEphemeris, StaticArrays

# Fitxategiak kargatu (ls-n agertzen diren izenekin)
include("integrator.jl")
include("equations.jl")
include("system.jl")

# Notebook-ean erabili nahi dituzun funtzioak esportatu
export RK4, f_all!, f_all_rel!, f_all_rel_J!, f_all_rel_J2_Sun!, f_all_rel_J2_Sun_Earth!,
       f_master!, dd_perturbation,
       planet_system, DEFAULT_PHYSICS, PLANET_IDS, ASTEROID_IDS, ASTEROID_GMS

end # module
