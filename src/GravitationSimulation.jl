module GravitationSimulation

using LinearAlgebra, SPICE, LittleEphemeris, StaticArrays

# Fitxategiak kargatu (ls-n agertzen diren izenekin)
include("integrator.jl")
include("equations.jl")

# Notebook-ean erabili nahi dituzun funtzioak esportatu
export RK4, f_all!, f_all_rel!, f_all_rel_J!, acceleration

end # module
