using Test
using GravitationSimulation
using StaticArrays
using LinearAlgebra

@testset "GravitationSimulation.jl" begin
    include("test_equations.jl")
    include("test_integrator.jl")
end
