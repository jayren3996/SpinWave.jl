using SpinWave
using LinearAlgebra
using Test

@testset "SpinWave.jl" begin
    include("lattice_tests.jl")
    include("qpath_tests.jl")
end
