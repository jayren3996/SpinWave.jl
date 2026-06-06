@testset "hamiltonian assembly" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    H0 = SpinWave.bosonic_hamiltonian(model, [0, 0, 0])
    @test size(H0) == (2, 2)
    @test ishermitian(H0)
    @test H0 ≈ zeros(ComplexF64, 2, 2) atol=1e-12

    Hpi = SpinWave.bosonic_hamiltonian(model, [0.5, 0, 0])
    @test ishermitian(Hpi)
    @test Hpi ≈ [4 0; 0 4]

    @test_throws ArgumentError SpinWave.bosonic_hamiltonian(model, [0, 0])
end
