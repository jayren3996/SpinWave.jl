@testset "diagonalization" begin
    H = ComplexF64[2 0; 0 2]
    result = SpinWave.diagonalize_bosonic(H)
    @test result.energies ≈ [2.0]
    @test size(result.modes) == (2, 1)
    @test result.metric_residual <= 1e-8

    Hzero = zeros(ComplexF64, 2, 2)
    @test SpinWave.diagonalize_bosonic(Hzero).energies ≈ [0.0]

    @test_throws ArgumentError SpinWave.diagonalize_bosonic([1.0 1.0; 0.0 1.0])
    @test_throws SpinWave.UnstableHamiltonianError SpinWave.diagonalize_bosonic([1.0 0.0; 0.0 -1.0])

    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    qs = [0.25, 0.5]
    energies = [
        only(SpinWave.diagonalize_bosonic(SpinWave.bosonic_hamiltonian(model, [q, 0, 0])).energies)
        for q in qs
    ]
    expected = @. 2 * (1 - cos(2pi * qs))
    @test energies ≈ expected atol=1e-8

    afm = SpinModel(lattice([1, 1, 1]))
    addsite!(afm, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addsite!(afm, :B, [0.5, 0, 0]; spin=1, moment=[0, 0, -1])
    addmatrix!(afm, :J, heisenberg(1.0))
    addbond!(afm, :J, :A, :B, [0, 0, 0])
    addbond!(afm, :J, :B, :A, [1, 0, 0])

    afm_result = SpinWave.diagonalize_bosonic(SpinWave.bosonic_hamiltonian(afm, [0.1, 0, 0]))
    @test afm_result.energies ≈ fill(2sin(0.1pi), 2) atol=1e-8
    @test afm_result.metric_residual <= 1e-8
end
