@testset "spinwave spectrum" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[3])
    spec = spinwave(model, path)
    @test size(spec.energies) == (1, 3)
    @test spec.energies[:, 1] ≈ [2.0]
    @test spec.energies[:, end] ≈ [4.0]
    @test size(spec.correlations) == (3, 3, 1, 3)
    @test spec.qpoints === path

    weights = intensity(spec)
    @test size(weights) == (1, 3)
    @test all(>=(0), weights)

    grid = broaden(spec, range(0, 4; length=25); eta=0.1)
    @test size(grid.intensity) == (25, 3)
    @test grid.omegas[1] == 0
    @test grid.qpoints === path

    @test_throws DomainError broaden(spec, range(0, 4; length=25); eta=0)
end
