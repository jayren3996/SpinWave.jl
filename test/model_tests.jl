@testset "model builders" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J1, heisenberg(-1.0))
    addbond!(model, :J1, :A, :A, [1, 0, 0])

    @test length(model.sites) == 1
    @test model.sites[1].label == :A
    @test haskey(model.matrices, :J1)
    @test model.matrices[:J1].matrix ≈ -I(3)
    @test length(model.bonds) == 1
    @test validate(model) === model

    dm_matrix = dm([1, 2, 3]).matrix
    @test dm_matrix ≈ [0 -3 2; 3 0 -1; -2 1 0]

    @test_throws ArgumentError addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    @test_throws DomainError addsite!(SpinModel(lattice([1, 1, 1])), :B, [0, 0, 0]; spin=0, moment=[0, 0, 1])
    @test_throws DomainError addsite!(SpinModel(lattice([1, 1, 1])), :B, [0, 0, 0]; spin=1, moment=[0, 0, 0])
    @test_throws ArgumentError addbond!(model, :missing, :A, :A, [1, 0, 0])
    @test_throws ArgumentError addbond!(model, :J1, :missing, :A, [1, 0, 0])
    @test_throws ArgumentError exchange_matrix([1 2; 3 4])
    @test_throws ArgumentError exchange_matrix([1 2 0; 0 1 0; 0 0 1])
    @test_throws ArgumentError dm([1, 2])

    malformed = SpinModel(lattice([1, 1, 1]))
    addsite!(malformed, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(malformed, :J, heisenberg(-1))
    push!(malformed.bonds, Bond(:J, 99, 1, [0, 0, 0], 1.0))
    @test_throws ArgumentError validate(malformed)
end
