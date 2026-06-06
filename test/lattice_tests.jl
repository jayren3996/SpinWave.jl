@testset "lattice" begin
    lat = lattice([1, 2, 3]; angles=(90, 90, 90))
    @test size(lat.vectors) == (3, 3)
    @test size(lat.reciprocal) == (3, 3)
    @test lat.reciprocal[:, 1] ≈ [2pi, 0, 0]
    @test lat.reciprocal[:, 2] ≈ [0, pi, 0]
    @test lat.reciprocal[:, 3] ≈ [0, 0, 2pi / 3]

    cubic = lattice(2, 2, 2)
    @test cubic.vectors ≈ 2I(3)
    @test reciprocal_cartesian(cubic, [0.5, 0, 0]) ≈ [pi / 2, 0, 0]

    @test_throws DomainError lattice([0, 1, 1])
    @test_throws ArgumentError lattice([1, 1])
    @test_throws ArgumentError lattice([1 0; 0 1])
end
