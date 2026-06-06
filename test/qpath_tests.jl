@testset "qpath" begin
    path = qpath([[0, 0, 0], [1, 0, 0], [1, 1, 0]]; points=[3, 4])
    @test size(path.q_rlu) == (3, 6)
    @test path.q_rlu[:, 1] ≈ [0, 0, 0]
    @test path.q_rlu[:, 3] ≈ [1, 0, 0]
    @test path.q_rlu[:, end] ≈ [1, 1, 0]
    @test path.q_cartesian === nothing
    @test path.ticks == [1, 3, 6]
    @test path.labels == ["", "", ""]

    labelled = qpath([[0, 0, 0], [0.5, 0, 0]]; points=[2], labels=["Γ", "X"])
    @test labelled.labels == ["Γ", "X"]

    cartesian = qpath([[0, 0, 0], [0.5, 0, 0]], lattice([1, 1, 1]); points=[2])
    @test cartesian.q_cartesian[:, 1] ≈ [0, 0, 0]
    @test cartesian.q_cartesian[:, end] ≈ [pi, 0, 0]

    @test_throws ArgumentError qpath([[0, 0, 0], [1, 0, 0]]; points=[1])
    @test_throws ArgumentError qpath([[0, 0, 0], [1, 0, 0]]; points=[2, 2])
    @test_throws ArgumentError qpath([[0, 0], [1, 0]]; points=[2])
end
