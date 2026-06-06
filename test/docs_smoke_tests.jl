@testset "docs scaffold" begin
    root = dirname(@__DIR__)
    @test isfile(joinpath(root, "docs", "Project.toml"))
    @test isfile(joinpath(root, "docs", "make.jl"))
    @test isfile(joinpath(root, "docs", "src", "index.md"))
end
