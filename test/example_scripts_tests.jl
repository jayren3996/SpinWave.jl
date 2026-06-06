@testset "example scripts" begin
    root = dirname(@__DIR__)
    expected_assets = Dict(
        "ferromagnetic_chain.jl" => [
            "ferromagnetic-chain-dispersion.svg",
            "ferromagnetic-chain-intensity.svg",
        ],
        "antiferromagnetic_chain.jl" => [
            "antiferromagnetic-chain-dispersion.svg",
        ],
        "square_lattice_antiferromagnet.jl" => [
            "square-lattice-antiferromagnet-dispersion.svg",
            "square-lattice-antiferromagnet-intensity.svg",
        ],
        "anisotropic_ferromagnetic_chain.jl" => [
            "anisotropic-ferromagnetic-chain-dispersion.svg",
            "anisotropic-ferromagnetic-chain-components.svg",
        ],
        "next_nearest_neighbor_chain.jl" => [
            "next-nearest-neighbor-chain-dispersion.svg",
        ],
    )

    for (script, assets) in sort(collect(expected_assets))
        @testset "$script" begin
            path = joinpath(root, "examples", script)
            @test isfile(path)
            mktempdir() do plot_dir
                withenv(
                    "SPINWAVE_WRITE_EXAMPLE_PLOTS" => "1",
                    "SPINWAVE_EXAMPLE_PLOT_DIR" => plot_dir,
                ) do
                    @test include(path) === nothing
                end
                for asset in assets
                    asset_path = joinpath(plot_dir, asset)
                    @test isfile(asset_path)
                    @test filesize(asset_path) > 500
                    @test occursin("<svg", read(asset_path, String))
                end
            end
        end
    end
end
