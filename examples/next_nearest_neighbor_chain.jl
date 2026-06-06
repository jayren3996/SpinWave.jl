using SpinWave

include(joinpath(@__DIR__, "plotting.jl"))

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J1, heisenberg(-1.0))
addmatrix!(model, :J2, heisenberg(-0.25))
addbond!(model, :J1, :A, :A, [1, 0, 0])
addbond!(model, :J2, :A, :A, [2, 0, 0])

path = qpath([[0, 0, 0], [0.5, 0, 0]]; points=101, labels=["Γ", "X"])
spec = spinwave(model, path)

samples = [1, 26, 51, 76, 101]
println("sampled mode energies:")
show(stdout, "text/plain", round.(spec.energies[:, samples]; digits=4))
println()

if should_write_example_plots()
    write_dispersion_svg(
        joinpath(example_plot_dir(), "next-nearest-neighbor-chain-dispersion.svg"),
        path,
        spec.energies;
        title="Ferromagnetic chain with first and second neighbors",
        ylabel="energy",
    )
end

nothing
