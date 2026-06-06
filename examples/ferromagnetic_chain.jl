using SpinWave

include(joinpath(@__DIR__, "plotting.jl"))

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0, 0, 0], [0.5, 0, 0]]; points=101, labels=["Γ", "X"])
spec = spinwave(model, path)
grid = broaden(spec, range(0, 4.4; length=120); eta=0.12)

samples = [1, 26, 51, 76, 101]
println("sampled q-points (rlu):")
show(stdout, "text/plain", path.q_rlu[:, samples])
println("\n\nsampled mode energies:")
show(stdout, "text/plain", round.(spec.energies[:, samples]; digits=4))
println()

if should_write_example_plots()
    out = example_plot_dir()
    write_dispersion_svg(
        joinpath(out, "ferromagnetic-chain-dispersion.svg"),
        path,
        spec.energies;
        title="Ferromagnetic chain dispersion",
        ylabel="energy",
    )
    write_heatmap_svg(
        joinpath(out, "ferromagnetic-chain-intensity.svg"),
        path,
        grid;
        title="Ferromagnetic chain broadened intensity",
        ylabel="energy",
    )
end

nothing
