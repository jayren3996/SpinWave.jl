using SpinWave

include(joinpath(@__DIR__, "plotting.jl"))

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addsite!(model, :B, [0.5, 0, 0]; spin=1, moment=[0, 0, -1])
addmatrix!(model, :J, heisenberg(1.0))
addbond!(model, :J, :A, :B, [0, 0, 0])
addbond!(model, :J, :B, :A, [1, 0, 0])

# Start slightly away from the exact Goldstone point for stable mode sorting.
path = qpath([[0.02, 0, 0], [0.5, 0, 0]]; points=81, labels=["near Γ", "X"])
spec = spinwave(model, path)

samples = [1, 21, 41, 61, 81]
println("sampled q-points (rlu):")
show(stdout, "text/plain", path.q_rlu[:, samples])
println("\n\nsampled mode energies:")
show(stdout, "text/plain", round.(spec.energies[:, samples]; digits=4))
println()

if should_write_example_plots()
    write_dispersion_svg(
        joinpath(example_plot_dir(), "antiferromagnetic-chain-dispersion.svg"),
        path,
        spec.energies;
        title="Antiferromagnetic chain dispersion",
        ylabel="energy",
    )
end

nothing
