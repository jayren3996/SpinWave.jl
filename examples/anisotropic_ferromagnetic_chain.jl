using SpinWave

include(joinpath(@__DIR__, "plotting.jl"))

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(
    model,
    :Jdiag,
    exchange_matrix([
        -3.0 0.0 0.0
        0.0 -4.0 0.0
        0.0 0.0 -5.0
    ]),
)
addbond!(model, :Jdiag, :A, :A, [1, 0, 0])

path = qpath([[0, 0, 0], [0.5, 0, 0]]; points=81, labels=["Γ", "X"])
spec = spinwave(model, path)
sx = real.(spec.correlations[1, 1, 1, :])
sy = real.(spec.correlations[2, 2, 1, :])

samples = [1, 21, 41, 61, 81]
println("mode energies:")
show(stdout, "text/plain", round.(spec.energies[:, samples]; digits=4))
println("\n\nSxx weight:")
show(stdout, "text/plain", round.(sx[samples]; digits=4))
println()

if should_write_example_plots()
    out = example_plot_dir()
    write_dispersion_svg(
        joinpath(out, "anisotropic-ferromagnetic-chain-dispersion.svg"),
        path,
        spec.energies;
        title="Anisotropic ferromagnetic chain dispersion",
        ylabel="energy",
    )
    write_series_svg(
        joinpath(out, "anisotropic-ferromagnetic-chain-components.svg"),
        collect(1:length(path)),
        [sx, sy];
        labels=["Sxx", "Syy"],
        title="Anisotropic chain component weights",
        ylabel="correlation weight",
        xticks=path.ticks,
        xticklabels=path.labels,
        ymin=0.0,
    )
end

nothing
