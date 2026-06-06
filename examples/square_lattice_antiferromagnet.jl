using SpinWave

include(joinpath(@__DIR__, "plotting.jl"))

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addsite!(model, :B, [0.5, 0.5, 0]; spin=1, moment=[0, 0, -1])
addmatrix!(model, :J, heisenberg(1.0))

addbond!(model, :J, :A, :B, [0, 0, 0])
addbond!(model, :J, :A, :B, [-1, 0, 0])
addbond!(model, :J, :A, :B, [0, -1, 0])
addbond!(model, :J, :A, :B, [-1, -1, 0])

# The path follows near-Γ -> X -> M -> near-Γ.
path = qpath(
    [[0.05, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0.05, 0, 0]];
    points=[41, 41, 41],
    labels=["near Γ", "X", "M", "near Γ"],
)
spec = spinwave(model, path)
grid = broaden(spec, range(0, 4.4; length=120); eta=0.16)

println("ticks:")
show(stdout, "text/plain", path.ticks)
println("\n\nenergies at path ticks:")
show(stdout, "text/plain", round.(spec.energies[:, path.ticks]; digits=4))
println()

if should_write_example_plots()
    out = example_plot_dir()
    write_dispersion_svg(
        joinpath(out, "square-lattice-antiferromagnet-dispersion.svg"),
        path,
        spec.energies;
        title="Square-lattice antiferromagnet dispersion",
        ylabel="energy",
    )
    write_heatmap_svg(
        joinpath(out, "square-lattice-antiferromagnet-intensity.svg"),
        path,
        grid;
        title="Square-lattice antiferromagnet broadened intensity",
        ylabel="energy",
    )
end

nothing
