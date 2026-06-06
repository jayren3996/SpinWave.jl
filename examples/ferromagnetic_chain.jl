using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0, 0, 0], [0.5, 0, 0]]; points=[6], labels=["Γ", "X"])
spec = spinwave(model, path)

println("q-points (rlu):")
show(stdout, "text/plain", path.q_rlu)
println("\n\nmode energies:")
show(stdout, "text/plain", spec.energies)
println()
