# Ferromagnetic Chain

This example computes the one-site nearest-neighbor ferromagnetic chain.

```@example chain
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0, 0, 0], [0.5, 0, 0]]; points=[6], labels=["Γ", "X"])
spec = spinwave(model, path)

round.(spec.energies; digits=4)
```

The zero at `Γ` is the ferromagnetic Goldstone mode. A lightweight energy grid
can be produced without any plotting dependency:

```@example chain
grid = broaden(spec, range(0, 4; length=21); eta=0.15)
maximum(grid.intensity)
```
