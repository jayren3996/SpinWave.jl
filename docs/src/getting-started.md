# Getting Started

This page builds a one-site ferromagnetic chain. The Hamiltonian convention is
that a negative Heisenberg coupling is ferromagnetic:

```@example getting-started
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

validate(model)
```

The bond above represents the interaction from site `:A` in one unit cell to
site `:A` in the neighboring cell displaced by `[1, 0, 0]`. The Hermitian
counterpart is added during Hamiltonian assembly.

Build a reciprocal-space path and solve:

```@example getting-started
path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[4])
spec = spinwave(model, path)

spec.energies
```

For this model the analytic dispersion is
``\omega(q) = 2JS(1-\cos(2\pi q))`` with ``J=1`` after taking the absolute
ferromagnetic scale. The values at ``q=1/4`` and ``q=1/2`` are therefore `2`
and `4`.

To build an energy grid:

```@example getting-started
grid = broaden(spec, range(0, 4; length=9); eta=0.1)
size(grid.intensity)
```
