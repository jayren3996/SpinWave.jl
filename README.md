<div align="center">

# SpinWave.jl

Typed, validated linear spin-wave calculations for Julia.

[![CI](https://github.com/jayren3996/SpinWave.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jayren3996/SpinWave.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/jayren3996/SpinWave.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/jayren3996/SpinWave.jl/actions/workflows/documentation.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://jayren3996.github.io/SpinWave.jl/)

**[Read the documentation](https://jayren3996.github.io/SpinWave.jl/)**

</div>

`SpinWave.jl` is a compact Julia package for linear spin-wave calculations on
explicit bilinear spin models. It now uses a typed model-building workflow:
define a lattice, add magnetic sites, register exchange matrices, add bonds,
choose a reciprocal-space path, and compute a named spectrum object.

The rebuilt core is intentionally focused. It supports dense commensurate
models with explicit bonds and validates inputs before solving. Space-group bond
generation, incommensurate structures, form factors, twinning, instrumental
resolution, and plotting are left for later layers.

## Installation

```julia
pkg> add https://github.com/jayren3996/SpinWave.jl
```

## Quick Start

```julia
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[5], labels=["q=1/4", "q=1/2"])
spec = spinwave(model, path)

spec.energies
```

`spec.energies` has shape `(mode, q)`. Component-resolved correlations are in
`spec.correlations[a, b, mode, q]`.

Create a broadened energy grid without plotting dependencies:

```julia
grid = broaden(spec, range(0, 4; length=100); eta=0.1)
size(grid.intensity) # (omega, q)
```

## Repository Layout

| Path | Purpose |
| --- | --- |
| `src/` | Core model, Hamiltonian, diagonalization, and spectrum code |
| `test/` | Deterministic unit and smoke tests |
| `docs/` | Documenter.jl manual and API reference |
| `examples/` | Runnable scripts that are not tests |

## References

The public workflow is inspired by the staged model construction used by
[SpinW](https://spinw.org/spinwdoc/), but the implementation is Julian: typed
objects, explicit validation, and stable result shapes.
