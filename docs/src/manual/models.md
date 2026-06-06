# Models

SpinWave.jl models are built with explicit, validated pieces.

## Sites

[`addsite!`](@ref) adds a magnetic basis site:

```@example models
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
length(model.sites)
```

The ordered moment defines the local spin frame used in the Holstein-Primakoff
expansion. The vector is normalized internally for frame construction; the
`spin` keyword carries the spin length.

## Exchange Matrices

Reusable matrices are registered by label:

```@example models
addmatrix!(model, :J, heisenberg(-1.0))
addmatrix!(model, :D, dm([0, 0, 0.1]))
sort(collect(keys(model.matrices)))
```

Use [`heisenberg`](@ref) for isotropic exchange, [`dm`](@ref) for an
antisymmetric Dzyaloshinskii-Moriya tensor, and [`exchange_matrix`](@ref) for a
full `3 x 3` matrix.

## Bonds

[`addbond!`](@ref) adds one directed physical bond:

```@example models
addbond!(model, :J, :A, :A, [1, 0, 0])
validate(model)
```

The Hermitian counterpart is added during Hamiltonian assembly, so do not add
both `[1, 0, 0]` and `[-1, 0, 0]` for the same physical nearest-neighbor shell
unless that is intentionally a different interaction.
