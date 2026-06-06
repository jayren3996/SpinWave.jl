# SpinWave.jl

`SpinWave.jl` builds and solves compact linear spin-wave models in Julia. It is
organized around an explicit staged workflow: define a lattice, add magnetic
sites, register exchange matrices, add bonds, choose a reciprocal-space path,
and compute the spectrum.

The current rebuilt core focuses on commensurate bilinear models. It validates
model dimensions and numerical assumptions before solving, and returns a named
[`SpinWaveSpectrum`](@ref) object instead of an ambiguous tuple.

## A Minimal Run

```@example home
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[5], labels=["q=1/4", "q=1/2"])
spec = spinwave(model, path)

round.(spec.energies; digits=3)
```

## Where To Start

- [Getting Started](getting-started.md) walks through the ferromagnetic chain.
- [Conventions](manual/conventions.md) defines lattice, q-point, and array
  shapes.
- [Models](manual/models.md) explains sites, exchange matrices, and bonds.
- [Spectra](manual/spectra.md) covers `SpinWaveSpectrum`, `intensity`, and
  `broaden`.
- [Examples](examples/ferromagnetic-chain.md) show runnable chain and lattice
  models, including rendered SVG spectra generated from the example scripts.
- [API Reference](reference/api.md) lists public docstrings.

## Current Scope

SpinWave.jl currently supports dense, commensurate, bilinear spin-wave models
with explicitly supplied bonds. It does not yet implement space-group bond
generation, incommensurate magnetic structures, form factors, twinning,
instrumental convolution, or plotting.
