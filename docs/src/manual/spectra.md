# Spectra

[`spinwave`](@ref) validates a [`SpinModel`](@ref), assembles the dense bosonic
Hamiltonian at each q-point, diagonalizes it, and returns a
[`SpinWaveSpectrum`](@ref).

```@example spectra
using SpinWave

model = SpinModel(lattice([1, 1, 1]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J, heisenberg(-1.0))
addbond!(model, :J, :A, :A, [1, 0, 0])

path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[3])
spec = spinwave(model, path)

size(spec.energies), size(spec.correlations)
```

[`intensity`](@ref) currently supports `components=:trace`, summing the diagonal
spin-correlation components for each mode:

```@example spectra
intensity(spec)
```

[`broaden`](@ref) applies Lorentzian broadening:

```@example spectra
grid = broaden(spec, range(0, 4; length=11); eta=0.1)
size(grid.intensity)
```

## Numerical Assumptions

The current dense solver rejects non-Hermitian Hamiltonians and Hermitian
Hamiltonians with clearly negative eigenvalues. Semidefinite Hamiltonians are
allowed when a positive-norm zero-mode branch can be identified numerically. It
diagonalizes the bosonic dynamical matrix, selects the nonnegative
positive-norm branch, and metric-normalizes those modes. Complex
dynamical-matrix eigenvalues above tolerance are treated as unstable or
inconsistent input.
