# SpinWave.jl Rebuild Design

## Context

SpinWave.jl is an old, compact Julia package for linear spin-wave calculations.
The current API exposes low-level arrays directly:

```julia
spinsystem(theta, phi, ind, vec, mat, t)
spinwave(system, kpoints, omegas)
```

That interface is hard to validate and does not map cleanly to how users think
about magnetic models. The current implementation also has correctness defects:

- the Hamiltonian onsite term is allocated with the number of bonds instead of
  the number of magnetic sites and is updated before initialization,
- the solver silently regularizes failed Cholesky factorizations by adding a
  fixed diagonal shift,
- site counts can disagree between spin directions and sublattice positions,
- invalid dimensions, bond indices, interaction matrices, k-point paths, and
  damping values are not rejected early,
- the returned correlation array has no stable documented orientation,
- `Pkg.test()` does not run because there is no `test/runtests.jl`.

The closest mature MATLAB API reference uses this workflow:

1. define a lattice,
2. add magnetic atoms,
3. generate or add couplings,
4. define reusable interaction matrices,
5. define a magnetic structure,
6. compute a rich spin-wave spectrum.

The rebuild should borrow that staged model while using Julia types, explicit
validation, and stable result objects rather than MATLAB-style mutable option
bags.

## Goals

- Replace the low-level public API with a logical model-building interface.
- Keep the first rebuilt solver focused on commensurate magnetic structures and
  bilinear spin interactions.
- Make invalid models fail with clear `ArgumentError` or `DomainError`
  exceptions before numerical work begins.
- Make bosonic dynamical-matrix diagonalization explicit and diagnostic: no
  hidden mutation or blind diagonal shifts.
- Return a named `SpinWaveSpectrum` object with stable dimensions and metadata.
- Add deterministic tests that verify constructors, validation, k-paths,
  Hamiltonian assembly, diagonalization failure modes, and at least one analytic
  model.
- Add a Documenter.jl site with a manual, examples, and API reference, following
  the style of the user's newer Julia packages.
- Add package hygiene: compat bounds, `Test` target, `.gitignore`, CI, docs
  workflow, README rewrite, and examples.

## Non-Goals

- Full feature parity with the reference MATLAB package.
- Crystallographic space-group generation.
- Incommensurate or multi-k magnetic structures.
- Automatic magnetic-structure optimization.
- Neutron cross-section, form-factor, twinning, instrumental convolution, or
  plotting APIs.
- Preserving old `spinsystem`, `kpoints`, or tuple-returning `spinwave`
  compatibility.

Those can be added later on top of a validated core.

## Recommended Approach

Use a compact but typed Julia model:

```julia
model = SpinModel(lattice([1, 1, 1], [90, 90, 90]))
addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
addmatrix!(model, :J1, heisenberg(1.0))
addbond!(model, :J1, :A, :A, [1, 0, 0])

path = qpath([[0, 0, 0], [1, 0, 0], [1, 1, 0]]; points=[100, 100])
spectrum = spinwave(model, path)
```

This is closer to the staged reference workflow than the existing array API, but
keeps the Julia surface straightforward:

- use typed structs instead of stringly typed modes,
- use mutating builder methods for model construction,
- validate before allocation-heavy numerical work,
- return data objects instead of tuples with ambiguous orientation.

## Public Types

### `Lattice`

Stores real-space lattice vectors and reciprocal lattice vectors.

Required constructors:

- `lattice(a::Real, b::Real, c::Real; angles=(90, 90, 90), degrees=true)`
- `lattice(lengths::AbstractVector; angles=(90, 90, 90), degrees=true)`
- `lattice(vectors::AbstractMatrix)`

Conventions:

- Real-space vectors are stored as columns in a `3 x 3` matrix.
- Reciprocal vectors include the `2pi` factor.
- Fractional coordinates are interpreted in lattice units.
- Reciprocal-lattice-unit q-points are converted with the reciprocal basis.

Validation:

- all lengths must be positive,
- all angles must define a non-singular cell,
- vector matrices must be `3 x 3`, finite, and non-singular.

### `SpinSite`

Stores one magnetic basis site:

- `label::Symbol`,
- `position::SVector{3,Float64}` or equivalent fixed-size vector,
- `spin::Float64`,
- `moment::SVector{3,Float64}`.

Validation:

- labels must be unique in a model,
- spin length must be positive,
- moment must be finite and nonzero.

The moment is normalized for local-frame construction. The spin length scales the
quadratic Hamiltonian and correlations.

### `ExchangeMatrix`

Stores a reusable `3 x 3` interaction tensor:

- `heisenberg(J)` creates `J * I`,
- `exchange_matrix(M)` stores a full matrix,
- `dm(D)` creates the antisymmetric Dzyaloshinskii-Moriya tensor.

Validation:

- matrices must be finite,
- non-Hermitian or nonsymmetric real matrices are accepted only through explicit
  constructors such as `dm`, so users cannot accidentally pass a typo as a
  general exchange matrix.

### `Bond`

Stores one directed bond:

- matrix label,
- source site label or index,
- target site label or index,
- integer cell offset,
- optional scale factor.

The Hamiltonian assembly adds the Hermitian conjugate partner automatically, so
users add each physical bond once.

Validation:

- source and target sites must exist,
- offset must be a length-3 integer vector,
- matrix label must exist,
- scale must be finite.

### `SpinModel`

Owns the lattice, sites, exchange matrices, and bonds.

The model is mutable during construction:

```julia
model = SpinModel(lat)
addsite!(model, ...)
addmatrix!(model, ...)
addbond!(model, ...)
```

Before solving, `validate(model)` checks cross-object invariants and returns the
model or throws a clear exception. `spinwave` always validates unless
`check=false` is passed.

### `QPath`

Stores reciprocal-space scan points and labels:

- `qpath(vertices; points)` accepts vertices in reciprocal lattice units and
  expands linearly,
- each segment includes its starting point and excludes duplicated interior
  vertices; the final endpoint is included exactly once,
- stored fields include `q_rlu`, `q_cartesian`, `ticks`, and `labels`.

Validation:

- path vertices must be finite `3`-vectors,
- segment point counts must be integers `>= 2`,
- the number of segment counts must equal number of vertices minus one.

### `SpinWaveSpectrum`

Result object returned by `spinwave`:

- `energies::Matrix{Float64}` with shape `(nmode, nq)`,
- `correlations::Array{ComplexF64,4}` with shape `(3, 3, nmode, nq)`,
- `qpoints::QPath`,
- `status::Vector{Symbol}` or diagnostics per q-point,
- optional `hamiltonians` and `transforms` when requested.

The object may also expose convenience methods:

- `intensity(spec; components=:trace)` returns `(nmode, nq)` mode weights,
- `broaden(spec, omegas; eta, components=:trace)` returns an energy grid with
  shape `(length(omegas), nq)`.

Post-processing stays separate from raw LSWT output.

## Solver Design

### Local Frames

For every spin direction, construct an orthonormal local frame:

- local z-axis points along the ordered moment,
- local x/y axes span the transverse plane,
- the complex transverse vector defines the Holstein-Primakoff creation and
  annihilation basis.

Frame construction must avoid singular behavior for moments close to global
axes by choosing a stable reference vector.

Tests should verify:

- frame vectors are orthonormal,
- `z` matches the normalized moment,
- right-handed orientation is consistent.

### Hamiltonian Assembly

For every bond `(i, j, R, J)`:

- rotate `J` into local spin frames,
- include spin-length factors,
- add particle-particle, particle-hole, hole-particle, and hole-hole blocks,
- apply phase `exp(-2pi * im * dot(q_rlu, r_j + R - r_i))`,
- add the conjugate counterpart once.

The assembled bosonic Hamiltonian must be Hermitian to numerical tolerance. If
not, throw an error with the q-point index and Hermiticity residual.

### Diagonalization

Use dense bosonic dynamical-matrix diagonalization for positive-semidefinite
quadratic forms:

1. verify `H` is Hermitian,
2. verify the Hermitian form is positive semidefinite within tolerance,
3. diagonalize `ηH`,
4. reject complex mode eigenvalues above tolerance,
5. sort the positive-norm branch by nonnegative energy,
6. metric-normalize the selected modes,
7. verify the selected modes satisfy metric normalization.

No automatic diagonal shift is applied. A later release can add a stricter
Cholesky/Colpa path or explicit `regularization` keyword, but the default must
not hide instability.

### Correlations

The rebuilt package should compute component-resolved mode correlations
`S^{ab}(q, omega)` in a stable shape:

```julia
correlations[a, b, mode, iq]
```

The initial implementation should prioritize correctness and diagnostics over
speed. It may be dense and allocation-heavy if the API and tests are clear.

## Tests

Create a standard `test/runtests.jl` and split tests into small files:

- `lattice_tests.jl`,
- `model_tests.jl`,
- `qpath_tests.jl`,
- `local_frame_tests.jl`,
- `hamiltonian_tests.jl`,
- `diagonalization_tests.jl`,
- `spinwave_tests.jl`,
- `docs_smoke_tests.jl`.

Required cases:

- invalid lattice lengths, angles, singular vector matrices,
- duplicate site labels, invalid moments, invalid spin lengths,
- unknown bond matrix/site labels, bad offsets,
- scalar Heisenberg and full-matrix exchange constructors,
- q-path endpoint de-duplication and tick positions,
- local frame orthonormality,
- Hermiticity of a simple ferromagnet Hamiltonian,
- rejection of non-Hermitian and indefinite test matrices,
- analytic ferromagnetic chain dispersion:
  `omega(q) = 2 J S (1 - cos(q))` for a one-site nearest-neighbor chain,
- spectrum object dimensions,
- energy-grid broadening dimensions.

Tests should not depend on plotting packages.

## Documentation

Add a Documenter.jl site:

- `docs/Project.toml`,
- `docs/make.jl`,
- `docs/src/index.md`,
- `docs/src/getting-started.md`,
- `docs/src/manual/conventions.md`,
- `docs/src/manual/models.md`,
- `docs/src/manual/spectra.md`,
- `docs/src/examples/ferromagnetic-chain.md`,
- `docs/src/examples/square-lattice-antiferromagnet.md`,
- `docs/src/reference/api.md`.

The docs should explain:

- lattice coordinate conventions,
- reciprocal lattice units vs Cartesian q-vectors,
- how to build sites, matrices, and bonds,
- current solver scope and limitations,
- shape conventions for energies, correlations, and broadened intensity,
- numerical assumptions behind bosonic dynamical-matrix diagonalization.

Examples must execute in the docs build when practical and avoid `PyPlot`.

## Repository Hygiene

Update or create:

- `Project.toml` with `[compat]`, `[extras]`, and `[targets]`,
- `.gitignore` for Julia build products, manifests, `.DS_Store`, docs build,
  coverage, and worktrees,
- `README.md` with installation, scope, quick start, docs links, and repository
  map,
- `.github/workflows/CI.yml`,
- `.github/workflows/documentation.yml`,
- `examples/` for runnable scripts that are not tests.

The current `Manifest.toml` and `.DS_Store` files should be removed from git
tracking as cleanup in the rebuild branch.

## Implementation Strategy

Use test-driven development:

1. Add package metadata and a failing test scaffold.
2. Implement lattice and q-path types.
3. Implement model-building types and validation.
4. Implement local frames and interaction constructors.
5. Implement Hamiltonian assembly for commensurate bilinear models.
6. Implement bosonic dynamical-matrix diagonalization with explicit failure
   modes.
7. Implement `spinwave`, `SpinWaveSpectrum`, intensity, and broadening.
8. Add docs and README.
9. Add CI/docs workflows.
10. Run package tests and docs build.

The implementation should prefer clarity and correctness over compatibility with
the old API.

## Success Criteria

- `Pkg.test()` runs and passes locally.
- `julia --project=docs docs/make.jl` builds without Documenter errors.
- The quick-start example in the README is executable.
- The analytic ferromagnetic chain test passes within a documented tolerance.
- Invalid inputs produce clear exceptions, not undefined memory reads or
  unrelated linear algebra stack traces.
- The Documenter site exposes manual pages and API docstrings.
