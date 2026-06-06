# SpinWave.jl Rebuild Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild SpinWave.jl around validated model objects, a stable spin-wave spectrum API, real tests, and a Documenter.jl manual.

**Architecture:** Replace the old low-level array API with typed model-building layers: lattice and q-path primitives, mutable model builders, local-frame and Hamiltonian assembly internals, explicit bosonic diagonalization, and a named `SpinWaveSpectrum` result. Documentation and examples sit on top of the same public API used by tests.

**Tech Stack:** Julia 1.10+, LinearAlgebra stdlib, Test stdlib, Documenter.jl v1, GitHub Actions.

---

## File Structure

- Modify: `Project.toml` for compat/test target.
- Create: `.gitignore` for Julia/docs/worktree artifacts.
- Replace: `src/SpinWave.jl` as the module entrypoint and export list.
- Create: `src/Lattices.jl` for `Lattice`, `lattice`, and coordinate conversion.
- Create: `src/QPaths.jl` for `QPath` and `qpath`.
- Create: `src/Models.jl` for `SpinSite`, `ExchangeMatrix`, `Bond`, `SpinModel`, and builder methods.
- Create: `src/Frames.jl` for local spin frames.
- Create: `src/Hamiltonians.jl` for bosonic Hamiltonian assembly.
- Replace: `src/Diagonalization.jl` with validated Colpa/dynamic diagonalization helpers.
- Create: `src/Spectra.jl` for `SpinWaveSpectrum`, `spinwave`, `intensity`, and `broaden`.
- Stop including: `src/HBlock.jl`, `src/WBlock.jl`, and `src/Solver.jl`; leave them deleted or unreferenced.
- Create: `test/runtests.jl` plus focused test files.
- Move: plotting/demo scripts to `examples/` if retained.
- Create: `docs/Project.toml`, `docs/make.jl`, and `docs/src/**`.
- Modify: `README.md`.
- Create: `.github/workflows/CI.yml` and `.github/workflows/documentation.yml`.

## Task 1: Package Hygiene And Test Scaffold

**Files:**
- Modify: `Project.toml`
- Create: `.gitignore`
- Create: `test/runtests.jl`
- Create: `test/lattice_tests.jl`
- Create: `test/qpath_tests.jl`

- [ ] **Step 1.1: Write the failing test scaffold**

Create `test/runtests.jl`:

```julia
using SpinWave
using Test

@testset "SpinWave.jl" begin
    include("lattice_tests.jl")
    include("qpath_tests.jl")
end
```

Create `test/lattice_tests.jl`:

```julia
@testset "lattice" begin
    lat = lattice([1, 2, 3]; angles=(90, 90, 90))
    @test size(lat.vectors) == (3, 3)
    @test size(lat.reciprocal) == (3, 3)
    @test lat.reciprocal[:, 1] ≈ [2pi, 0, 0]
    @test_throws DomainError lattice([0, 1, 1])
end
```

Create `test/qpath_tests.jl`:

```julia
@testset "qpath" begin
    path = qpath([[0, 0, 0], [1, 0, 0], [1, 1, 0]]; points=[3, 4])
    @test size(path.q_rlu) == (3, 6)
    @test path.q_rlu[:, 1] ≈ [0, 0, 0]
    @test path.q_rlu[:, 3] ≈ [1, 0, 0]
    @test path.q_rlu[:, end] ≈ [1, 1, 0]
    @test path.ticks == [1, 3, 6]
    @test_throws ArgumentError qpath([[0, 0, 0], [1, 0, 0]]; points=[1])
end
```

- [ ] **Step 1.2: Run tests to verify they fail**

Run:

```bash
/Users/ren/.local/bin/julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because `lattice` and `qpath` are not defined.

- [ ] **Step 1.3: Add package metadata and ignores**

Update `Project.toml`:

```toml
name = "SpinWave"
uuid = "5d8f5327-b1cb-4fc0-964a-d944ba32196b"
authors = ["JayRen <jieren3996@foxmail.com>"]
version = "0.3.0"

[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
LinearAlgebra = "1"
julia = "1.10"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```

Create `.gitignore`:

```gitignore
.DS_Store
*.cov
*.ji
*.mem
/Manifest.toml
/docs/Manifest.toml
/docs/build/
/.worktrees/
```

- [ ] **Step 1.4: Implement lattice and q-path primitives**

Create `src/Lattices.jl` with `Lattice`, `lattice`, and `reciprocal_cartesian`.
Create `src/QPaths.jl` with `QPath` and `qpath`.
Modify `src/SpinWave.jl` to include and export those names.

- [ ] **Step 1.5: Run tests to verify Task 1 passes**

Run:

```bash
/Users/ren/.local/bin/julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Expected: PASS for lattice and qpath testsets.

## Task 2: Model Builders And Validation

**Files:**
- Create: `src/Models.jl`
- Modify: `src/SpinWave.jl`
- Create: `test/model_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 2.1: Write failing model tests**

Add `include("model_tests.jl")` to `test/runtests.jl`.

Create `test/model_tests.jl`:

```julia
@testset "model builders" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J1, heisenberg(-1.0))
    addbond!(model, :J1, :A, :A, [1, 0, 0])

    @test length(model.sites) == 1
    @test haskey(model.matrices, :J1)
    @test length(model.bonds) == 1
    @test validate(model) === model

    @test_throws ArgumentError addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    @test_throws DomainError addsite!(SpinModel(lattice([1, 1, 1])), :B, [0, 0, 0]; spin=0, moment=[0, 0, 1])
    @test_throws ArgumentError addbond!(model, :missing, :A, :A, [1, 0, 0])
    @test_throws ArgumentError exchange_matrix([1 2; 3 4])
end
```

- [ ] **Step 2.2: Run tests to verify they fail**

Run:

```bash
/Users/ren/.local/bin/julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because model API is not defined.

- [ ] **Step 2.3: Implement model types and builders**

Implement:

```julia
mutable struct SpinModel
    lattice::Lattice
    sites::Vector{SpinSite}
    matrices::Dict{Symbol,ExchangeMatrix}
    bonds::Vector{Bond}
end
```

Add:

- `SpinModel(lat::Lattice)`,
- `addsite!`,
- `heisenberg`,
- `exchange_matrix`,
- `dm`,
- `addmatrix!`,
- `addbond!`,
- `validate(model)`.

- [ ] **Step 2.4: Run tests to verify Task 2 passes**

Run package tests. Expected: PASS.

## Task 3: Local Frames And Hamiltonian Assembly

**Files:**
- Create: `src/Frames.jl`
- Create: `src/Hamiltonians.jl`
- Modify: `src/SpinWave.jl`
- Create: `test/local_frame_tests.jl`
- Create: `test/hamiltonian_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 3.1: Write failing frame and Hamiltonian tests**

Add includes to `test/runtests.jl`.

Create `test/local_frame_tests.jl`:

```julia
@testset "local frames" begin
    frame = SpinWave.local_frame([0, 0, 2])
    @test frame.z ≈ [0, 0, 1]
    @test dot(frame.x, frame.y) ≈ 0 atol=1e-12
    @test dot(frame.x, frame.z) ≈ 0 atol=1e-12
    @test dot(frame.y, frame.z) ≈ 0 atol=1e-12
    @test cross(frame.x, frame.y) ≈ frame.z
    @test_throws DomainError SpinWave.local_frame([0, 0, 0])
end
```

Create `test/hamiltonian_tests.jl`:

```julia
@testset "hamiltonian assembly" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    H = SpinWave.bosonic_hamiltonian(model, [0, 0, 0])
    @test size(H) == (2, 2)
    @test ishermitian(H)

    Hpi = SpinWave.bosonic_hamiltonian(model, [0.5, 0, 0])
    @test ishermitian(Hpi)
end
```

- [ ] **Step 3.2: Run tests to verify they fail**

Expected: FAIL because internal frame and Hamiltonian functions are absent.

- [ ] **Step 3.3: Implement frames and dense Hamiltonian assembly**

Implement:

- `LocalFrame`,
- `local_frame(moment)`,
- `local_basis(frame)`,
- `bosonic_hamiltonian(model, q_rlu)`.

Use one physical directed bond and add its Hermitian counterpart during assembly.
Reject a non-Hermitian final matrix with an `ArgumentError`.

- [ ] **Step 3.4: Run tests to verify Task 3 passes**

Run package tests. Expected: PASS.

## Task 4: Diagonalization And Analytic Energies

**Files:**
- Replace: `src/Diagonalization.jl`
- Modify: `src/SpinWave.jl`
- Create: `test/diagonalization_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 4.1: Write failing diagonalization tests**

Create `test/diagonalization_tests.jl`:

```julia
@testset "diagonalization" begin
    H = [2.0 + 0im 0; 0 2.0 + 0im]
    result = SpinWave.diagonalize_bosonic(H)
    @test result.energies ≈ [2.0]

    @test_throws ArgumentError SpinWave.diagonalize_bosonic([1.0 1.0; 0.0 1.0])
    @test_throws SpinWave.UnstableHamiltonianError SpinWave.diagonalize_bosonic([1.0 0.0; 0.0 -1.0])

    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    qs = [0.25, 0.5]
    energies = [only(SpinWave.diagonalize_bosonic(SpinWave.bosonic_hamiltonian(model, [q, 0, 0])).energies) for q in qs]
    expected = @. 2 * (1 - cos(2pi * qs))
    @test energies ≈ expected atol=1e-8
end
```

- [ ] **Step 4.2: Run tests to verify they fail**

Expected: FAIL because diagonalization result/errors are absent or incorrect.

- [ ] **Step 4.3: Implement diagonalization**

Implement:

- `struct UnstableHamiltonianError <: Exception`,
- `struct DiagonalizationResult`,
- `diagonalize_bosonic(H; tol=1e-9)`.

Behavior:

- reject odd-sized matrices,
- reject non-Hermitian matrices,
- reject clearly indefinite Hermitian matrices,
- compute positive magnon energies from the bosonic dynamical matrix,
- reject complex eigenvalues larger than tolerance.

- [ ] **Step 4.4: Run tests to verify Task 4 passes**

Run package tests. Expected: PASS.

## Task 5: Spectrum API And Broadening

**Files:**
- Create: `src/Spectra.jl`
- Modify: `src/SpinWave.jl`
- Create: `test/spinwave_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 5.1: Write failing spectrum tests**

Create `test/spinwave_tests.jl`:

```julia
@testset "spinwave spectrum" begin
    model = SpinModel(lattice([1, 1, 1]))
    addsite!(model, :A, [0, 0, 0]; spin=1, moment=[0, 0, 1])
    addmatrix!(model, :J, heisenberg(-1.0))
    addbond!(model, :J, :A, :A, [1, 0, 0])

    path = qpath([[0.25, 0, 0], [0.5, 0, 0]]; points=[3])
    spec = spinwave(model, path)
    @test size(spec.energies) == (1, 3)
    @test size(spec.correlations) == (3, 3, 1, 3)
    @test spec.qpoints === path

    weights = intensity(spec)
    @test size(weights) == (1, 3)

    grid = broaden(spec, range(0, 4; length=25); eta=0.1)
    @test size(grid.intensity) == (25, 3)
    @test grid.omegas[1] == 0
end
```

- [ ] **Step 5.2: Run tests to verify they fail**

Expected: FAIL because `spinwave`, `SpinWaveSpectrum`, `intensity`, and `broaden` are absent.

- [ ] **Step 5.3: Implement spectrum objects**

Implement:

- `SpinWaveSpectrum`,
- `EnergyGrid`,
- `spinwave(model, path; check=true, correlations=true)`,
- `intensity(spec; components=:trace)`,
- `broaden(spec, omegas; eta, components=:trace)`.

Initial correlations may use a conservative component-resolved approximation, but
must have the documented shape and finite values.

- [ ] **Step 5.4: Run tests to verify Task 5 passes**

Run package tests. Expected: PASS.

## Task 6: Documentation, README, Examples, And Workflows

**Files:**
- Create: `docs/Project.toml`
- Create: `docs/make.jl`
- Create: `docs/src/index.md`
- Create: `docs/src/getting-started.md`
- Create: `docs/src/manual/conventions.md`
- Create: `docs/src/manual/models.md`
- Create: `docs/src/manual/spectra.md`
- Create: `docs/src/examples/ferromagnetic-chain.md`
- Create: `docs/src/reference/api.md`
- Modify: `README.md`
- Create: `examples/ferromagnetic_chain.jl`
- Create: `.github/workflows/CI.yml`
- Create: `.github/workflows/documentation.yml`
- Create: `test/docs_smoke_tests.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 6.1: Add docs smoke tests**

Create `test/docs_smoke_tests.jl`:

```julia
@testset "docs scaffold" begin
    root = dirname(@__DIR__)
    @test isfile(joinpath(root, "docs", "Project.toml"))
    @test isfile(joinpath(root, "docs", "make.jl"))
    @test isfile(joinpath(root, "docs", "src", "index.md"))
end
```

- [ ] **Step 6.2: Run tests to verify docs smoke test fails**

Expected: FAIL because docs files do not exist.

- [ ] **Step 6.3: Create Documenter site**

Use `Documenter = "1"` in `docs/Project.toml` and a `docs/make.jl` patterned on
`GaussianFermions.jl`, with pages for home, getting started, manual, examples,
and API.

- [ ] **Step 6.4: Rewrite README and add example script**

README must include installation, scope, quick start, output shape conventions,
docs link, and repository map. The quick start must run without plotting.

- [ ] **Step 6.5: Add CI and docs workflows**

Use branch `master` in workflow triggers and Julia versions `1.10` and `1`.

- [ ] **Step 6.6: Run tests and docs build**

Run:

```bash
/Users/ren/.local/bin/julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
/Users/ren/.local/bin/julia --project=docs --startup-file=no -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
/Users/ren/.local/bin/julia --project=docs --startup-file=no docs/make.jl
```

Expected: all commands exit 0.

## Task 7: Cleanup And Final Review

**Files:**
- Delete or untrack: `Manifest.toml`, `.DS_Store`, `src/.DS_Store`
- Move or delete: `demo.jl`, `test/KHzigzag.jl` depending on whether they can be
  turned into useful examples.
- Review: all source, tests, docs, workflows.

- [ ] **Step 7.1: Remove tracked generated/local files**

Run non-destructive git removals:

```bash
git rm --cached Manifest.toml .DS_Store src/.DS_Store
git rm test/KHzigzag.jl
```

Move useful demo content to `examples/` before deleting `demo.jl`.

- [ ] **Step 7.2: Run final verification**

Run:

```bash
/Users/ren/.local/bin/julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
/Users/ren/.local/bin/julia --project=docs --startup-file=no docs/make.jl
git status --short
```

Expected: tests and docs pass; git status contains only intentional rebuild
changes.

- [ ] **Step 7.3: Commit implementation**

Commit with:

```bash
git add .
git commit -m "feat: rebuild SpinWave model and docs"
```
