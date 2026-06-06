# Conventions

## Lattice Vectors

[`Lattice`](@ref) stores real-space basis vectors as columns of a `3 x 3`
matrix. Reciprocal vectors include the `2π` convention:

```julia
lat.vectors' * lat.reciprocal ≈ 2π * I
```

Fractional site positions and bond offsets are expressed in lattice units.

## Q-Points

[`qpath`](@ref) accepts vertices in reciprocal lattice units. A segment point
count includes both endpoints, and shared vertices are stored once.

For example, two segments with counts `[3, 4]` produce `3 + 4 - 1 == 6`
q-points. The [`QPath`](@ref) object stores:

- `q_rlu[component, iq]`,
- `q_cartesian[component, iq]` when constructed as `qpath(vertices, lat; ...)`,
  otherwise `nothing`,
- `ticks`, the vertex positions in the expanded path,
- `labels`, one label per input vertex.

## Spectrum Arrays

[`SpinWaveSpectrum`](@ref) uses stable dimensions:

- `energies[mode, iq]`,
- `correlations[a, b, mode, iq]`.

The broadened [`EnergyGrid`](@ref) returned by [`broaden`](@ref) uses:

- `intensity[iω, iq]`.

These shapes are chosen so raw modes are indexed by mode first, while plotted
energy grids are indexed by energy first.
