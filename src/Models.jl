"""
    SpinSite

One magnetic site in the crystallographic basis.
"""
struct SpinSite
    label::Symbol
    position::Vector{Float64}
    spin::Float64
    moment::Vector{Float64}
end

"""
    ExchangeMatrix

Reusable `3 x 3` bilinear spin interaction matrix.
"""
struct ExchangeMatrix
    matrix::Matrix{Float64}
end

"""
    Bond

Directed bond from `source` to `target` with integer unit-cell `offset`.
The Hermitian counterpart is added during Hamiltonian assembly.
"""
struct Bond
    matrix::Symbol
    source::Int
    target::Int
    offset::Vector{Int}
    scale::Float64
end

"""
    SpinModel(lat)

Mutable builder for a commensurate bilinear spin-wave model.
"""
mutable struct SpinModel
    lattice::Lattice
    sites::Vector{SpinSite}
    matrices::Dict{Symbol,ExchangeMatrix}
    bonds::Vector{Bond}
end

SpinModel(lat::Lattice) = SpinModel(lat, SpinSite[], Dict{Symbol,ExchangeMatrix}(), Bond[])

"""
    addsite!(model, label, position; spin, moment)

Add a magnetic basis site to `model`.
"""
function addsite!(
    model::SpinModel,
    label::Symbol,
    position::AbstractVector{<:Real};
    spin::Real,
    moment::AbstractVector{<:Real},
)
    _has_site(model, label) && throw(ArgumentError("site label $label already exists"))
    pos = _finite_vector3(position, "site position")
    m = _finite_vector3(moment, "site moment")
    isfinite(spin) || throw(ArgumentError("spin length must be finite"))
    spin > 0 || throw(DomainError(spin, "spin length must be positive"))
    norm(m) > sqrt(eps(Float64)) || throw(DomainError(moment, "site moment must be nonzero"))
    push!(model.sites, SpinSite(label, pos, Float64(spin), m))
    return model
end

"""
    heisenberg(J)

Create an isotropic exchange matrix `J * I`.
"""
function heisenberg(J::Real)
    isfinite(J) || throw(ArgumentError("Heisenberg exchange must be finite"))
    return ExchangeMatrix(Float64(J) * Matrix{Float64}(I, 3, 3))
end

"""
    exchange_matrix(M)

Create a full `3 x 3` exchange matrix.
"""
function exchange_matrix(M::AbstractMatrix{<:Real})
    size(M) == (3, 3) || throw(ArgumentError("exchange matrix must be 3 x 3"))
    A = Matrix{Float64}(M)
    all(isfinite, A) || throw(ArgumentError("exchange matrix entries must be finite"))
    return ExchangeMatrix(A)
end

"""
    dm(D)

Create the antisymmetric exchange matrix for Dzyaloshinskii-Moriya vector `D`.
"""
function dm(D::AbstractVector{<:Real})
    d = _finite_vector3(D, "DM vector")
    Dx, Dy, Dz = d
    return ExchangeMatrix([
        0.0 -Dz Dy
        Dz 0.0 -Dx
        -Dy Dx 0.0
    ])
end

"""
    addmatrix!(model, label, matrix)

Register a named exchange matrix.
"""
function addmatrix!(model::SpinModel, label::Symbol, matrix::ExchangeMatrix)
    haskey(model.matrices, label) && throw(ArgumentError("exchange matrix label $label already exists"))
    model.matrices[label] = matrix
    return model
end

"""
    addbond!(model, matrix, source, target, offset; scale=1)

Add one directed physical bond. `source` and `target` may be site labels or
1-based site indices.
"""
function addbond!(
    model::SpinModel,
    matrix::Symbol,
    source,
    target,
    offset::AbstractVector;
    scale::Real=1,
)
    haskey(model.matrices, matrix) || throw(ArgumentError("unknown exchange matrix label $matrix"))
    src = _site_index(model, source)
    dst = _site_index(model, target)
    off = _integer_vector3(offset, "bond offset")
    isfinite(scale) || throw(ArgumentError("bond scale must be finite"))
    push!(model.bonds, Bond(matrix, src, dst, off, Float64(scale)))
    return model
end

"""
    validate(model)

Check cross-object model invariants and return `model`.
"""
function validate(model::SpinModel)
    isempty(model.sites) && throw(ArgumentError("model must contain at least one magnetic site"))
    labels = [site.label for site in model.sites]
    length(unique(labels)) == length(labels) || throw(ArgumentError("site labels must be unique"))
    for bond in model.bonds
        checkbounds(model.sites, bond.source)
        checkbounds(model.sites, bond.target)
        haskey(model.matrices, bond.matrix) || throw(ArgumentError("bond references unknown exchange matrix $(bond.matrix)"))
    end
    return model
end

function _finite_vector3(v::AbstractVector{<:Real}, name::AbstractString)
    length(v) == 3 || throw(ArgumentError("$name must contain exactly 3 entries"))
    out = Float64.(collect(v))
    all(isfinite, out) || throw(ArgumentError("$name must be finite"))
    return out
end

function _integer_vector3(v::AbstractVector, name::AbstractString)
    length(v) == 3 || throw(ArgumentError("$name must contain exactly 3 entries"))
    all(x -> x isa Integer, v) || throw(ArgumentError("$name must contain integers"))
    return Int.(collect(v))
end

_has_site(model::SpinModel, label::Symbol) = any(site -> site.label == label, model.sites)

function _site_index(model::SpinModel, label::Symbol)
    idx = findfirst(site -> site.label == label, model.sites)
    idx === nothing && throw(ArgumentError("unknown site label $label"))
    return idx
end

function _site_index(model::SpinModel, index::Integer)
    1 <= index <= length(model.sites) || throw(ArgumentError("site index $index is out of range"))
    return Int(index)
end

_site_index(::SpinModel, site) = throw(ArgumentError("site reference must be a Symbol label or Integer index, got $(typeof(site))"))
