"""
    Lattice

Real-space and reciprocal-space lattice bases.

`vectors` is a `3 x 3` matrix whose columns are the real-space basis vectors.
`reciprocal` is the reciprocal basis with the `2π` convention, so
`vectors' * reciprocal == 2πI`.
"""
struct Lattice
    vectors::Matrix{Float64}
    reciprocal::Matrix{Float64}

    function Lattice(vectors::AbstractMatrix{<:Real})
        size(vectors) == (3, 3) || throw(ArgumentError("lattice vectors must be a 3 x 3 matrix"))
        A = Matrix{Float64}(vectors)
        all(isfinite, A) || throw(ArgumentError("lattice vectors must be finite"))
        abs(det(A)) > sqrt(eps(Float64)) || throw(ArgumentError("lattice vectors must be non-singular"))
        B = 2pi * inv(transpose(A))
        return new(A, B)
    end
end

"""
    lattice(a, b, c; angles=(90, 90, 90), degrees=true)
    lattice(lengths; angles=(90, 90, 90), degrees=true)
    lattice(vectors)

Construct a [`Lattice`](@ref). Length/angle constructors use conventional cell
parameters ``a,b,c,α,β,γ`` where `α` is the angle between `b` and `c`, `β`
between `a` and `c`, and `γ` between `a` and `b`.
"""
function lattice(a::Real, b::Real, c::Real; angles=(90, 90, 90), degrees::Bool=true)
    return lattice([a, b, c]; angles=angles, degrees=degrees)
end

function lattice(lengths::AbstractVector{<:Real}; angles=(90, 90, 90), degrees::Bool=true)
    length(lengths) == 3 || throw(ArgumentError("lattice lengths must contain exactly 3 entries"))
    all(isfinite, lengths) || throw(ArgumentError("lattice lengths must be finite"))
    any(<=(0), lengths) && throw(DomainError(lengths, "lattice lengths must be positive"))

    α, β, γ = _angles_radians(angles, degrees)
    sinγ = sin(γ)
    abs(sinγ) > sqrt(eps(Float64)) || throw(ArgumentError("lattice angle γ makes the cell singular"))

    a, b, c = Float64.(lengths)
    cx = c * cos(β)
    cy = c * (cos(α) - cos(β) * cos(γ)) / sinγ
    cz2 = c^2 - cx^2 - cy^2
    cz2 > 100eps(Float64) || throw(ArgumentError("lattice parameters define a singular cell"))

    vectors = [
        a b*cos(γ) cx
        0 b*sinγ cy
        0 0 sqrt(max(cz2, 0.0))
    ]
    return Lattice(vectors)
end

lattice(vectors::AbstractMatrix{<:Real}) = Lattice(vectors)

"""
    reciprocal_cartesian(lat, q)

Convert a reciprocal-lattice-unit vector `q` to Cartesian coordinates using
`lat.reciprocal`.
"""
function reciprocal_cartesian(lat::Lattice, q::AbstractVector{<:Real})
    length(q) == 3 || throw(ArgumentError("q must contain exactly 3 entries"))
    all(isfinite, q) || throw(ArgumentError("q must be finite"))
    return lat.reciprocal * Float64.(q)
end

function _angles_radians(angles, degrees::Bool)
    length(angles) == 3 || throw(ArgumentError("lattice angles must contain exactly 3 entries"))
    all(isfinite, angles) || throw(ArgumentError("lattice angles must be finite"))
    vals = Float64.(collect(angles))
    degrees && (vals .= deg2rad.(vals))
    any(a -> a <= 0 || a >= pi, vals) && throw(DomainError(angles, "lattice angles must lie between 0 and π"))
    return vals
end
