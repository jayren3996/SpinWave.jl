"""
    LocalFrame

Right-handed orthonormal frame attached to an ordered spin moment.
"""
struct LocalFrame
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

"""
    local_frame(moment)

Construct a stable right-handed local frame whose local z-axis follows
`moment`.
"""
function local_frame(moment::AbstractVector{<:Real})
    z = _finite_vector3(moment, "spin moment")
    nz = norm(z)
    nz > sqrt(eps(Float64)) || throw(DomainError(moment, "spin moment must be nonzero"))
    z ./= nz

    reference = abs(z[3]) < 0.9 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    x = cross(reference, z)
    x ./= norm(x)
    y = cross(z, x)
    y ./= norm(y)
    return LocalFrame(x, y, z)
end

_transverse(frame::LocalFrame) = complex.(frame.x, frame.y)
