"""
    QPath

Expanded reciprocal-space scan path.

`q_rlu` is a `3 x nq` matrix of reciprocal-lattice-unit coordinates.
`q_cartesian` is either `nothing` or a `3 x nq` Cartesian matrix when the path
is constructed with a [`Lattice`](@ref). `ticks` contains the column indices of
the original vertices after endpoint de-duplication.
"""
struct QPath
    q_rlu::Matrix{Float64}
    q_cartesian::Union{Nothing,Matrix{Float64}}
    ticks::Vector{Int}
    labels::Vector{String}
end

"""
    qpath(vertices; points, labels=nothing)

Expand piecewise-linear path vertices in reciprocal lattice units.

Each segment includes its starting point and endpoint, but shared interior
vertices are stored once. Therefore segment counts `[3, 4]` produce
`3 + 4 - 1 == 6` total q-points.
"""
function qpath(vertices; points, labels=nothing)
    verts = _qvertices(vertices)
    length(verts) >= 2 || throw(ArgumentError("qpath requires at least two vertices"))
    nseg = length(verts) - 1
    counts = _segment_counts(points, nseg)
    labs = _path_labels(labels, length(verts))

    columns = Vector{Vector{Float64}}()
    ticks = Int[1]
    for iseg in 1:nseg
        n = counts[iseg]
        start = verts[iseg]
        stop = verts[iseg + 1]
        ts = range(0.0, 1.0; length=n)
        first_index = iseg == 1 ? 1 : 2
        for i in first_index:n
            t = ts[i]
            push!(columns, @. (1 - t) * start + t * stop)
        end
        push!(ticks, length(columns))
    end

    q_rlu = reduce(hcat, columns)
    return QPath(q_rlu, nothing, ticks, labs)
end

"""
    qpath(vertices, lat::Lattice; points, labels=nothing)

Expand a q-path and also store Cartesian q-vectors using `lat.reciprocal`.
"""
function qpath(vertices, lat::Lattice; points, labels=nothing)
    path = qpath(vertices; points=points, labels=labels)
    q_cartesian = lat.reciprocal * path.q_rlu
    return QPath(path.q_rlu, q_cartesian, path.ticks, path.labels)
end

Base.length(path::QPath) = size(path.q_rlu, 2)

function _qvertices(vertices)
    verts = [Float64.(collect(v)) for v in vertices]
    for v in verts
        length(v) == 3 || throw(ArgumentError("each qpath vertex must contain exactly 3 entries"))
        all(isfinite, v) || throw(ArgumentError("qpath vertices must be finite"))
    end
    return verts
end

function _segment_counts(points::Integer, nseg::Integer)
    points >= 2 || throw(ArgumentError("each qpath segment must contain at least 2 points"))
    return fill(Int(points), nseg)
end

function _segment_counts(points, nseg::Integer)
    counts = collect(points)
    length(counts) == nseg || throw(ArgumentError("points must have one entry per qpath segment"))
    all(n -> n isa Integer, counts) || throw(ArgumentError("qpath point counts must be integers"))
    all(n -> n >= 2, counts) || throw(ArgumentError("each qpath segment must contain at least 2 points"))
    return Int.(counts)
end

function _path_labels(labels, nverts::Integer)
    if labels === nothing
        return fill("", nverts)
    end
    labs = string.(collect(labels))
    length(labs) == nverts || throw(ArgumentError("labels must have one entry per qpath vertex"))
    return labs
end
