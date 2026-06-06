"""
    UnstableHamiltonianError

Raised when a quadratic bosonic Hamiltonian is Hermitian but not positive
semidefinite within the requested tolerance.
"""
struct UnstableHamiltonianError <: Exception
    min_eigenvalue::Float64
    tolerance::Float64
end

function Base.showerror(io::IO, err::UnstableHamiltonianError)
    print(
        io,
        "bosonic Hamiltonian is unstable: minimum Hermitian eigenvalue ",
        err.min_eigenvalue,
        " is below tolerance ",
        err.tolerance,
    )
end

"""
    DiagonalizationResult

Result of dense bosonic dynamical-matrix diagonalization.
"""
struct DiagonalizationResult
    energies::Vector{Float64}
    eigenvalues::Vector{ComplexF64}
    modes::Matrix{ComplexF64}
    metric_residual::Float64
end

"""
    diagonalize_bosonic(H; tol=1e-9)

Validate and diagonalize a dense quadratic bosonic Hamiltonian.

The Hamiltonian must be Hermitian and positive semidefinite within `tol`. The
returned energies are the nonnegative positive-norm branch of the eigenvalues of
`ηH`, where `η = Diagonal([ones(N); -ones(N)])`.
"""
function diagonalize_bosonic(H::AbstractMatrix; tol::Real=1e-9)
    size(H, 1) == size(H, 2) || throw(ArgumentError("bosonic Hamiltonian must be square"))
    iseven(size(H, 1)) || throw(ArgumentError("bosonic Hamiltonian dimension must be even"))
    all(isfinite, H) || throw(ArgumentError("bosonic Hamiltonian entries must be finite"))

    Hc = Matrix{ComplexF64}(H)
    hermitian_residual = norm(Hc - Hc')
    hermitian_residual <= tol || throw(ArgumentError("bosonic Hamiltonian must be Hermitian; residual=$hermitian_residual"))

    hvals = eigvals(Hermitian(Hc))
    mineig = minimum(real.(hvals))
    mineig >= -tol || throw(UnstableHamiltonianError(mineig, Float64(tol)))

    N = size(Hc, 1) ÷ 2
    metric = Diagonal(vcat(ones(Float64, N), -ones(Float64, N)))
    dyn = Matrix(metric * Hc)
    eig = eigen(dyn)
    max_imag = maximum(abs, imag.(eig.values); init=0.0)
    max_imag <= sqrt(tol) || throw(ArgumentError("bosonic dynamical matrix has complex modes; max imaginary part=$max_imag"))

    candidates = _positive_metric_modes(real.(eig.values), eig.vectors, metric, Float64(tol))
    length(candidates) >= N || throw(ArgumentError("could not identify $N positive-norm bosonic modes"))

    energies = zeros(Float64, N)
    modes = Matrix{ComplexF64}(undef, 2N, N)
    for mode in 1:N
        energy, vector = candidates[mode]
        energies[mode] = energy
        modes[:, mode] .= vector
    end

    metric_residual = norm(modes' * metric * modes - I(N))
    metric_residual <= sqrt(tol) || throw(ArgumentError("positive modes are not paraunitary-normalized; residual=$metric_residual"))

    return DiagonalizationResult(energies, ComplexF64.(eig.values), modes, metric_residual)
end

function _positive_metric_modes(values::AbstractVector{<:Real}, vectors::AbstractMatrix{<:Complex}, metric, tol::Float64)
    cluster_tol = sqrt(tol)
    order = sortperm(values)
    candidates = Tuple{Float64,Vector{ComplexF64}}[]
    cursor = firstindex(order)
    while cursor <= lastindex(order)
        first_value = values[order[cursor]]
        last = cursor
        while last < lastindex(order)
            next_value = values[order[last + 1]]
            abs(next_value - first_value) <= cluster_tol * max(1.0, abs(first_value)) || break
            last += 1
        end

        if first_value >= -tol
            group = order[cursor:last]
            V = vectors[:, group]
            gram = Hermitian(Matrix(V' * metric * V))
            signature = eigen(gram)
            for idx in eachindex(signature.values)
                norm = real(signature.values[idx])
                norm > tol || continue
                vector = ComplexF64.(V * signature.vectors[:, idx]) ./ sqrt(norm)
                push!(candidates, (max(first_value, 0.0), vector))
            end
        end
        cursor = last + 1
    end
    return candidates
end
