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

    candidates = Tuple{Float64,Int,Float64}[]
    for (idx, value) in pairs(real.(eig.values))
        value >= -tol || continue
        vec = eig.vectors[:, idx]
        metric_norm = real(dot(vec, metric * vec))
        metric_norm > tol || continue
        push!(candidates, (max(value, 0.0), idx, metric_norm))
    end
    sort!(candidates; by=first)
    length(candidates) >= N || throw(ArgumentError("could not identify $N positive-norm bosonic modes"))

    energies = zeros(Float64, N)
    modes = Matrix{ComplexF64}(undef, 2N, N)
    for mode in 1:N
        energy, idx, metric_norm = candidates[mode]
        energies[mode] = energy
        modes[:, mode] .= eig.vectors[:, idx] ./ sqrt(metric_norm)
    end

    metric_residual = norm(modes' * metric * modes - I(N))
    metric_residual <= sqrt(tol) || throw(ArgumentError("positive modes are not paraunitary-normalized; residual=$metric_residual"))

    return DiagonalizationResult(energies, ComplexF64.(eig.values), modes, metric_residual)
end
