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
    eigenvectors::Matrix{ComplexF64}
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

    real_values = real.(eig.values)
    positives = sort!(collect(v for v in real_values if v > tol))
    zeros_needed = N - length(positives)
    if zeros_needed < 0
        positives = positives[end-N+1:end]
        zeros_needed = 0
    end
    energies = vcat(zeros(Float64, zeros_needed), positives)
    length(energies) == N || throw(ArgumentError("could not identify $N nonnegative bosonic modes"))

    return DiagonalizationResult(energies, ComplexF64.(eig.values), ComplexF64.(eig.vectors))
end
