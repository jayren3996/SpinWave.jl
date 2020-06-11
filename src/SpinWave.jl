module SpinWave
using LinearAlgebra
export spinsystem, pushk!, pushω!, kpoints
#--- UV matrix
function UV(θ,ϕ)
    N = length(θ)
    U = Matrix{ComplexF64}(undef, N, 3)
    V = Matrix{Float64}(undef, N, 3)
    cosθ, sinθ, cosϕ, sinϕ = cos.(θ), sin.(θ), cos.(ϕ), sin.(ϕ)
    @. U[:, 1] = cosθ * cosϕ + 1im * sinϕ
    @. U[:, 2] = cosθ * sinϕ - 1im * cosϕ
    @. U[:, 3] = -sinθ
    @. V[:, 1] = sinθ * cosϕ
    @. V[:, 2] = sinθ * sinϕ
    @. V[:, 3] = cosθ
    U, V
end
#--- HBlock
struct HBlock
    BondH::Matrix{ComplexF64} # BondH[4, NJ]
    DiagH::Vector{Float64}    # DiagH[N]
    Ind  ::Matrix{Int64}      # Ind[NJ, 2]
    Vec  ::Matrix{Int64}      # Vec[NJ, d]
end

sandwich(v1,m,v2) = sum(v1 .* (m * v2))
function HBlock(U,V,ind,vec,mat,N,NJ)
    UC = conj.(U)
    BH = Array{ComplexF64}(undef, 4, NJ)
    DH = zeros(Float64,N)
    for i = 1:NJ
        n, m = ind[i,:]
        BH[1,i] = sandwich(UC[n,:], mat[i,:,:],  U[m,:])
        BH[2,i] = sandwich(UC[n,:], mat[i,:,:], UC[m,:])
        BH[3,i] = sandwich( U[n,:], mat[i,:,:],  U[m,:])
        BH[4,i] = sandwich( U[n,:], mat[i,:,:], UC[m,:])
        dhmn = sandwich(V[n,:], mat[i,:,:], V[m,:])
        DH[n] -= dhmn
        DH[m] -= dhmn
    end
    HBlock(BH, DH, ind, vec)
end
#--- H(k)
function hamiltonian!(Hk,k,hblock::HBlock,N,NJ)
    BH  = hblock.BondH
    DH  = hblock.DiagH
    ind = hblock.Ind
    vec = hblock.Vec
    fill!(Hk, 0.0)
    for i = 1:NJ
        n, m = ind[i, :]
        expkv = exp(-1im * dot(k, vec[i,:])) / 2
        H11 = expkv * BH[1,i]
        H12 = expkv * BH[2,i]
        H21 = expkv * BH[3,i]
        H22 = expkv * BH[4,i]
        Hk[n  , m  ] += H11
        Hk[n  , m+N] += H12
        Hk[n+N, m  ] += H21
        Hk[n+N, m+N] += H22
        Hk[m  , n  ] += conj(H11)
        Hk[m  , n+N] += conj(H21)
        Hk[m+N, n  ] += conj(H12)
        Hk[m+N, n+N] += conj(H22)
    end
    for i=1:N Hk[i,i] += DH[i]; Hk[i+N,i+N] += DH[i] end
end
#--- WBlock
struct WBlock
    UU        ::Array{ComplexF64,3}
    Sublattice::Matrix{Float64}
end

function WBlock(U,t,αβ,N,ND)
    UU = Array{ComplexF64}(undef, ND, 2*N, 2*N)
    UC = conj.(U)
    for i=1:ND, j=1:N, k=1:N
        α,β = αβ[i,:]
        UU[i,j  ,k  ] = UC[j,α] * U[k,β]
        UU[i,j  ,k+N] = UC[j,α] * UC[k,β]
        UU[i,j+N,k  ] =  U[j,α] * U[k,β]
        UU[i,j+N,k+N] =  U[j,α] * UC[k,β]
    end
    WBlock(UU, t)
end

function wmat!(Wk,k,wblock::WBlock,N,ND)
    t = wblock.Sublattice
    UU = wblock.UU
    expm = [exp(-1im * dot(k, t[n,:])) for n=1:N]
    expp = conj.(expm)
    for n=1:N, m=1:N
        expkt = expm[n] * expp[m]
        for i=1:ND
            Wk[i,n  ,m  ] = UU[i,n  ,m  ] * expkt
            Wk[i,n  ,m+N] = UU[i,n  ,m+N] * expkt
            Wk[i,n+N,m  ] = UU[i,n+N,m  ] * expkt
            Wk[i,n+N,m+N] = UU[i,n+N,m+N] * expkt
        end
    end
end
#--- Spectrum
function diagonalization(Hk,N)
    C = cholesky(Hk)
    Kd = C.L
    K = C.U
    Im = Diagonal([ones(N); -ones(N)])
    temp = Hermitian(K * Im * Kd)
    D, U = eigen(temp)
    D[1:N] *= -1
    T = inv(K) * U * Diagonal(sqrt.(D))
    D, T
end

function spectrum!(Hk,N)
    Hk[N+1:end,:] *= -1
    real(eigvals(Hk))
end
#--- Green Function
function green!(Gω,ω,D,η,N)
    @. Gω[1:N] = -η / ((ω + D[1:N])^2 + η^2)
    @. Gω[N+1:end] = η / ((ω - D[N+1:end])^2 + η^2)
end
#--- F Matrix
function fmat!(Fk,Wk,T,N,ND)
    for i=1:ND, j=1:2*N
        Fk[i,j] = real(dot(T[:,j], Wk[i,:,:], T[:,j]))
    end
end
#---correlation
correlation(Fk, Gω, N) = sum(Fk * Gω)/(2N)
#--- Spin Wave Solver
mutable struct SpinSystem
    hblock::HBlock
    wblock::WBlock
    Hk    ::Array{ComplexF64,2}
    Wk    ::Array{ComplexF64,3}
    Dk    ::Array{Float64,1}
    Tk    ::Array{ComplexF64,2}
    Gkω   ::Array{Float64,1}
    Fk    ::Array{Float64,2}
    η     ::Float64
    N     ::Int64
    NJ    ::Int64
    ND    ::Int64
end
function spinsystem(θ,ϕ,ind,vec,mat,αβ,t,η)
    N      = length(θ)
    NJ     = size(ind, 1)
    ND     = size(αβ, 1)
    U, V   = UV(θ, ϕ)
    hblock = HBlock(U, V, ind, vec, mat, N, NJ)
    wblock = WBlock(U, t, αβ, N, ND)
    Hk     = Array{ComplexF64}(undef, 2*N, 2*N)
    Wk     = Array{ComplexF64}(undef, ND, 2*N, 2*N)
    Dk     = Array{Float64}(undef, 2*N)
    Tk     = Array{ComplexF64}(undef, 2*N, 2*N)
    Gkω    = Array{Float64}(undef, 2*N)
    Fk     = Array{Float64}(undef, ND, 2*N)
    SpinSystem(hblock, wblock, Hk, Wk, Dk, Tk, Gkω, Fk, η, N, NJ, ND)
end
function pushk!(s::SpinSystem, k)
    hamiltonian!(s.Hk, k, s.hblock, s.N, s.NJ)
    wmat!(s.Wk, k, s.wblock, s.N, s.ND)
    Dk, Tk = diagonalization(s.Hk, s.N)
    s.Dk = Dk
    s.Tk = Tk
    fmat!(s.Fk, s.Wk, s.Tk, s.N, s.ND)
end
function pushω!(s::SpinSystem, ω)
    green!(s.Gkω, ω, s.Dk, s.η, s.N)
    correlation(s.Fk, s.Gkω, s.N)
end
#--- Kpoints
struct KPoints
    N ::Int64
    kp::Matrix{Float64}
end

function kpoints(k,n)
    ln = length(n)
    N  = sum(n) - ln + 1
    d  = size(k,2)
    kp = Array{Float64,2}(undef,N,d)
    ni = 1
    for i = 1:ln
        for j = 1:d
            kp[ni:ni+n[i]-1,j] .= 2*π*range(k[i,j],k[i+1,j],length=n[i])
        end
        ni += n[i]-1
    end
    KPoints(N,kp)
end
end
