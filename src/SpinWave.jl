module SpinWave
using LinearAlgebra
export spinsystem, kpoints, spinwave
#--- Types
const Umat{T} = Array{T,2} where T<:Number
const Vmat{T} = Array{T,2} where T<:Real
const Indmat{T} = Array{T,2} where T<:Integer
const Vecmat{T} = Array{T,2} where T<:Integer
const Matmat{T} = Array{T,3} where T<:Number
const tmat{T} = Array{T,2} where T<:Real
const αβmat{T} = Array{T,2} where T<:Integer
struct KPoints{T<:Real}
    N ::Int64
    kp::Matrix{T}
end
#--- Helper functions
sandwich(v1,m,v2) = sum(v1 .* (m * v2))
function UV(θ::AbstractVector,
            ϕ::AbstractVector,
            T1::DataType=ComplexF64,
            T2::DataType=Float64)
    N = length(θ)
    cosθ, sinθ, cosϕ, sinϕ = cos.(θ), sin.(θ), cos.(ϕ), sin.(ϕ)
    U = Array{T1}(undef, N,3)
    V = Array{T2}(undef, N,3)
    @. U[:, 1] = cosθ * cosϕ + 1im * sinϕ
    @. U[:, 2] = cosθ * sinϕ - 1im * cosϕ
    @. U[:, 3] = -sinθ
    @. V[:, 1] = sinθ * cosϕ
    @. V[:, 2] = sinθ * sinϕ
    @. V[:, 3] = cosθ
    U, V
end
function kpoints(k::Matrix,
                 n::Vector;
                 ϵ::Real=0.0)
    ln = length(n)
    N  = sum(n) - ln + 1
    d  = size(k,2)
    kp = Array{Float64,2}(undef,N,d)
    ni = 1
    for i = 1:ln
        for j = 1:d
            kp[ni:ni+n[i]-1,j] .= 2π * range(k[i,j],k[i+1,j],length=n[i])
        end
        ni += n[i]-1
    end
    if ϵ != 0.0
        kp .+= ϵ
    end
    KPoints(N, kp)
end
#--- include
include("HBLock.jl")
include("WBlock.jl")
include("Diagonalization.jl")
include("Solver.jl")

end
