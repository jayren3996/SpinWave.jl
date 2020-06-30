# Struct
struct SpinSystem{T1<:HBlock,T2<:WBlock, Thk,Twk,Tdk,Ttk,Tgω,Tfk,Tη}
    hblock::T1
    wblock::T2
    Hk    ::Array{Thk,2}
    Wk    ::Array{Twk,3}
    Dk    ::Array{Tdk,1}
    Tk    ::Array{Ttk,2}
    Gkω   ::Array{Tgω,1}
    Fk    ::Array{Tfk,2}
    η     ::Tη
    N     ::Int64
    NJ    ::Int64
    ND    ::Int64
end
#--- Initiation
function spinsystem(θ::Vector,
                    ϕ::Vector,
                    ind::Indmat,
                    vec::Vecmat,
                    mat::Matmat,
                    t::tmat,
                    αβ::αβmat=[1 1;2 2;3 3],
                    η::Real=0.2;
                    Thk::DataType=ComplexF64,
                    Twk::DataType=ComplexF64,
                    Tdk::DataType=Float64,
                    Ttk::DataType=ComplexF64,
                    Tgω::DataType=Float64,
                    Tfk::DataType=Float64)
    N    = length(θ)
    NJ   = size(ind, 1)
    ND   = size(αβ, 1)
    U, V = UV(θ, ϕ)
    hb   = hblock(U, V, ind, vec, mat)
    wb   = WBlock(U, t, αβ)
    Hk   = Array{Thk}(undef, 2*N, 2*N)
    Wk   = Array{Twk}(undef, ND, 2*N, 2*N)
    Dk   = Array{Tdk}(undef, 2*N)
    Tk   = Array{Ttk}(undef, 2*N, 2*N)
    Gkω  = Array{Tgω}(undef, 2*N)
    Fk   = Array{Tfk}(undef, ND, 2*N)
    SpinSystem(hb, wb, Hk, Wk, Dk, Tk, Gkω, Fk, η, N, NJ, ND)
end
#--- push k
function pushk!(s::SpinSystem, k::Real)
    hamiltonian!(s.Hk, s.hblock, k)
    wmat!(s.Wk, s.wblock, k)
    diagonalization!(s.Dk, s.Tk, s.Hk)
    fmat!(s.Fk, s.Wk, s.Tk, s.N, s.ND)
end
#--- push ω
function pushω!(s::SpinSystem, ω::Real)
    green!(s.Gkω, ω, s.Dk, s.η, s.N)
    correlation(s.Fk, s.Gkω, s.N)
end
#--- solve all
function solve(s::SpinSystem, ks::Vector, ωs::Vector)
    nk = length(ks)
    nω = length(ωs)
    corr = Array{Float64}(undef, nk,nω)
    spec = Array{Float64}(undef, nk,s.N)
    for i=1:nk
        pushk!(s,ks[i])
        spec[i,:] .= s.Dk[1:s.N]
        for j=1:nω
            corr[i,j] = pushω!(s,ωs[i])
        end
    end
end
