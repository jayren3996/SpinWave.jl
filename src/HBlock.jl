# Struct
struct HBlock{T1<:Number, T2<:Real, Tind, Tvec}
    N    ::Int64
    NJ   ::Int64
    Ind  ::Indmat{Tind}     # Ind[NJ, 2]
    Vec  ::Vecmat{Tvec}     # Vec[NJ, d]
    BondH::Matrix{T1}       # BondH[4, NJ]
    DiagH::Vector{T2}       # DiagH[N]
end
#--- HBlock
function hblock(U::Umat,
                V::Vmat,
                ind::Indmat,
                vec::Vecmat,
                mat::Matmat;
                T1::DataType=ComplexF64,
                T2::DataType=Float64)
    N = size(U,1)
    NJ = size(ind,1)
    UC = conj.(U)
    BH = Matrix{T1}(undef, 4,NJ)
    DH = Vector{T2}(undef, NJ)
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
    HBlock(N,NJ,ind,vec,BH,DH)
end
#--- hamiltonian
function hamiltonian!(Hk::Matrix,
                      hb::HBlock,
                      k::Real)
    BH  = hblock.BondH
    DH  = hblock.DiagH
    ind = hblock.Ind
    vec = hblock.Vec
    fill!(Hk, 0.0)
    for i = 1:hb.NJ
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
    for i = 1:hb.N
        Hk[i,i] += DH[i]
        Hk[i+N,i+N] += DH[i]
    end
end
