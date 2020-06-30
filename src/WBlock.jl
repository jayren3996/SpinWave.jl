# Struct
struct WBlock{T1<:NUmber,T2<:Real}
    N   ::Int64
    ND  ::Int64
    UU  ::Array{T1,3}
    SubL::Matrix{T2}
end
#--- Wblock
function wblock(U::Umat,
                t::tmat,
                αβ::αβmat;
                T1::DataType=ComplexF64)
    N = size(t,1)
    ND = size(αβ,1)
    UU = Array{T1}(undef, ND, 2*N, 2*N)
    UC = conj.(U)
    for i=1:ND, j=1:N, k=1:N
        α,β = αβ[i,:]
        UU[i,j  ,k  ] = UC[j,α] * U[k,β]
        UU[i,j  ,k+N] = UC[j,α] * UC[k,β]
        UU[i,j+N,k  ] =  U[j,α] * U[k,β]
        UU[i,j+N,k+N] =  U[j,α] * UC[k,β]
    end
    WBlock(N, NJ, UU, t)
end
#--- wmat
function wmat!(Wk::Array, wb::WBlock, k)
    t = wb.SubL
    UU = wb.UU
    N = wb.N
    ND = wb.ND
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
