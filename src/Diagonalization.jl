#--- Spectrum
function diagonalization!(D::Vector,
                          T::Matrix,
                          Hk::Matrix)
    N = length(D) ÷ 2
    C = cholesky!(Hk)
    Kd = C.L
    K = C.U
    Im = Diagonal([ones(N); -ones(N)])
    temp = Hermitian(K * Im * Kd)
    vals, vecs = eigen(temp)
    D[1:N] .= -1 * vals[1:N]
    D[N+1:2N] .= vals[N+1:2N]
    T .= inv(K) * vecs * Diagonal(sqrt.(D))
end

function spectrum!(Hk::Matrix, N::Int64)
    Hk[N+1:end,:] *= -1
    real(eigvals!(Hk))
end
#--- Green Function
function green!(Gω::Vector,
                ω::Real,
                D::Vector,
                η::Real,
                N::Int64)
    @. Gω[1:N] = -η / ((ω + D[1:N])^2 + η^2)
    @. Gω[N+1:end] = η / ((ω - D[N+1:end])^2 + η^2)
end
#--- F Matrix
function fmat!(Fk::Matrix,
               Wk::Array,
               Tk::Matrix,
               N::Int64,
               ND::Int64)
    for i=1:ND, j=1:2*N
        Fk[i,j] = real(dot(Tk[:,j], Wk[i,:,:], Tk[:,j]))
    end
end
#---correlation
function correlation(Fk::Matrix,
                     Gω::Vector,
                     N::Int64)
    sum(Fk * Gω)/(2N)
end
