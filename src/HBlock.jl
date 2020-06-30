#--- UV matrix
function UV!(U::Umat,
             V::Vmat,
             θ::AbstractVector,
             ϕ::AbstractVector)
    cosθ, sinθ, cosϕ, sinϕ = cos.(θ), sin.(θ), cos.(ϕ), sin.(ϕ)
    @. U[:, 1] = cosθ * cosϕ + 1im * sinϕ
    @. U[:, 2] = cosθ * sinϕ - 1im * cosϕ
    @. U[:, 3] = -sinθ
    @. V[:, 1] = sinθ * cosϕ
    @. V[:, 2] = sinθ * sinϕ
    @. V[:, 3] = cosθ
end
