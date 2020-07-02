include("src/SpinWave.jl")
using .SpinWave
using LinearAlgebra
using PyPlot
#--- 2D anti-ferromagnetic Heisenberg Model
ind = [1 2; 2 1; 3 4; 4 3; 1 2; 3 1; 2 4; 4 2] # interaction index.
vec = [0 0; 1 0; 0 0; 1 0; 0 0; 0 1; 0 0; 0 1] # distance between two interacting sites.
mat = begin                                    # interaction term.
    H0 = Diagonal([1.0, 1.0, 1.0])
    m = Array{Float64}(undef, 8, 3, 3)
    for i=1:8 m[i,:,:] .= H0 end
    m
end

θ = [0, π, π, 0]                               # spin direction.
ϕ = [0, 0, 0, 0]                               # spin direction.
t = [0 0; 0.5 0; 0 0.5; 0.5 0.5]               # sub-lattice index.

s = spinsystem(θ,ϕ,ind,vec,mat,t)              # cunstruct spin system
#--- spin wave
ωs = range(0, 5, length=200)
kp = kpoints([0 0;1 0;1 1;0 0],[100, 100, 141])

corr, spec = spinwave(s, kp, ωs)
#--- Plot
@. corr[corr>3.0] = 3.0
pygui(true)
figure(dpi=200)
contourf(1:kp.N,ωs,corr',levels=50,cmap="jet")
plot(1:kp.N,spec,"--",lw=0.7,color="grey")
colorbar()
show()
