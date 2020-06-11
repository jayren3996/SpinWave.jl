# SpinWave.jl
Julia package for spin-wave calculation.

Compute spin-wave spectrum and spin-spin correlation in momentum space (imaginary part).

## Installation

Run the following script in the ```Pkg REPL``` :

```julia
pkg> add https://github.com/jayren3996/SpinWave.jl
```

## Examples 

Example on 2D anti-ferromagnetic Heissenberg model. The input is the following
```julia
include("./SpinWave.jl")
using .SpinWave
using LinearAlgebra
using PyPlot
# 2D anti-ferromagnetic Heisenberg Model
ind = [1 2; 2 1; 3 4; 4 3; 1 2; 3 1; 2 4; 4 2] # interaction index.
vec = [0 0; 1 0; 0 0; 1 0; 0 0; 0 1; 0 0; 0 1] # distance between two interacting sites.
H0 = Diagonal([1.0, 1.0, 1.0])
mat = Array{Float64}(undef, 8, 3, 3) # interaction term.
for i=1:8 mat[i,:,:] .= H0 end
θ = [0, π, π, 0] # spin direction.
ϕ = [0, 0, 0, 0] # spin direction.
t = [0 0; 0.5 0; 0 0.5; 0.5 0.5] # sub-lattice index.
αβ = [1 1;2 2;3 3] # spin summation indeces.
η = 0.1 # quasi-particle lifetime.
s = spinsystem(θ,ϕ,ind,vec,mat,αβ,t,η) # cunstruct spin system
```
The computation part is
```julia
ωmax = 5.0
ωl = range(0,ωmax,length=200)
ε = 1e-5 # avoid zero energy mode
kp = kpoints([ε ε; 1+ε ε; 1+ε 1+ε; ε ε],[100, 100, 141]); kl = kp.kp

C = zeros(length(ωl),kp.N)
ek = zeros(kp.N,4)
# Calculate spin wave
for i=1:kp.N
    pushk!(s,kl[i,:])
    ek[i,:] .= s.Dk[5:8]
    for j=1:length(ωl)
        C[j,i] += pushω!(s,ωl[j])
    end
end
# Plot
cmax = 3.0
@. C[C>cmax] = cmax
pygui(true)
contourf(1:kp.N,ωl,C,levels=50,cmap="jet")
plot(1:kp.N, ek, "--", lw=0.7, color="grey")
colorbar()
show()
```
