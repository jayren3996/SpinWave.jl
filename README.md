# SpinWave.jl
Julia package for spin-wave calculation.

Compute spin-wave spectrum and spin-spin correlation in momentum space (imaginary part).

## Installation

Run the following script in the ```Pkg REPL``` :

```julia
pkg> add https://github.com/jayren3996/SpinWave.jl
```

## Algorithm

![Algorithm](https://raw.github.com/jayren3996/SpinWave.jl/master/Algorithm.jpg)

## Examples 

Example on 2D anti-ferromagnetic Heissenberg model. The input is the following
```julia
using SpinWave
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
θ = [0, π, π, 0] # spin direction.
ϕ = [0, 0, 0, 0] # spin direction.
t = [0 0; 0.5 0; 0 0.5; 0.5 0.5] # sub-lattice index.
αβ = [1 1;2 2;3 3] # spin summation indeces.
η = 0.1 # quasi-particle lifetime.

s = spinsystem(θ,ϕ,ind,vec,mat,t,αβ,η) # cunstruct spin system
```
The computation part is
```julia
#--- spin wave calculation
ωs = range(0, 5, length=200)
kp = kpoints([0 0;1 0;1 1;0 0],[100, 100, 141], ϵ=1e-5) # avoid zero energy mode

corr, spec = spinwave(s, kp, ωs)
```

Plot the result

```julia
# Plot
cmax = 3.0
@. corr[corr>cmax] = cmax
pygui(true)
contourf(1:kp.N, ωs, corr, levels=50, cmap="jet")
plot(1:kp.N, spec, "--", lw=0.7, color="grey")
colorbar()
show()
```



