include("../src/SpinWave.jl")
using .SpinWave
using LinearAlgebra
using PyPlot
#--- K-H-Γ-G-J2-J3 Model
include("KH4.jl")
function spins(ck,ga,cg,c2,c3)
    xx = [ck 0 0;0 cg ga;0 ga cg]
    yy = [cg 0 ga;0 ck 0;ga 0 cg]
    zz = [cg ga 0;ga cg 0;0 0 ck]
    j2 = Array(Diagonal([1.0,1.0,1.0]))*c2
    j3 = Array(Diagonal([1.0,1.0,1.0]))*c3
    θ = begin
        y = 2 * sqrt(2) * ga
        x = -2 * ck + 2*cg + ga
        θi = atan(y/x)/2 + π/2
        θm = π - θi
        [θi, θm, θm, θi]
    end
    ϕ = begin
        ϕi = π/4
        ϕm = -3*π/4
        [ϕi, ϕm, ϕm, ϕi]
    end
    model(xx,yy,zz,j2,j3,θ,ϕ)
end

ωmax = 10.0
ωl = range(0,ωmax,length=1000)
rx,ry = 0,0

kp1 = kpoints([ 1.0  0.0;0 0; 1.0 1.0],[173,200]);kl1 = kp1.kp
kp2 = kpoints([ 0.5  1.5;0 0; 0.0 2.0],[173,200]);kl2 = kp2.kp
kp3 = kpoints([-0.5  1.5;0 0;-1.0 1.0],[173,200]);kl3 = kp3.kp

K,Γ,G,J₂,J₃ = -1.5, 0.05, 0.01, 1.0, 0.5
s = spins(K,Γ,G,J₂,J₃)
C = zeros(length(ωl),kp1.N)
ek = zeros(kp1.N,4)
#--- Calculate Spin Wave
corr, spec = solve(s, kl1, ωl)
#--- Plot
cmax = 3.0
@. corr[corr>cmax] = cmax
pygui(true)
figure(dpi=200)
contourf(1:kp1.N,ωl,corr',levels=50,cmap="jet")
plot(1:kp1.N,spec,"--",lw=0.7,color="grey")
title("K=$K,Γ=$Γ,G=$G,J₂=$J₂,J₃=$J₃")
colorbar()
show()
