using Pkg
Pkg.activate("/Users/pietro/code/software/Lattice.jl")

using Revise, TimerOutputs, Random
using Lattice, Printf


##%

space = Grid{2}((32,32),(8,8))


ϕ = ScalarField(Float64,space)
heatup!(ϕ)

pic = Array{Float64,2}(undef,space.iL...)
for b in 1:space.bsz
    for r in 1:space.rsz
        I = point_coord((b,r),space)
        pic[I] = ϕ.conf[b,1,r]        
    end
end