using Pkg
Pkg.activate("/Users/pietro/code/software/Lattice.jl")

using Revise, TimerOutputs, Random
using Lattice, Printf


##%

space = Grid{2}((32,32),(8,8))

u = GaugeField(U1{Float64},space)
heatup!(u)

