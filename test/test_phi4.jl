using Pkg, Revise, UnicodePlots, TimerOutputs, Random
Pkg.activate("/Users/pietro/code/software/Lattice.jl")

using Lattice





##%

space = Grid{4}((32,32,32,32),(8,8,8,8))

pars = Phi4_params(0.19, 1.3282)
phi4 = Phi4_workspace(space)

ϕ = ScalarField(Float64,space)
bin = ScalarField(Float64,space)
freeze!(ϕ)

int = leapfrog(Float64,0.05,20)
HMC!(ϕ,pars,int,phi4)


print_timer()
