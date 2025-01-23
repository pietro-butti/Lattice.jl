using Pkg
Pkg.activate("/Users/pietro/code/software/Lattice.jl")

using Revise, TimerOutputs, Random
using Lattice, Printf


##%

# space = Grid{4}((32,32,32,32),(8,8,8,8))
space = Grid{2}((32,32),(8,8))

pars = Phi4_params(0.19, 1.3282)
phi4 = Phi4_workspace(space)

ϕ = ScalarField(Float64,space)
bin = ScalarField(Float64,space)
freeze!(ϕ)

int = omf2(Float64,0.01,10)


accumulator = 0
for _ in 1:10000

    dH, acc = HMC!(ϕ,pars,int,phi4)
    global accumulator += acc

    mag = abs(sum(ϕ.conf))/prod(space.iL)
    
    out = @sprintf("dH=%6.3f (%s)     %10.8f", dH, acc ? "yes" : "noo", mag)

    println(out)
end

print_timer()
