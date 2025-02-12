using Lattice

space = Grid{2}((32,32),(8,8))
pars = Phi4_params(0.3, 0.025)
phi4 = Phi4_workspace(space)
ϕ = ScalarField(Float64,space)

heatup!(ϕ)

save(ϕ,"scemo.hdf5")