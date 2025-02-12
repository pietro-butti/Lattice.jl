module GaugeFields
    using ..Geometry
    using ..Fields
    using ..MolecularDynamics
    using ..Groups
    using TimerOutputs, Random, JLD2
    
    GaugeField(::Type{G},lattice::Grid{D,M,B,F}) where {G,D,M,B,F} = Field{G,D,D}(lattice)
    export GaugeField

    include("GaugeFields_actions.jl")
        export gauge_params
        export plaquette_krnl!, compute_plaquette! 

    include("GaugeFields_hmc.jl")
        export gauge_workspace
        export force_wilson_krnl!, compute_force_wilson!

end