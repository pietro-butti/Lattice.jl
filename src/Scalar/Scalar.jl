module Scalar
    using ..Geometry
    using ..Fields
    using ..MolecularDynamics
    using TimerOutputs, Random, JLD2

    NScalarField(::Type{T}, N, lattice::Grid{D,M,B,F}) where {T,D,M,B,F} = Field{T,D,N}(lattice)
    ScalarField(::Type{T},     lattice::Grid{D,M,B,F}) where {T,D,M,B,F} = Field{T,D,1}(lattice)
    export NScalarField, ScalarField, scemo

    include("Scalar_Phi4.jl")
        export Phi4_params
        export action_krnl!, compute_action!, force_krnl!, compute_force!
        export Phi4_workspace, leapfrog!, HMC!

    include("Scalar_maps.jl")
        export to_pic

end

