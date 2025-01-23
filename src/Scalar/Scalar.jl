module Scalar
    using ..Geometry
    using ..MolecularDynamics
    using TimerOutputs, Random, JLD2

    include("Scalar_types.jl")
        export NScalarField, ScalarField
        export heatup!, freeze!
        export save, read!

    include("Scalar_Phi4.jl")
        export Phi4_params
        export action_krnl!, compute_action!, force_krnl!, compute_force!
        export Phi4_workspace, leapfrog!, HMC!

end

