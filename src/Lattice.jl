module Lattice
    using TimerOutputs

    include("Geometry.jl")
        using .Geometry
        export Grid
        export up, dw, updw, point_coord
        export BC_PERIODIC

    include("HMC.jl")
        using .MolecularDynamics
        export IntrScheme, omf4, leapfrog, omf2
            
    include("Scalar/Scalar.jl")
        using .Scalar
        export NScalarField, ScalarField, heatup!, freeze!
        export Phi4_params, action_krnl!, compute_action!, force_krnl!, compute_force!, Phi4_workspace, leapfrog!, HMC!

    # include("Ising.jl")
    #     using .IsingModel
    #     export SpinField, behold
    #     export heatup!, freeze!
    #     export staple_sum, localE, magnetization
    #     export accept, update_ranflips!


end # module Lattice
