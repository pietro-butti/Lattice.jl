module Lattice
    using TimerOutputs

    include("Geometry.jl")
        using .Geometry
        export Grid
        export up, dw, updw, point_coord, point_index
        export BC_PERIODIC

    include("Fields.jl")
        using .Fields 
        export Field, heatup!, freeze!, save, read

    # include("Groups/Groups.jl")
    #     using .Groups
    #     export Group, Algebra
    #     export U1, U1alg
    #     export dot, expm, exp, dag, unitarize, inverse, tr, projalg, norm, norm2, isgroup, alg2mat, dev_one
    
    include("HMC.jl")
        using .MolecularDynamics
        export IntrScheme, omf4, leapfrog, omf2
            
    include("Scalar/Scalar.jl")
        using .Scalar
        export NScalarField, ScalarField, scemo
        export Phi4_params, action_krnl!, compute_action!, force_krnl!, compute_force!, Phi4_workspace, leapfrog!, HMC!
        export to_pic

    # include("GaugeFields/GaugeFields.jl")
    #     using .GaugeFields
    #     export GaugeField, heatup!
    #     export gauge_params, plaquette_krnl!, compute_plaquette!
    #     export gauge_workspace, force_wilson_krnl!, compute_force_wilson!
        
    # include("Ising.jl")
    #     using .IsingModel
    #     export SpinField, behold
    #     export heatup!, freeze!
    #     export staple_sum, localE, magnetization
    #     export accept, update_ranflips!


end # module Lattice
