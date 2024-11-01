module Lattice

    include("Geometry.jl")
        using .Geometry
        export Grid
        export up, dw, updw, point_coord
        export BC_PERIODIC
        
    include("Ising.jl")
        using .IsingModel
        export SpinField, behold
        export heatup!, freeze!
        export staple_sum, localE, magnetization
        export accept, update_ranflips!

    include("Scalar/Scalar.jl")
        using .Scalar
        export NScalarField, ScalarField, heatup!, freeze!

end # module Lattice
