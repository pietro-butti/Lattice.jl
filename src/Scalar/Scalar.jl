module Scalar
    using ..Geometry
    using ..MolecularDynamics
    using TimerOutputs, Random

    # ============================= TYPES =============================
    struct NScalarField{T,D,N} # {precision, dimensions, components}
        conf::Array{T,3}
        NScalarField{T,N}(lattice::Grid{D,M,B,F}) where {T,N,D,M,B,F} = new{T,D,N}(
            Array{T,3}(undef,lattice.bsz, N, lattice.rsz)
        )
    end
    NScalarField{T}(n,lattice) where T = NScalarField{T,n}(lattice)
    ScalarField(::Type{T},lattice::Grid{D,M,B,F}) where {T,D,M,B,F} = NScalarField{T}(1,lattice)
    # =================================================================
    
    # ========================= INITIALIZATION =========================
    function heatup!(ϕ::NScalarField{T,D,N}) where {T,D,N}
        rand!(ϕ.conf)
    end

    function freeze!(ϕ::NScalarField{T,D,N}) where {T,D,N}
        ϕ.conf .= ones(size(ϕ.conf)...)
    end
    # ==================================================================
    
    export NScalarField, ScalarField
    export heatup!, freeze!

    include("Scalar_Phi4.jl")
        export Phi4_params
        export action_krnl!, compute_action!, force_krnl!, compute_force!
        export Phi4_workspace, leapfrog!, HMC!




end

