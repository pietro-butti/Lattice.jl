module Scalar
    using ..Lattice

    # ============================= TYPES =============================
    struct NScalarField{T,D,N}
        conf::Array{T,3}
        NScalarField{T,N}(lattice::Grid{D,M,B,F}) where {T,N,D,M,B,F} = new{T,D,N}(
            Array{T,3}(undef,lattice.bsz, N, lattice.rsz)
        )
    end
    NScalarField{T}(n,lattice)     where T = NScalarField{T,n}(lattice)
    ScalarField(::Type{T},lattice) where T = NScalarField{T}(1,lattice)
    # =================================================================
    
    # ========================= INITIALIZATION =========================
    function heatup!(ϕ::NScalarField{T,D,N}) where {T,D,N}
        ϕ.conf .= rand(size(ϕ.conf)...)
    end

    function freeze!(ϕ::NScalarField{T,D,N}) where {T,D,N}
        ϕ.conf .= ones(size(ϕ.conf)...)
    end
    # ==================================================================
    

    # include("Scalar_observables.jl")


    export NScalarField, ScalarField
    export heatup!, freeze!
end

