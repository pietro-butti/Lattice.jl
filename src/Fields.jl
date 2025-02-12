module Fields
    using ..Geometry
    using Random, JLD2

    struct Field{T,D,N} # dof type, dimensions, components
        conf::Array{T,3}
        function Field{T,D,N}(lattice::Grid{D,M,B,F}) where {T,N,D,M,B,F} 
            return new{T,D,N}(
                Array{T,3}(undef,lattice.bsz, N, lattice.rsz)
            )
        end
    end
    export Field

    # ========================= INITIALIZATION =========================
    function heatup!(ϕ::Field{T,D,N}) where {T,D,N}
        rand!(ϕ.conf)
    end
    
    function freeze!(ϕ::Field{T,D,N}) where {T,D,N}
        ϕ.conf .= ones(T,size(ϕ.conf)...)
    end

    export heatup!, freeze!
    # ==================================================================
    
    # =============================== I/O ==============================
    function save(ϕ::Field{T,D,N}, name) where {T,D,N}
        mkpath(dirname(name))
        jldsave("$name.jld2"; ϕ.conf)
        return nothing
    end
    
    
    function read!(ϕ::Field{T,D,N}, from) where {T,D,N}
        buffer = load_object(from)
        @assert eltype(buffer)==T
        @assert size(ϕ.conf)==size(buffer)
        ϕ.conf .= buffer
        return nothing
    end
    export save, read!
    # ==================================================================

end