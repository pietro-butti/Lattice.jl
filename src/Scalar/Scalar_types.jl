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

# =============================== I/O ==============================
function save(ϕ::NScalarField{T,D,N}, name) where {T,D,N}
    jldsave("$name.jld2"; ϕ.conf)
    return nothing
end


function read!(ϕ::NScalarField{T,D,N}, from) where {T,D,N}
    buffer = load_object(from)

    @assert eltype(buffer)==T
    @assert size(ϕ.conf)==size(buffer)

    ϕ.conf .= buffer

    return nothing
end
# ==================================================================