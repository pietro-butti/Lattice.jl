struct gauge_workspace{G,A,T,D}
    _u::Field{G,D,D}
    _bin1::Field{G,D,D}
    _bin2::Field{G,D,D}
    _mom::Field{A,D,D}
    _force::Field{A,D,D}

    compute_action!::Function
    compute_force!::Function

    function gauge_workspace(::Type{U1{T}}, lattice::Grid{D,M,B,F}) where {T,D,M,B,F}
        new{U1{T},U1alg{T},T,D}(
            
        )
    end
end


################################################################################################################
############################################## ONLY VALID FOR U(1) #############################################
################################################################################################################
function force_wilson_krnl!(
    acc::Field{U1alg{T},D,D}, 
    U::Field{U1{T},D,D}, 
    coord::NTuple{2}, 
    μ::Int,
    lattice::Grid{D,M,B,F}
    ) where {T,D,B,F,M}

    # Get local coordinates
    (b,r) = coord
    loc = zero(U1alg{T})

    b_pμ,r_pμ = up((b,r),μ,lattice)
    for ν in 1:D
        if ν!=μ
            # upper staple
            b_pν,r_pν = up((b,r),ν,lattice)
            loc += projalg((U.conf[b,μ,r]*U.conf[b_pμ,ν,r_pμ]) / (U.conf[b,ν,r] * U.conf[b_pν,μ,r_pν]))
            
            # lower staple
            b_mν,r_mν = dw((b,r),ν,lattice)
            b_mν_pμ,r_mν_pμ = dw((b_mν,r_mν),μ,lattice)
            loc += projalg( adj(U.conf[b_mν,μ,r_mν] * U.conf[b_mν_pμ,μ,r_mν_pμ]) / U.conf[b_mν,ν,r_mν] )
        end
    end
    acc.conf[b,μ,r] = loc

    return nothing
end

function compute_force_wilson!(acc::Field{T,D,1}, U::Field{G,D,D}, lattice::Grid{D,M,B,F}) where {G,T,D,B,F,M}
    @timeit "Force [wilson action] computation" begin
        for b in 1:lattice.bsz
            for r in 1:lattice.rsz
                action_krnl!(acc,U,(b,r),lattice)
            end
        end
    end
    return  nothing    
end



