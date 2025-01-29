struct YM_workspace{G,A,T,D}
    _u::Field{G,D,D}
    _bin1::Field{G,D,D}
    _bin2::Field{G,D,D}
    _mom::Field{A,D,D}

    compute_action!::Function
    compute_force!::Function
end




####### ONLY DEFINED FOR U(1)
function force_YMN1_krnl!(
    acc::Field{A,D,D}, 
    U::Field{G,D,D}, 
    coord::NTuple{2}, 
    dir::Int,
    lattice::Grid{D,M,B,F}
    ) where {G,A,D,B,F,M}

    # Get local coordinates
    (b,r) = coord

    # Initialize staple
    staple = zero(G)

    b_pμ,r_pμ = up((b,r),dir,lattice)
    for ν in 1:D
        if ν!=dir
            b_pν,r_pν = up((b,r),ν,lattice)
            staple += U.conf[b_pμ,ν,r_pμ] / (U.conf[b,ν,r] * U.conf[b_pν,dir,r_pν])
        end
    end

    # for ν in 1:D
    #     if ν!=dir
    #         b_pν_mμ,r_pν_mμ = up((b_pν,r_pν),ν  ,lattice)
    #         b_mμ   ,r_mμ    = up((b   ,r   ),dir,lattice)
    #         staple += dag(U.conf[b_mμ,dir,r_mμ] * U.conf[b_pν_mμ,ν,r_pν_mμ]) * U.conf[b_mμ,ν ,r_mμ]
    #     end
    # end


    return nothing
end

function compute_force!(acc::Field{T,D,1}, U::Field{G,D,D}, lattice::Grid{D,M,B,F}) where {G,T,D,B,F,M}
    @timeit "Plaquette computation" begin
        for b in 1:lattice.bsz
            for r in 1:lattice.rsz
                action_krnl!(acc,U,(b,r),lattice)
            end
        end
    end
    return  nothing    
end



