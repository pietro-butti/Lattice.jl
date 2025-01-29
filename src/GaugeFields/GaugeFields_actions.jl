struct gauge_params{T}
    β::T
end

"""
    accumulates ∑[μ<ν] Re Tr[ Uᵤ(n) Uᵥ(n+μ) Uᵤ(n+ν)' Uᵥ(n)' ]
"""
function plaquette_krnl!(
    acc::Field{T,D,1}, 
    U::Field{G,D,D}, 
    coord::NTuple{2}, 
    lattice::Grid{D,M,B,F}
    ) where {G,T,D,B,F,M}

    # Get local coordinates
    (b,r) = coord

    plaq = zero(T)
    for μ in 1:D
        for ν in (μ+1):D
            # Get hop coordinates
            b_pμ,r_pμ = up((b,r),μ,lattice)
            b_pν,r_pν = up((b,r),ν,lattice)

            # build plaquette
            plaq += real(tr( # for U(1) Re(g) = complex(re(g),0)
                (U.conf[b,μ,r] * U.conf[b_pμ,μ,r_pμ]) / 
                (U.conf[b,ν,r] * U.conf[b_pν,ν,r_pν])
            ))
        end
    end
    acc.conf[b,1,r] = plaq

    return nothing
end

function compute_plaquette!(acc::Field{T,D,1}, U::Field{G,D,D}, lattice::Grid{D,M,B,F}) where {G,T,D,B,F,M}
    @timeit "Plaquette computation" begin
        for b in 1:lattice.bsz
            for r in 1:lattice.rsz
                plaquette_krnl!(acc,U,(b,r),lattice)
            end
        end
    end
    return  nothing    
end



