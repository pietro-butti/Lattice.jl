
module IsingModel
    using ..Lattice

    # ============================= TYPES =============================
    struct SpinField{T,N}
        conf::Array{T,2}
        SpinField{T}(space::Grid{N,M,B,D}) where {T,N,M,B,D} = new{T,N}(Array{T,2}(undef, space.bsz, space.rsz))
    end

    function behold(spin::SpinField{T,2}, space::Grid{2,M,B,D}) where {T,M,B,D}
        tmp = Matrix{Float64}(undef, space.blk .* space.rbk)
        
        for idx in CartesianIndices(spin.conf)
            coord = point_coord(Tuple(idx),space)
            tmp[coord] = spin.conf[idx]
        end
        
        return tmp
    end
    # ==================================================================

    # ========================= INITIALIZATION =========================  
    function heatup!(spin::SpinField{T,N}, space::Grid{N,M,B,D}) where {T,N,M,B,D}
        spin.conf .= rand([-1,1], space.bsz, space.rsz)
        return
    end
    function freeze!(spin::SpinField{T,N}, space::Grid{N,M,B,D}) where {T,N,M,B,D}
        spin.conf .= ones(T, space.bsz, space.rsz)
        return
    end
    # ===================================================================

    # =========================== OBSERVABLES ===========================  
    function staple_sum(spin::SpinField{T,N}, space::Grid{N,M,B,D}, coord::NTuple{2}) where {T,N,M,B,D}
        staple = 0.
        for μ in 1:N
            b_up, r_up, b_dw, r_dw = updw(coord, μ, space)
            staple = staple + 
                spin.conf[b_up,r_up] + 
                spin.conf[b_dw,r_dw]
        end    
        return staple
    end
    localE(spin,space,coord) = - spin.conf[coord...] * staple_sum(spin,space,coord)

    magnetization(spin::SpinField{T,N}, space::Grid{N,M,B,D}) where {T,N,M,B,D} = sum(spin.conf)/(space.bsz*space.rsz)
    # ===================================================================


    # =========================== MONTE-CARLO ===========================
    @inline accept(dS) = dS<=0 ? true : rand()<exp(-dS)


    function update_ranflips!(spin::SpinField{T,N}, copy::SpinField{T,N}, space::Grid{N,M,B,D}, temp, Nflips) where {T,N,M,B,D}
        # Make a buffer copy
        copy.conf[:] .= spin.conf[:]
        
        # Update conf
        dE = 0.
        for _ in 1:Nflips
            ib = rand(1:space.bsz)
            ir = rand(1:space.rsz)
            
            spin.conf[ib,ir] = -spin.conf[ib,ir]
            dE += 2. * localE(spin,space,(ib,ir))
        end

        # Metropolis acc/rej step
        acc = accept(dE/temp)
        if !acc
            spin.conf[:] .= copy.conf[:]
        end
        
        return acc, dE/temp
    end
    # ===================================================================

    export SpinField, behold
    export randomize!, freeze!
    export staple_sum, localE, magnetization
    export accept, update_ranflips!

end