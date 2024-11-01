using Lattice 

space = Grid{2}((32,32),(8,8))

field = NScalarField{Float64}(2,space)
bin = NScalarField{Float64}(2,space)
acc = ScalarField(Float64,space)
heatup!(field)


function local_interaction_kernel(
    phi::NScalarField{T,D,N}, 
    coupling, 
    order::Int64, 
    coords::NTuple{2}
    ) where {T,N,D}
    
    b,r = coords

    # Compute inline modulus square
    bin = zero(T)
    for i in 1:N
        bin += phi.conf[b,i,r]*phi.conf[b,i,r]
    end

    # Compute inline nth power
    for _ in 1:order
        bin *= bin
    end

    return coupling * bin
end

function local_interaction!(
    acc::NScalarField{T,D,1}, 
    phi::NScalarField{T,D,N}, 
    coupling, 
    order::Int64, 
    lattice::Grid{D,M,B,F}
    ) where {T,N,D,M,B,F}

    # Cycle over all local elements
    for b in 1:lattice.bsz
        for r in 1:lattice.rsz
            acc.conf[b,1,r] = local_interaction_kernel(phi,coupling,order,(b,r))
        end
    end

    return nothing
end

local_interaction2(phi, coupling) = sum(phi.conf.^2) * coupling


function staple_sum!(
    staple::NScalarField{T,D,N}, 
    phi::NScalarField{T,D,N}, 
    lattice::Grid{D,M,B,F},
    ) where {T,N,D,M,B,F}

    # Cyclec over all spacetime points
    for b in 1:lattice.bsz
        for r in 1:lattice.rsz
            # Cycle over directions
            for μ in 1:D
                b_up, r_up, b_dw, r_dw = updw(coord, μ, space)
                # Cycle over vector dof
                for i in 1:N
                    staple.conf[b,i,r] = phi.conf[b_up,i,r_up] + phi.conf[b_dw,i,r_dw]
                end
            end    
        end
    end

    return staple
end










local_interaction!(acc,field,0.5,2,space)





##

# spin = SpinField{Float64}(space)
# copy = SpinField{Float64}(space)
# freeze!(spin,space)


# obs = []
# counter = 0
# @showprogress for _ in 1:1_000_000
#     acc,dS = update_ranflips!(spin,copy,space,temperature,1)
#     push!(obs, magnetization(spin,space))
    
#     counter += acc
# end

# ##
# plot(obs,title="$(counter/1_000_000)")