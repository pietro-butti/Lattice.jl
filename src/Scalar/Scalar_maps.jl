function to_pic(ϕ::Field{T,D,1}, lattice::Grid{D,M,B,F}) where {T,D,M,B,F}
    pic = Array{Float64,D}(undef, lattice.iL...)
    for b in 1:space.bsz
        for r in 1:space.rsz
            I = point_coord((b,r),space)
            pic[I] = ϕ.conf[b,1,r]        
        end
    end    
    return pic
end