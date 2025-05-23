function conf_to_pic(ϕ::Field{T,D,1}, lattice::Grid{D,M,B,F}) where {T,D,M,B,F}
    pic = Array{Float64,D}(undef, lattice.iL...)
    for b in 1:lattice.bsz
        for r in 1:lattice.rsz
            I = point_coord((b,r),lattice)
            pic[I] = ϕ.conf[b,1,r]        
        end
    end    
    return pic
end

function pic_to_conf!(ϕ::Field{T,D,1}, pic::AbstractArray{T,D}, lattice::Grid{D,M,B,F}) where {T,D,M,B,F}
    @assert size(pic)==lattice.iL
    for I in CartesianIndices(axes(pic))
        (b,r) = point_index(I, lattice)
        ϕ.conf[b,1,r] = pic[I]
    end
    return
end
pic_to_conf(pic, lattice) = begin
    phi = ScalarField(eltype(pic),lattice)
    pic_to_conf!(phi,pic,lattice)
    phi
end

