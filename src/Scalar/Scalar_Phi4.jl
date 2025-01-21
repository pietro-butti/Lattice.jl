###                               
### author:  pietro.butti
### file:    Scalar_Phi4.jl
### created: Tue 21 Jan 2023 16:03:00 CEST
### 

struct Phi4_params{T}
    κ::T
    λ::T
end

# DEFINE ACTION=================================================
    function action_krnl!(action::NScalarField{T,D,1}, ϕ::NScalarField{T,D,1}, coord::NTuple{2}, params::Phi4_params{T}, lattice::Grid{D,M,B,F}) where {T,D,B,F,M}
        # Fetch local field component
        (b,r) = coord
        ϕₙ  = ϕ.conf[b,1,r]
        ϕₙ² = ϕₙ*ϕₙ
        
        # Compute local action
        Sₙ =  ϕₙ² + params.λ * (ϕₙ² - one(T))*(ϕₙ² - one(T))
        staple = zero(T)
        for μ in 1:D
            b_up, r_up = up(coord, μ, lattice)
            staple += ϕ.conf[b_up, 1, r_up]
        end
        Sₙ += -(params.κ+params.κ) * ϕₙ * staple 

        # Update local value
        action.conf[b,1,r] = Sₙ
        
        return nothing
    end

    function compute_action!(act::NScalarField{T,D,1}, ϕ::NScalarField{T,D,1}, params::Phi4_params{T}, lattice::Grid{D,M,B,F}) where {T,D,B,F,M}
        @timeit "ϕ⁴ action computation" begin
            for b in 1:lattice.bsz
                for r in 1:lattice.rsz
                    action_krnl!(act,ϕ,(b,r),params,lattice)
                end
            end
        end
        return  nothing
    end
# DEFINE FORCE =================================================
    function force_krnl!(force::NScalarField{T,D,1}, ϕ::NScalarField{T,D,1}, coord::NTuple{2}, params::Phi4_params{T}, lattice::Grid{D,M,B,F}) where {T,D,B,F,M}
        # Fetch local field component
        (b,r) = coord
        ϕₙ  = ϕ.conf[b,1,r]
        ϕₙ² = ϕₙ*ϕₙ
        
        # Compute local action
        Fₙ = convert(T,2)*ϕₙ + convert(T,4)* params.λ * (ϕₙ² - one(T)) * ϕₙ
        staple = zero(T)
        for μ in 1:D
            b_up, r_up, b_dw, r_dw = updw(coord, μ, lattice)
            staple += ϕ.conf[b_up, 1, r_up] + ϕ.conf[b_dw, 1, r_dw]  
        end
        Fₙ += -convert(T,2)*params.κ * staple 

        # Update local value
        force.conf[b,1,r] = Fₙ
        
        return nothing
    end

    function compute_force!(fc::NScalarField{T,D,1}, ϕ::NScalarField{T,D,1}, params::Phi4_params{T}, lattice::Grid{D,M,B,F}) where {T,D,B,F,M}
        @timeit "ϕ⁴ force computation" begin
            for b in 1:lattice.bsz
                for r in 1:lattice.rsz
                    force_krnl!(fc,ϕ,(b,r),params,lattice)
                end
            end
        end
        return nothing
    end
# ==============================================================

# DEFINE WORKSPACE =============================================
    struct Phi4_workspace{T,D}
        _phi::NScalarField{T,D,1}
        _scalar1::NScalarField{T,D,1}
        _scalar2::NScalarField{T,D,1}
        mom::NScalarField{T,D,1}

        compute_action!::Function
        compute_force!::Function
        
        function Phi4_workspace(::Type{T},lattice::Grid{D,M,B,F}) where {T,D,M,B,F}
            _phi  = ScalarField(T,lattice)
            _bin1 = ScalarField(T,lattice)
            _bin2 = ScalarField(T,lattice)
            _mom =  ScalarField(T,lattice)
            _S(act,phi,params) = compute_action!(act,phi,params,lattice)
            _F(frc,phi,params) = compute_force!(frc,phi,params,lattice) 
            return new{T,D}(_phi,_bin1,_bin2,_mom,_S,_F)
        end
        Phi4_workspace(lattice::Grid) = Phi4_workspace(Float64,lattice::Grid)
        
    end
# ==============================================================

# DEFINE HMC FUNCTIONS =========================================
    function leapfrog!(ϕ::NScalarField{T,D,1}, params::Phi4_params{T}, int::IntrScheme{NI,T}, ws::Phi4_workspace{T,D}) where {T,D,NI}
        @timeit "MD evolution" begin
            for leap in 1:int.Nleaps
                # update conf
                ϕ.conf .= ϕ.conf .+ int.eps/convert(T,2) .* ws.mom.conf
                
                # update mom
                ws.compute_force!(ws._scalar1, ϕ, params)
                ws.mom.conf .= ws.mom.conf .- int.eps .* ws._scalar1.conf
                
                # update conf
                ϕ.conf .= ϕ.conf .+ int.eps/convert(T,2) .* ws.mom.conf
                println("cra cra") 
            end
        end
        return nothing
    end

    function HMC!(ϕ::NScalarField{T,D,1}, params::Phi4_params{T}, int::IntrScheme{NI,T}, ws::Phi4_workspace{T,D}) where {T,D,NI}
        @timeit "HMC trajectory" begin
            # copy the configuration
            ws._phi.conf .= ϕ.conf 

            # generate momenta
            randn!(ws.mom.conf)

            # compute hamiltonian (force used as bin)
            ws.compute_action!(ws._scalar1,ϕ,params)
            H = sum(ws._scalar1.conf) + convert(T,0.5)*sum(ws.mom.conf.^2)
            
            # molecular dynamics 
            leapfrog!(ϕ,params,int,ws)
            
            # compute dH 
            ws.compute_action!(ws._scalar1,ϕ,params)
            dH = sum(ws._scalar1.conf) + convert(T,0.5)*sum(ws.mom.conf.^2)
            dH -= H

            # acc/rej
            acc = true
            pacc = exp(-dH)
            if (pacc < 1.0)
                if (pacc < rand()) 
                    ϕ.conf .= ws._phi.conf # reject
                    acc = false
                end
            end

            return dH, acc

        end
    end
# ==============================================================


