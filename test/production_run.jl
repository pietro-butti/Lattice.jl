# using Pkg
# Pkg.activate("/Users/pietro/code/software/Lattice.jl")

using Lattice
using Printf, TimerOutputs, JLD2, ArgParse

function parse_commandline()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--size", "-L"
            nargs = '*'
        "--blocks", "-b"
            nargs = '*'
        "--kappa", "-k"
            arg_type = Float64
        "--lambda", "-l"
            arg_type = Float64
        "--eps"
            arg_type = Float64
            default  = 0.01
        "--Nleaps"
            arg_type = Int
            default  = 10
        "--Ntraj"
            arg_type = Int
            default  = 10000
        "--start_from"
            default = "cold"
        "--save_each"
            default = 123456
        "--save_folder"
            default = "./"
        "--save_as"
            default = "___"
    end

    return parse_args(args)
end

function main()
    pargs = parse_commandline()

    # Geometry parameters
    size = Tuple(parse(Int,x) for x in pargs["size"])
    blks = Tuple(parse(Int,x) for x in pargs["blocks"])
    dims = length(size)
    
    # Theory parameters    
    λ = pargs["lambda"]
    κ = pargs["kappa"]
    
    # HMC parameters
    start_from = pargs["start_from"]
    eps    = pargs["eps"] 
    Nleaps = pargs["Nleaps"] 
    Ntraj  = pargs["Ntraj"] 

    # Output parameters
    save_each = pargs["save_each"]
    save_as   = pargs["save_as"]

    # =============================================
    # =============================================
    ensemble_code = "./d$(dims)L$(size[1])kappa$(κ)lambda$(λ)"

    # Set geomtery and workspace
    space = Grid{dims}(size,blks)
    pars = Phi4_params(κ, λ)
    phi4 = Phi4_workspace(space)

    # initialize field
    ϕ = ScalarField(Float64,space)
    if start_from=="cold"
        freeze!(ϕ)
    elseif start_from=="hot"
        heatup!(ϕ)
    else
        read!(ϕ,start_from)
    end
    
    int = omf2(Float64,eps,Nleaps) # hmc integrator
    
    # =============================================
    counter = 0
    for mctime in 1:Ntraj
        # HMC update
        dH, acc = HMC!(ϕ,pars,int,phi4)
        counter += acc

        # Measurements
        mag  = sum(ϕ.conf)/prod(space.iL)
        mag2 = sum(ϕ.conf.^2)/prod(space.iL)
        
        # I/O
        if save_each!=123456
            if save_each%mctime==0
                prefix = save_as=="___" ? ensemble_code : save_as
                conf_name = "$prefix.cfg$(lpad(mctime,5,'0'))"
                save(ϕ,joinpath([save_folder,conf_name]))
            end
        end

        # Output
        accs =  acc ? "yes" : "no " 
        out = @sprintf("hmc: %4i %10.6f (%s)     %10.8f      %10.8f",mctime,dH,accs,mag,mag2)
        println(out)
    end

    print_timer()
    println("acc/rej = $(counter/Ntraj)")
    println(ensemble_code)
end

main()