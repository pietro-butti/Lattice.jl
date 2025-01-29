module Groups
    using Random

    import Base.:*, Base.:+, Base.:-,Base.:/,Base.:\,Base.exp,Base.one,Base.zero
    import Random.rand
    
    abstract type Group end
    abstract type Algebra end
    export Group, Algebra

    include("Groups_U1.jl")
        export U1, U1alg
        export dot, expm, exp, dag, unitarize, inverse, tr, projalg, norm, norm2, isgroup, alg2mat, dev_one

end