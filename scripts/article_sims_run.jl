cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###


###
#generic parameters close to stable
params = []
        pert = "rand"
        T  = 2.0
        χ = 0.25
        L = 128
        using Roots
        f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        root = find_zero(f, (0.6,  0.8))
for ρ in (-collect(0.0:0.01:0.05).+root), Pe in [20.], Dθ in [100.]
        name = "article_sim"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = 128, Δt = 0.001
                )
                #
                push!(params,param)
end
# 
T  = 2.0
        χ  = 0.75
        ρ  = 0.8
        Pe = 20.
        Dθ = 100.
        name = "article_sim"
param = sim_param_fraction(; name = name, 
        ρ = ρ, Pe = Pe, χ = χ, T = T, 
        Dθ = Dθ, L = 128, Δt = 0.001
)
push!(params,param)
#run sims
pmap(run_sim, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_sim_vid, params; distributed = true, batch_size=1, on_error=nothing,)
#
###