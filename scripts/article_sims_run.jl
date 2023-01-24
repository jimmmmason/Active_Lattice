cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###


###
#generic parameters close to stable
params = []
        pert = "rand"
        T  = 4.0
        χ = 1.0
        L = 128
        #using Roots
        #f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        #root = find_zero(f, (0.6,  0.8))
for ρ in [0.7], Pe in [24.], Dθ in [4.]
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
        Pe = 40.
        Dθ = 100.
        name = "article_sim"
param = sim_param_fraction(; name = name, 
        ρ = ρ, Pe = Pe, χ = χ, T = T, 
        Dθ = Dθ, L = 128, Δt = 0.001
)
push!(params,param)
# attempt #2 
# lower Dθ as to not break jump rate 
# for Dθ = 4. , χ = 0.5 
# minimium of stab -> ρ = 0.92 , Pe = 17.546606636771152
# for Dθ = 4. , χ = 0.25 
# minimium of stab -> ρ = 0.946 , Pe = 29.68393713002393
params = []
        pert = "rand"
        T  = 2.0
        L = 128
χ = 0.5
for ρ in [0.92], Pe in [15.,20.,25.], Dθ in [4.]
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
χ = 0.25
for ρ in [0.946], Pe in [25.,30.,35.], Dθ in [4.]
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
#run sims
pmap(run_sim, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_sim_vid, params; distributed = true, batch_size=1, on_error=nothing,)
#
###
# attempt #3
# try extra mixing to encourage waves? 
γ = 0.1
χ = 0.5
params = []
        pert = "rand"
        T  = 4.0
        L = 128
for ρ in [0.8], Pe in [250.], Dθ in [4.]
        name = "article_sim_γ=$(γ))"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, γ= γ, Dθ = Dθ,
                        T = T, L = 128, Δt = 0.001, 
                )
                #
                push!(params,param)
end
#run sims
pmap(run_sim, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_sim_vid, params; distributed = true, batch_size=1, on_error=nothing,)
#
pmap(make_sim_vid_lite, params; distributed = true, batch_size=1, on_error=nothing,)
#
make_sim_vid_lite(param)
###
param = params[2]