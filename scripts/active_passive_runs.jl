cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###


### 1d pde wave runs
params = []
        pert = "n=1"
        T  = 24.0
        δ  = 1e-2
        k = 20
        save_interval = 0.01
        χ = 0.1
        Dθ = 4.0
        Nx = 50
        Nθ = 20
name = "ap_wave_1d_δ=$(δ)_χ=$(χ)_l=$(1/sqrt(Dθ))"
#
for ρ in 0.8:0.01:0.99, χ in [χ], Pe in 40:1:50, Dθ in [Dθ]
                local param
                #
                param = pde_param_fraction(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                                save_interval = save_interval, max_steps = 1e7,
                                pert = pert, k =k, δ = δ,
                        )
                #
                push!(params,param)
end
#
#run pdes
pmap(load_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
###
#make videos
pmap(pmap_make_phase_video_1d , params; distributed = true, batch_size=1, on_error=nothing,)
#