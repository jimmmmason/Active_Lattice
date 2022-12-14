cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###


###
#varying parameters
params = []
        pert = "n=1"
        T  = 1.0
        δ = 0.01
for ρ in collect(0.4:0.01:0.95), χ in [0.25,0.5,0.75], Pe in collect(5.:5.:100.), Dθ in [100.]
        local param
        #
        name = "article_stability_1d_δ=$(δ)"
        param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
                )
        #
        push!(params,param)
end
#run pdes
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
###


###
#travelling wave parameters
params = []
        pert = "n=1"
        T  = 1.0
        χ = 0.25
        using Roots
        f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        root = find_zero(f, (0.6,  0.8))
for ρ in (-collect(0.0:0.01:0.05).+root), Pe in [20.], Dθ in [100.], δ in [1e-6]
        Nx = 128
        Nθ = 64
        name = "article_wave_1d_δ=$(δ)"
        local param
                #
                param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40
                )
                #
                push!(params,param)
end
#run pdes
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)
###


###
#generic 2d wave parameters
params = []
        pert = "rand"
        T  = 2.0
        χ = 0.25
        using Roots
        f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        root = find_zero(f, (0.6,  0.8))
for ρ in (-collect(0.0:0.01:0.05).+root), Pe in [20.], Dθ in [100.], δ in [1e-2]
        Nx = 50
        Nθ = 20
        name = "article_rand_2d_δ=$(δ)"
        local param
                #
                param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40
                )
                #
                push!(params,param)
end
#run pdes
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)
###

###
#generic 2d unstable parameters
pert = "rand"
        T  = 2.0
        χ  = 0.75
        ρ  = 0.8
        Pe = 20. 
        Dθ = 100. 
        δ  = 1e-2
        Nx = 50
        Nθ = 20
        name = "article_rand_2d_δ=$(δ)"
param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40
        )
        params = [param]
#run pdes
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)
###