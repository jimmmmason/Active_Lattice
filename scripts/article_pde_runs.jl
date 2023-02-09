cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###


###
#varying parameters
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
#
for ρ in [0.55], χ in [1.0], Pe in [12., 14., 16.], Dθ in [4.]
                local param
                #
                name = "article_stability_active_5_1d_δ=$(δ)"
                param = pde_param_fraction(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                                save_interval = 0.01, max_steps = 1e7,
                                pert = pert, k =40, δ = δ,
                        )
                #
                push!(params,param)
end
#for ρ in append!(collect(0.425:0.025:0.9),collect(0.95:0.001:0.999)), χ in [1.0], Pe in collect(2.:2.0:40.), Dθ in [4.]
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
#
name = "article_stability_10_1d_δ=$(δ)"
for ρ in collect(0.45:0.01:0.56), χ in [1.0], Pe in collect(10.:2.0:40.), Dθ in [4.]
        local param
        #
        param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
                )
        #
        push!(params,param)
end
for ρ in collect(0.57:0.01:0.9), χ in [1.0], Pe in collect(6.:2.0:14.), Dθ in [4.]
        local param
        #
        param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
                )
        #
        push!(params,param)
end
#
for ρ in collect(0.91:0.01:0.99), χ in [1.0], Pe in collect(6.:2.0:20.), Dθ in [4.]
        local param
        #
        param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
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
        T  = 2.0
        χ = 0.5
        γ = 0.1
        using Roots
        #f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        #root = find_zero(f, (0.6,  0.8))
        # ρ in (-collect(0.0:0.01:0.05).+root)
for ρ in [0.9], Pe in [150.,125.,200.], Dθ in [4.], δ in [1e-2]
        Nx = 50
        Nθ = 20
        name = "article_wave_1d_δ=$(δ)_γ=$(γ)"
        local param
                #
                param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40, γ = 0.1
                )
                #
                push!(params,param)
end
#run pdes
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)

param = params[2]
perturb_pde_run_1d(param)
###


###
#generic 2d wave parameters
params = []
        pert = "rand"
        T  = 4.0
        χ = 1.0
for ρ in [0.7], Pe in [8., 10., 12., 20.], Dθ in [4.], δ in [1e-2]
        Nx = 128
        Nθ = 64
        name = "article_rand_actually_2d_δ=$(δ)"
        local param
                #
                param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-2, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40
                )
                #
                push!(params,param)
end
#run pdes
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#



pmap(load_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)
###



###
#generic 2d unstable parameters
pert = "rand"
        T  = 4.0
        χ  = 0.75
        ρ  = 0.8
        Pe = 40. 
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
pmap(load_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#make video
pmap(make_video, params; distributed = true, batch_size=1, on_error=nothing,)
#make_video_1d(params[1]; frames = 1000)
###