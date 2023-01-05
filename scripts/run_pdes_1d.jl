cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "stability_1d_actpass_2"
pert = "n=1"
T  = 1.0
for ρ in collect(0.4:0.05:0.6), ρp in collect(0.2:0.2:0.4), Pe in collect(50.:5.:100.), Dθ in [100.]
        local param
        ρa = ρ - ρp
        if ρa>0.
                #
                param = pde_param_1d(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert 
                )
                #
                push!(params,param)
        end
end
for ρ in collect(0.4:0.05:0.95), ρp in collect(0.0:0.1:0.5), Pe in collect(5.:5.:100.), Dθ in [1.]
        local param
        ρa = ρ - ρp
        if ρa>0.
                #
                param = pde_param_1d(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert 
                )
                #
                push!(params,param)
        end
end
params = []
name = "stability_1d_actpass_3"
pert = "rand"
T  = 2.0
for ρa in [0.225], ρp in [0.5], Pe in [377.1], Dθ in [100.]
        local param
                #
                param = pde_param(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-6, Nx = 50, Nθ = 20, 
                        save_interval = 1e-4, max_steps = 1e10,
                        pert = pert, δ = 0.001
                )
                #
                push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
length(params)/nworkers()
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#


params = []
pert = "n=1"
T  = 1.0
ρp = 0.4
using Roots
f(x) = ap_lin_stab_line(x,ρp; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
root = find_zero(f, (0.1,  0.3))
for ρa in [root], ρp in [ρp], Pe in [20.], Dθ in [100.], δ in [1e-6],k in [40]
        Nx = 128
        Nθ = 64
        name = "stability_1d_actpass_3_δ=$(δ)"
        local param
                #
                param = pde_param_k(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k=k
                )
                #
                push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
###
#length(params)/nworkers()
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#
###
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
pmap(make_video_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#
make_video_1d(params[1]; frames = 1000)