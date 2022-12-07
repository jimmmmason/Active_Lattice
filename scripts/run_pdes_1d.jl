cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
#include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "stability_1d_actpass_2"
pert = "n=1"
T  = 1.0
for ρ in collect(0.4:0.05:0.95), ρp in collect(0.1:0.1:0.5), Pe in collect(5.:5.:50.), Dθ in [100.]
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
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
length(params)/nworkers()
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#