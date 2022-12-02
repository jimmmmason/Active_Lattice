cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
#include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "stability_1d_rand"
pert = "rand"
T  = 1.0
for ρa in collect(0.4:0.05:0.95), Pe in collect(1.:0.5:10.), Dθ in [100.]
        local param
        #
        param = pde_param_1d(; name = name, 
                ρa = ρa, Pe = Pe, ρp = 0., T = T, 
                Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, 
                save_interval = 0.01, max_steps = 1e7,
                pert = pert 
        )
        #
        push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
length(params)/nworkers()
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
#