cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "high_density_stability_v101"
        Pes = collect(4.:0.1:6.)
        T  = 1.0
        max_runs = 5
for ρa in collect(0.45:0.05:0.95), Dθ in [100., 10000.]
        local param
        λ_step = (0.1)*sqrt(Dθ)
        λs = Pes*sqrt(Dθ)
        λmax = (11.)*sqrt(Dθ)
        #
        param = pde_param(; name = name, ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = λ_step, λmax = λmax, λs=λs)
        #
        push!(params,param)
end
params = []
name = "stability_2d_actpass"
pert = "rand"
T  = 1.0
for ρa in [0.4], ρp in [0.,0.2,0.4], Pe in [20.], Dθ in [100.]
        local param
                #
                param = pde_param_1d(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, 
                        save_interval = 0.001, max_steps = 1e7,
                        pert = pert 
                )
                #
                push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
length(params)/nworkers()
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#