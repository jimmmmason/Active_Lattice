cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "high_density_stability_Dθ=100"
λs = 20.:10.:250.
max_runs = 10
λ_step = 20.
λmax = 250.
Dθ = 100. 
T  = 1.0
for ρa in collect(0.05:0.05:0.95)
        local param
        #
        param = pde_param(; name = name, λ = λ , ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = λ_step, λmax = λmax, λs=λs)
        #
        push!(params,param)
end
name = "high_density_stability_Dθ=100"
λs = 20.:100.:2000.
max_runs = 10
λ_step = 200.
λmax = 2040.
Dθ = 10000. 
T  = 1.0
λ = 0.
for ρa in collect(0.05:0.05:0.95)
        local param
        #
        param = pde_param(; name = name, λ = λ , ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = λ_step, λmax = λmax, λs=λs)
        #
        push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
length(params)/nworkers()
pmap(run_stab_search_stupid_mode, params; distributed = true, batch_size=1, on_error=nothing,)
