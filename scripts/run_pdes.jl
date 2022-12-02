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
name = "high_density_stability_v101"
        Pes = collect(4.:0.5:11.)
        T  = 1.0
        max_runs = 5
for ρa in collect(0.2:0.05:0.45), Dθ in [100., 10000.]
        local param
        λ_step = (0.5)*sqrt(Dθ)
        λs = Pes*sqrt(Dθ)
        λmax = (11.)*sqrt(Dθ)
        #
        param = pde_param(; name = name, ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = λ_step, λmax = λmax, λs=λs)
        #
        push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
length(params)/nworkers()
pmap(run_stab_search_stupid_mode, params; distributed = true, batch_size=1, on_error=nothing,)
#