#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed
addprocs([("refract",2)])
addprocs([("radius",2)])
addprocs([("terbium",2)])
addprocs([("kenku",3)])
addprocs([("leona",1)])
addprocs([("saz",2)])
@sync @distributed for i in 1:nworkers()
    println(gethostname())
end
nworkers()
Distributed.interrupt()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
@everywhere include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
# Figure 1


param = get_grid_param_wide(1,1)
param["T"] = 2000.0
param["ϕa"] = 0.5
param["ϕp"] = 0.3
param["save_interval"] = 1.0
run_new_pde(param)







# Create parameters
# params = []
# ϕas = [0.5, 0.3]
# ϕps = [0.1, 0.6]
# DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 20.0, 1.0, 100, 0.01, 2.0, 0.5, 1e-5, 0.01);
# T, save_interval, param_name = (30.0, 0.0001, "fig_4")
# map(ϕas, ϕps) do ϕa, ϕp
#     param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
#     push!(params,param)
# end
# # Figure 2
# #create parameters
# params = []
# # ϕas = fill(0.6,5)
# # ϕps = collect(0.01:0.01:0.05)
# ϕas = [0.45,0.5,0.55]
# ϕps = fill(0.25,3)
# DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 20.0, 0.5, 1e-4, 0.1);
# T, save_interval, param_name = (4000.0, 1.0, "periodic_plot")
# map(ϕas, ϕps) do ϕa, ϕp
#     param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
#     push!(params,param)
# end
# # Run Sims
# #pmap(load_and_run_pde, params)
# pmap(run_new_pde, params)
# # (Running on radius[my office] session 16)
# #periodic load 
# #create parameters
# #create parameters
# Params = []
# ϕas = collect(0.34:0.02:0.38)
# ϕps = fill(0.3,3)
# DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.1);
# #T, save_interval, param_name = (19.8, 0.01, "fig_1")
# T, save_interval, param_name, pert = (2000.0, 1.0, "periodic_pert_plot", "lin")
# map(ϕas, ϕps) do ϕa, ϕp
#     param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
#     push!(Params,param)
# end
# ###
# pmap(run_and_pert_pde, Params)
# ###

# param = Params[3]
# t, f = load_pde(param,2000.)
# param["save_interval"] = 0.001
# run_current_pde(param,1000., f,t)





