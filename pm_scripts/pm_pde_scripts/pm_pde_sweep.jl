#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed
addprocs([("refract",2)])
addprocs([("radius",2)])
addprocs([("heart",1)])
addprocs([("kenku",4)])
addprocs([("leona",1)])
addprocs([("saz",2)])
nworkers()
@sync @distributed for i in 1:nworkers()
    println(gethostname())
end
Distributed.interrupt()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
@everywhere include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
# Figure 1
# Create parameters
params = []
ϕas = [0.5, 0.3]
ϕps = [0.1, 0.6]
DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 20.0, 1.0, 100, 0.01, 2.0, 0.5, 1e-5, 0.01);
T, save_interval, param_name = (30.0, 0.0001, "fig_4")
map(ϕas, ϕps) do ϕa, ϕp
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
    push!(params,param)
end
# Figure 2
#create parameters
params = []
ϕas = fill(0.5,8)
ϕps = collect(0.025:0.025:0.2)
DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 7.5, 1.0, 100, 0.02, 20.0, 0.5, 1e-4, 2.0);
T, save_interval, param_name = (10000.0, 10.0, "fig_2")
map(ϕas, ϕps) do ϕa, ϕp
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
    push!(params,param)
end
# Run Sims
#pmap(load_and_run_pde, params)
pmap(run_new_pde, params)
# (Running on radius[my office] session 16)






