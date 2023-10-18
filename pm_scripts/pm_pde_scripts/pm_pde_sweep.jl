#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed
addprocs([("adrastea",1)])
addprocs([("refract",1)])
addprocs([("radius",1)])
addprocs([("heart",1)])
addprocs([("kenku",1)])
addprocs([("leona",1)])
addprocs([("saz",1)])
nworkers()
Distributed.interrupt()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
# Create parameters
params = []
ϕas = [0.5, 0.3]
ϕps = [0.1, 0.6]
DT, v0, DR, N, Lx, ϕa, ϕp, δt, δ = (1.0, 20.0, 1.0, 200, 2.0, 0.3, 0.6, 1e-5, 0.01);
T = 20.0
pde_name = "pde_run_2"
save_interval = 0.01
map(ϕas, ϕps) do ϕa, ϕp
    pde_param = new_pde_param(DT, v0, DR, N, Lx, ϕa, ϕp, δt, δ; T = T, name = pde_name, save_interval = save_interval, save_on = true);  
    push!(params,pde_param)
end
# Run Sims 
pmap(load_and_run_pde, params)
# (Running on session 14)






