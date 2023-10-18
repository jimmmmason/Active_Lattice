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
@everywhere include("/home/jm2386/Active_Lattice/src/pm_sims.jl");
# Create parameters
params = []
ϕas = [0.5, 0.3]
ϕps = [0.1, 0.6]
DT, v0, DR, N, Lx, Ly = (1.0, 20.0, 1.0, 100, 2.0, 0.5);
T = 20.
sim_name = "sim_run_2"
save_interval = 0.01
map(ϕas, ϕps) do ϕa, ϕp
    sim_param = new_sim_param(DT, v0, DR, N, Lx, Ly, ϕa, ϕp; T = T, name = sim_name, save_interval = save_interval, save_on = true);   
    push!(params,sim_param)
end
# Run Sims 
pmap(load_and_run_sim, params)
# (Running on session 14)






