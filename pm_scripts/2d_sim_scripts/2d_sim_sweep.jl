#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed
addprocs([("adrastea",1)]) 
addprocs([("refract",1)])
addprocs([("radius",1)])
addprocs([("heart",1)])
addprocs([("kenku",4)])
addprocs([("leona",1)])
addprocs([("saz",1)])
nworkers()
Distributed.interrupt()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/2d_sims.jl");
#
# Create parameters
params = []
DT, DR, N, Δx, Lx, Ly, δt, δ, ϕa, ϕp = (1.0, 1.0, 32, 0.01, 2.0, 2.0, 1e-5, 0.01, 0.5, 0.0);
T, save_interval, param_name = (20.0, 0.01, "test")
map([10., 20., 30. , 40.]) do v0
    param = _2d_new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
    push!(params,param)
end
#
pmap(load_and_run_sim, params)
#






