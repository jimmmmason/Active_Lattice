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
@everywhere include("/home/jm2386/Active_Lattice/src/2d_pdes.jl");
#
# Create parameters
#
# MIPS
#
Params = []
DT, DR, N, Nθ, Δx, Lx, Ly, δt, δ, ϕa, ϕp = (1.0, 1.0, 32, 20, 0.04, 2.0, 2.0, 1e-5, 0.1, 0.45, 0.0);
T, save_interval, param_name = (20.0, 0.01, "test")
map([50., 60. , 70., 80.]) do v0
    param = _2d_new_param(DT, v0, DR, N, Nθ, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
    push!(Params,param)
end
#
pmap(load_and_run_pde, Params)
#
# NRIPS
#
DT, DR, N, Nθ, Δx, Lx, Ly, δt, δ, ϕa, v0 = (1.0, 1.0, 32, 20, 0.2, 10.0, 10.0, 1e-5, 0.1, 0.5, 10.0);
T, save_interval, param_name = (20.0, 0.1, "test")
Params = []
map(0.34:0.02:0.4) do ϕp
    param = _2d_new_param(DT, v0, DR, N, Nθ, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true)
    push!(Params,param)
end
#
pmap(load_and_run_pde, Params)
#





