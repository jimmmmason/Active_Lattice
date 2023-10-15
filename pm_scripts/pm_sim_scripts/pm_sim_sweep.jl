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
#
params = []
γs = [ 2.5, #unstable complex no bin
        2.25, #unstable real no bin 
        2.0, #unstable real unstable bin
        1.5, #unstable real stable bin
]
ϕas = fill(0.7, length(γs))
ϕps = (γs .-1 ).*(-ϕas .+1)./ γs
DT, v0, DR, N, Lx, Ly = (1.0, 7.5, 1.0, 100, 8, 1);
T = 1.
sim_name = "sim_run_1"
save_interval = 0.001
map(ϕas, ϕps) do ϕa, ϕp
    sim_param = new_sim_param(DT, v0, DR, N, Lx, Ly, ϕa, ϕp; T = T, name = sim_name, save_interval = save_interval, save_on = true);   
    push!(params,sim_param)
end
#
pmap(load_and_run_sim, params)
#






