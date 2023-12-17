#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed, ClusterManagers

# addprocs_slurm(1,job_name="slurm_test", time="00:00:30", exeflags = "/home/jm2386/Active_Lattice/", exename="/home/jm2386/julia-1.8.1/bin/julia")



for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println(id, " " , pid, " ", host)
end

for i in workers()
    rmprocs(i)
end

@sync @distributed for i in 1:nworkers()
    println("hello from $(myid()):$(gethostname())")
end







addprocs([("refract",2)])
addprocs([("radius",2)])
addprocs([("terbium",2)])



Distributed.interrupt()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
@everywhere include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
# Repeating Soleton
New_Params = []
DT, v0, DR, N, Δx, Lx, Ly, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 20.0, 0.5, 1e-5, 0.1);
# T, save_interval, param_name, pert = (0.1, 0.001, "periodic_solon_test", "lin")
T, save_interval, param_name, pert = (5000.0, 1.0, "periodic_solon_test", "lin")
Lxs = [2*Lx, 3*Lx, 4*Lx]
ϕas = [0.35, 0.37,0.38]
ϕps = fill(0.3,3)
map(ϕas, ϕps) do ϕa, ϕp
    for Lx in Lxs
        param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
        push!(New_Params,param)
    end
end
###
pmap(force_load_and_pert_pde, New_Params)
###






