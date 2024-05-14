cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_sims.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");

# Load initial wave
param = get_grid_param(21,11)
loaded, f, t = load_last_pde(param)

param["save_interval"] = 10.0
param["name"] = "soliton_stretch"
@unpack Nx = param
for Lx in 19:(-1):1
    global f,t,param,Nx
    local T
    #change Lx
    print("running Lx = $(Lx)")
    param["Lx"] = Float64(Lx)
    param["Δx"] = Float64(Lx/Nx)

    # allow equilbriation 
    T = 40.0
    t = 0.0
    t, f =  run_current_pde(param,T, f,t)
end
# export PATH="/home/jm2386/julia-1.8.1/bin:${PATH}
# nice -n 19 julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/pm_stretcher.jl