#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
# Load relevant fuctions
include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
include("/home/jm2386/Active_Lattice/src/plot_functions.jl");
#load data 
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)_compressed.jld2"
saves = load(filename)["saves"]

#compress saves if loaded
for (param, t_saves, fa_saves, fp_saves) in saves[2:1:4]
    @unpack T, save_interval, Nx, Nθ= param
    start_time = maximum(t_saves)
    T = 20.0
    save_interval = 0.0001 # one deeper to pick up offest of 1e-4 
    t_s, fa_s, fp_s = load_pdes_pm(param,T; save_interval = save_interval, start_time = start_time)
    append!(fa_saves,fa_s )
    append!(fp_saves,fp_s )
    append!(t_saves,t_s )
    push!(saves, (param, t_saves, fa_saves, fp_saves))

    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)_compressed.jld2"
    save(filename, "saves", saves)

    println("completed anther")
end

filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)_compressed_T=20.jld2"
save(filename, "saves", saves_shortened)
