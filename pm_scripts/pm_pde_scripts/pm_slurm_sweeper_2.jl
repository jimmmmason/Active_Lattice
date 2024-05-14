#Set up
cd("/home/jm2386/Active_Lattice/");
# import Pkg; Pkg.add("DrWatson")
using DrWatson;
@quickactivate "Active_Lattice";
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");

# Access the command-line argument 'i'
i = parse(Int64, ARGS[1]);

# run solution i 
param = get_active_param(i)
loaded, f, t = load_last_pde(param)
t = 0.0
param["save_interval"] = 10.0
param["name"] = "soliton_vertical_sweep"

if loaded
    dt = 40.0
    print("running solution $(param["ϕa"]) $(param["ϕp"]) ",gethostname())
    t , f = run_current_pde(param,dt, f,t)
else
    print("invalid load $(param["ϕa"]) $(param["ϕp"]) ",gethostname())
end

# load_and_run_pde(param)
# comp_wave_speed(param)
