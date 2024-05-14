#Set up
cd("/home/jm2386/Active_Lattice/");
# import Pkg; Pkg.add("DrWatson")
using DrWatson;
@quickactivate "Active_Lattice";
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");

# Access the command-line argument 'i'
i = parse(Int64, ARGS[1]);

# run solution
param = get_soliton_param(i)
# change name in get_active_param(i)
param["name"] = "soliton_vertical_sweep"

loaded, f, t = load_last_pde(param)
if loaded
    global t, f
    param["Nx"] = 400
    param["save_interval"] = 5.0
    Δphi = 0.01
    dt = 40.0
    t = 0.0
    while (param["ϕp"]+param["ϕa"]) < 0.95
        global t, f
        param["ϕp"] = round(param["ϕp"] + Δphi; digits = 2)
        f[:,3] = f[:,3] .+ Δphi
        if is_valid(f, param)
            global t, f
            print("running solution $(param["ϕa"]) $(param["ϕp"]) ",gethostname())
            t, f = run_current_pde(param,dt, f,t)
        else
            global t, f
            print("f invalid  ",gethostname())
        end
    end
else
    print("load failed  ",gethostname())
end

# srun -J "incr phip" -c 1 -t 24:00:00 --nice=19 -o out.txt -e err.txt  julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/alter_param_scripts/pm_increase_phip.jl &
