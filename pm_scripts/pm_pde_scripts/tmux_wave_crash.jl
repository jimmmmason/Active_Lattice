cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
# Load relevant fuctions
# include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");

# Load initial wave
    param = get_grid_param(12,4) # get_dense_param(40, 0.025)
    T = 50000.0
    @pack! param = T
    @unpack Nx = param
    loaded, f, t = quiet_load_last_pde(param)
#

# create new param
    param = get_grid_param(12,4)
    param["Nx"] = 2*param["Nx"]
    param["name"] = "wave_crash"
    param["save_interval"] = 1.0
    param["T"] = 200.0;
# rotate f 
    f = circshift(f,(275,0));
# create new f
    f = vcat(f,f[end:-1:1,[2,1,3]]);
#

# loaded, g, t = quiet_load_last_pde(param)

run_current_pde_sym(param, 200.0, f,0.0)

# export PATH="/home/jm2386/.local/ /bin:${PATH}"
# nice -n 19 julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/tmux_wave_crash.jl