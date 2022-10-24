cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")

#varying parameters
params = []
name = "high_density_stability_v4"
for ρa in [0.9, 0.92, 0.94, 0.96, 0.98]
for λ in [10., 20., 30., 40., 50.]
        local param
        #
        param = pde_param(; name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        #
        push!(params,param)
end
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#plot pde 
perturb_pde_run(param)
# /store/DAMTP/jm2386/Active_Lattice/data/pde_raw/
t_saves, fa_saves, fp_saves = load_pdes(param,0.2; save_interval = 0.001)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = 4
t, fa, fp = t_saves[n], fa_saves[n], fp_saves[n]

display(fig)
#plot_pde_mag(fig,ax,param,density)
#display(fig)
name = "high_density_stability_v3"
ρa = 0.9
λ = 0.
t = 0.
param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
        data = wload(filename)
        fa = data["fa"]
        fp = data["fp"]
        t = data["t"]
        @pack! density = fa, fp, t
        fig, ax = PyPlot.subplots(figsize =(10, 10))
        #plot_pde_mass(fig,ax,param,density)
        plot_pde_mag(fig,ax,param,density)
        display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
param = pde_param(; name = name, λ = λ , ρa = ρ, ρp = 0., T = 0.4, Dθ = 10., δt = 1e-5)
t_saves, fa_saves, fp_saves = load_pdes(param,0.3; save_interval = 0.003)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(t_saves,dist_saves)
display(fig)
#
PyPlot.close()
fig = figure(figsize=(10,10))
i = 1 
name = "high_density_stability_v3"
for ρa in [0.9, 0.92, 0.94, 0.96, .98]
for λ in [10., 20., 30., 40., 50.]
        ax = fig[:add_subplot](5,5,i)
        param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.001)
        plot_error(fig,ax, param, t_saves, fa_saves, fp_saves; t_max = 1.0, y_max = 0.2)
        if λ == 10.
                ax.set_ylabel("e")
        end
        if ρa == 0.98
                ax.set_xlabel("t")
        end
        i += 1
end
end
fig.tight_layout()
display(fig)
PyPlot.savefig("/store/DAMTP/jm2386/Active_Lattice/data/time_stability_test.pdf",dpi = 300, format = "pdf")


