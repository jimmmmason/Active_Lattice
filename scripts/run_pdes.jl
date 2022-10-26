cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")

#varying parameters
params = []
name = "high_density_stability_v4"
for ρa in collect(0.90:0.02:0.98)
for λ in [60.]
        local param
        #
        param = pde_param(; name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7)
        #
        push!(params,param)
end
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
pmap(perturb_pde_run, main_pool, params; distributed = true, batch_size=1, on_error=nothing,)
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
name = "high_density_stability_v4"
ρa = 0.92
λ = 40.
t = 0.10
param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
        data = wload(filename)
        fa = data["fa"]
        fp = data["fp"]
        t = data["t"]
        @pack! density = fa, fp, t
        fig, ax = PyPlot.subplots(figsize =(10, 10))
        plot_pde_mass(fig,ax,param,density)
        #plot_pde_mag(fig,ax,param,density)
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
fig = figure(figsize=(30,10))
i = 1 
name = "high_density_stability_v4"
for λ in [50., 40., 30., 20., 10.]
for ρa in collect(0.05:0.05:0.85)
        ax = fig[:add_subplot](5,17,i)
        param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.001)
        plot_error(fig,ax, param, t_saves, fa_saves, fp_saves; t_max = 1.0, y_max = 0.2)
        if λ == 10.
                ax.set_xlabel("t")
        end
        if ρa == 0.05
                ax.set_ylabel("e")
        end
        i += 1
end
end
fig.tight_layout()
display(fig)
PyPlot.savefig("/store/DAMTP/jm2386/Active_Lattice/data/time_stability_test.pdf",dpi = 300, format = "pdf")

###
stabdata = refresh_stab_data(; ρs = 0.05:0.05:0.95,   Dθ = 10., Nx = 50, Nθ = 20, λs = 5.:5.:100., name = "high_density_stability_v4")
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.05:0.05:0.95 ,xs = collect(0.01:0.01:0.99), xtic = 0:0.2:1, ytic = 0:10:100, axlim = [0., 1., 15., 100.])
display(fig)


#
ρ = 0.6
param = pde_param(; name = name, λ = λ , ρa = ρ, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7)
#
unfinished, next, param = next_param(stabdata, param; λmax = 1e8, λ_step = 5.)

