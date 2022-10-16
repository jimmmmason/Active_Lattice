cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")

#varying parameters
params = []
name = "high_density_stability"
for ρ in [0.95, 0.96, 0.97, 0.98, 0.99]
for λ in [20, 30, 40, 50, 60]
        local param
        param = pde_param(; name = name, λ = λ , ρa = ρ, ρp = 0., T = 0.3 )
        push!(params,param)
end
end
#run pdes
pmap(pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#plot pde 
# /store/DAMTP/jm2386/Active_Lattice/data/pde_raw/
t_saves, fa_saves, fp_saves = load_pdes(param,0.2; save_interval = 0.001)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = 4
t, fa, fp = t_saves[n], fa_saves[n], fp_saves[n]
@pack! density = fa, fp, t
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_pde_mass(fig,ax,param,density)
display(fig)
#plot_pde_mag(fig,ax,param,density)
#display(fig)
