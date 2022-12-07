cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
#=

individual parameters
v101 runs with deterministic pert
v102 runs with additional random pert
v103 runs from an extremely seperated state
=#
###
name = "high_density_stability_v101"
ρa = 0.75
Dθ = 100.
Pe = 4.7
    ρp = 0.
    δt = 1e-5
    Nx = 50
    Nθ = 20
    λ  = Pe*sqrt(Dθ)
###
#=

Plot visual at T 

=#
###
t = 0.018
    @unpack ρa,ρp,Dθ,Pe,δt,λ,Nx,Nθ =param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2"
    data = wload(filename)
    fa = data["fa"]
    fp = data["fp"]
    @pack! density = fa, fp, t
    @pack! param = ρa,ρp,Dθ,Pe,δt,λ,Nx,Nθ
fig, axs = plt.subplots(1, 2, figsize=(12,5))
plot_pde(fig,axs,param,density)
display(fig)
#t += 0.1
#for i in 1:100 pde_step!(param,density) end
###
#=

Plot time series of stability

=#
###
param = pde_param(;name = name, Pe = Pe , ρa = ρa, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
    t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.01)
    dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    ax.plot(t_saves,dist_saves)
    ax.set_xlabel("t")
    ax.set_ylabel("‖ρ-ϕ‖₂")
    ax.set_title("Dθ = $(Dθ) ρ = $(ρa) Pe = $(Pe)")
display(fig)
###
###
#=

Create video of pde

=#
###
param = pde_param(;name = name, Pe = Pe , ρa = ρa, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
param = params[2]
    t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.0001)
    animate_pdes(param,t_saves,fa_saves,fp_saves)
###
