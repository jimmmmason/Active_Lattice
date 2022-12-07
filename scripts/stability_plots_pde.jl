cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
#=

Stability parameters

=#
###
name = "stability_1d_actpass_2"
T  = 1.0
Dθ = 1.
ρp = 0.3
ρs  = collect(0.4:0.05:0.95)
Pes = collect(5.:5.:100.)
#Pes = collect(0.:0.5:10.)
    Nx = 50
    Nθ = 20
    λs = Pes*sqrt(Dθ)
###
###
name = "stability_1d_actpass_2"
T  = 1.0
Dθ = 100.
ρp = 0.0
ρs  = collect(0.4:0.05:0.95)
Pes = collect(0.:5.0:100.)
    Nx = 50
    Nθ = 20
    λs = Pes*sqrt(Dθ)
###
#=

Get stability data from saves

=#
###
stabdata = refresh_stab_data_1d(; ρs = ρs, Dθ = Dθ, Nx = Nx, Nθ = Nθ, λs = λs, name = name, ρp = ρp)
###
#=

Load previously saved stability

=#
###
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ).jld2"
stabdata = wload(filename)
###
###
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_ρp=$(ρp).jld2"
stabdata = wload(filename)
###
#=

Plot stability

=#
###
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = ρs ,xs = collect(0.01:0.01:0.99), xtic = (minimum(ρs)-ρp):0.2:(1-ρp), ytic = 0:20:maximum(Pes), axlim = [minimum(ρs)-ρp, 1-ρp, minimum(Pes), maximum(Pes)], Dθ = Dθ,ρp=ρp)
display(fig)
###


