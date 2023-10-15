#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
# Load relevant fuctions
include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
# set parameters
Pe = 7.5
pert = "pm_lin"
    T  = 4.0
    δ  = 1e-2
    save_interval = 0.01
    Dx = 1. 
    Dθ = 160000.0
    Nx = 512
    Nθ = 2
name = "pm_pde_run_δ=$(δ)_l=$(1/sqrt(Dθ))"
#
params = []
Δx = 0.1
ϕs = Δx:Δx:(1-Δx)
for ϕa in ϕs, ϕp in ϕs
    ρ = ϕa+ϕp
    if ρ < 1
        χ = ϕa/ρ
        param = pde_param_pm(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                                save_interval = save_interval, max_steps = 1e7,
                                pert = pert, δ = δ,
        )
        push!(params,param)
    end
end
# Add remote computers
using Distributed
addprocs([("adrastea",3)])
addprocs([("refract",3)])
addprocs([("radius",3)])
addprocs([("heart",3)])
addprocs([("kenku",3)])
addprocs([("leona",3)])
addprocs([("saz",3)])
nworkers()
# Load relevant fuctions
@everywhere include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
#
pmap(perturb_pde_run_pm, params; distributed = true, batch_size=1, on_error=nothing,)
###