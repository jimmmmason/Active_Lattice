#Set up DrWatson
cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
# Load relevant fuctions
include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
# set parameters
params = []
Pe = 7.5
name = "pm_pde_binodal+pert_run_δ=$(δ)_l=$(1/sqrt(Dθ))"
pert = "pm_lin"
    T = 20.0
    save_interval = 0.001
    δ  = 1e-3
    Dx = 1. 
    Dθ = 400.0
    Nx = 2^10
    Nθ = 2
    δt = 1e-7
#create params
params = []
γs = [ 2.5, #unstable complex no bin
        2.25, #unstable real no bin 
        2.0, #unstable real unstable bin
        1.5, #unstable real stable bin
    ]
ϕas = fill(0.7, length(γs))
ϕps = (γs .-1 ).*(-ϕas .+1)./ γs
map(ϕas, ϕps) do ϕa, ϕp
        ρ = ϕa + ϕp
        χ = ϕa/ρ
        param = pde_param_pm(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                                save_interval = save_interval, max_steps = 1e9,
                                pert = pert, δ = δ,
        )
        push!(params,param)
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
@everywhere include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
#
pmap(pert_pde_run_pm, params; distributed = true, batch_size=1, on_error=nothing,)
###  
pmap(load_pde_run_pm, params; distributed = true, batch_size=1, on_error=nothing,)
###