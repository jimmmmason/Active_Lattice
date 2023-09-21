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
    T = 0.05
    save_interval = 0.001
    δ  = 1e-4
    Dx = 1. 
    Dθ = 400.0
    Nx = 2^10
    Nθ = 2
    δt = 1e-8
#
#load binodal values
name = "pm_pde_binodal_run_δ=$(δ)_l=$(1/sqrt(Dθ))"
#
filename = "/store/DAMTP/jm2386/Active_Lattice/data/binodal/Pe=$(Pe).jld2"
data = wload(filename)
@unpack Pe, γs, ϕ1s, ϕ2s, average_ϕs, χs = data
#create params
params = []
γ_length = length(γs)
n_lines = 40
interval = Int(round(γ_length/n_lines))
γs          = γs[interval:(interval):Int64(round(γ_length))]
average_ϕs  = average_ϕs[interval:(interval):Int64(round(γ_length))]
map(average_ϕs, γs) do ρ, γ
        χ = (1-γ*(1-ρ))/ρ
        param = pde_param_pm(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                                save_interval = save_interval, max_steps = 1e7,
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
pmap(binodal_pde_run_pm, params; distributed = true, batch_size=1, on_error=nothing,)
###