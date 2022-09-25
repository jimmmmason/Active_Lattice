cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

#varying parameters
params = []
name = "high_density_active"
for ρ in [0.5,0.8,0.9,0.95]
for λ in [1,8,16,32]
        local param
        param = uniform_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=64, Δt = 0.01)
        push!(params,param)
end
end
#run sims
@distributed for param ∈ params
    T = 1.0
    #param = uniform_initial_param(L=32, λ = 20, ρa = 0.5, ρp = 0.0, Δt = 0.01);
    model = initialize_model(param);
    run_model_until!(param, model, T; save_on = true);
    println("success")
end