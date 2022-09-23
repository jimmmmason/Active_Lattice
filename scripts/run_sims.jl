cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("./src/sim_functions.jl")
@everywhere include("./src/sim_functions.jl")

#varying parameters
params = []
name = "high_density_active"
for ρ in [0.5,0.8,0.9,0.95]
for λ in [1,8,16,32]
        local param
        param = uniform_initial_param(; name = name, λ = λ ,ρa = 0.0, ρp = ρ, L=64)
        push!(params,param)
end
end
#run sims
@distributed for param in params
    try
        model = initialize(param)
        run_model_intervals!(param,model,1.; interval = 0.001, save_on =true);
        println("success")
    catch
        println("fail")
    end
end