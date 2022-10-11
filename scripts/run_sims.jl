cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

#varying parameters
params = []
name = "high_density_mixing"
for ρ in [0.7, 0.8, 0.9, 0.95, 0.99]
for γ in [0., 0.01, 0.1]
for λ in [4, 8, 12, 16]
        local param
        param = extra_mixing_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=64, Δt = 0.01, γ = γ)
        push!(params,param)
end
end
end
#run sims
@distributed for param ∈ params
    T = 5.0
    model = initialize_model(param);
    run_model_until!(param, model, T; save_on = true);
    println("success")
end
#plot symmetry
fig, ax = PyPlot.subplots(figsize =(10, 10))
for param ∈ params
    local y
    @unpack L, Ω = param
    T = 1.0
    t_saves, η_saves = load_etas(param, T)
    y = translation_invariance.(η_saves;Ω = Ω,L= L)
    ax.plot(y)
end
display(fig)
#create vids
@distributed for param ∈ params
    T = 5.0
    t_saves, η_saves = load_etas(param, T)
    animate_etas(param,t_saves,η_saves)
end