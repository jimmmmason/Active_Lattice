cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

#varying parameters
params = []
name = "high_density_active"
for ρ in [0.96,0.97,0.98,0.99]
for λ in [1,8,16,32]
        local param
        param = uniform_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=64, Δt = 0.01)
        push!(params,param)
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
    println("success")
end