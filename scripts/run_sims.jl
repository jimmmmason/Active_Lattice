cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

#varying parameters
params = []
name = "high_density_mixing"
for ρ in [0.95, 0.99]
for γ in [0., 0.1]
for λ in [8, 10, 12, 14]
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
    for i in 1:5
        try 
            t = 0.99+1.0*(i-1)
            @unpack name, L, λ, γ, ρa, ρp, Δt, = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).jld2";
            η, t = wload(filename, "η", "t")
            model["η"] = η
            model["t"] = t
        catch
            println("no preload $(i)")
        end
    end
    println("starting run: γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])")
    run_model_until!(param, model, T; save_on = true);
    println("success")
end
#plot symmetry
clf()
ρ = 0.95
fig, ax = PyPlot.subplots(figsize =(10, 10))
    params = []
    name = "high_density_mixing"
    for γ in [0., 0.1]
    for λ in [8, 10, 12, 14]
                local param
                param = extra_mixing_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=64, Δt = 0.01, γ = γ)
                push!(params,param)
    end
    end
for param ∈ params
    local y
    @unpack L, Ω = param
    T = 5.0
    t_saves, η_saves = load_etas(param, T)
    if t_saves == []
        println(param["γ"])
        println(param["ρa"])
        println(param["λ"])
    else
        y = translation_invariance.(η_saves;Ω = Ω,L= L)
        ax.plot(y; label = "γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])" )
    end
end
legend()
display(fig)
PyPlot.savefig("translation_invariance_ρ=$(ρ).png",dpi = 300, format = "png")
#create vids
@distributed for param ∈ params
    T = 5.0
    t_saves, η_saves = load_etas(param, T)
    animate_etas(param,t_saves,η_saves)
end
#load final
name = "high_density_mixing"
param = extra_mixing_initial_param(; name = name, λ = 12 ,ρa = 0.95, ρp = 0., L=64, Δt = 0.01, γ = 0.)
t = 4.0
@unpack name, L, λ, γ, ρa, ρp, Δt, = param
filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).jld2";
η, t = wload(filename, "η", "t")
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_eta(fig,ax,param, t, η)
display(fig)
PyPlot.savefig("time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).png",dpi = 300, format = "png")



