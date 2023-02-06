cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

@everywhere include("/home/jm2386/Active_Lattice/src/sim_functions.jl")

#varying parameters
params = []
name = "presentation2"
Dθ =100.
Pes = [20.]
λs = Pes*sqrt(Dθ)/2
for ρa in [0.4], ρp in [0.2]
for λ in λs
        local param
        param = extra_mixing_initial_param(; name = name, λ = λ ,ρa = ρa, ρp = ρp, L = 128, Δt = 0.001, γ = 0., T = 0.5, Dθ =Dθ)
        push!(params,param)
end
end
param = params[1]
#run sims
pmap(run_sim, params; distributed = true, batch_size=1, on_error=nothing,)
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
PyPlot.savefig("translation_invariance_ρ=$(ρ).pdf",dpi = 300, format = "pdf")
#create vids
T = 0.499
@distributed for param ∈ params
    T = 0.5
    t_saves, η_saves = load_etas(param, T; dump_interval = 0.1)
    animate_etas(param,t_saves,η_saves)
end
t = 0.4
η = η_saves[length(t_saves)]
#load plots 
PyPlot.close("all")
fig, ax = figure(figsize=(10,10))  
i =1
name = "comparison"
t = 0.499
L = 128
fig, ax = PyPlot.subplots(figsize =(10, 10))
r = 7
clf()
for ρ in [0.7]
for λ in λs
    ax = fig[:add_subplot](2,2,i)
    param = extra_mixing_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=L, Δt = 0.01, γ = 0., T = 1.0,Dθ =100.)
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
    η, t = wload(filename, "η", "t")
    plot_eta(fig,ax,param, t, η; title = false, r = r)
    #h = density_hist(fig,ax,param, t, η; title = false)
    i += 1
    #=
    if λ == 8
        ax.yaxis.set_ticks(0:0.5:1)
        ax.set_ylabel("ρ = $(ρ)")
    end
    =#
    if ρ == 0.7
        ax.xaxis.set_ticks(0:0.5:1)
        ax.set_xlabel("Pe = $(λ*sqrt(Dθ)/2), ρ = $(ρ), Dθ = $(Dθ)")
    end
end
end
display(fig)
PyPlot.savefig("active_passive_sim_example_2.pdf",dpi = 300, format = "pdf")
# hist plot 
PyPlot.close()
fig = figure(figsize=(10,10))
i = 1 
name = "cross_density_hist_data_4"
t_start = 1.5
t_end = 1.7
L = 128
r = 5
for ρ in [0.7]
for λ in λs
    ax = fig[:add_subplot](2,2,i)
    param = extra_mixing_initial_param(; name = name, λ = λ ,ρa = ρ, ρp = 0., L=L, Δt = 0.01, γ = 0., T = 3.0,Dθ =100.)
    t_saves , η_saves = load_etas(param, t_end; dump_interval = 0.1, start_time = t_start);
    time_density_hist(fig, ax, param, t_saves, η_saves; r = r, bins = (2*r+1)^2 )
    #=
    if λ == 12
        ax.set_ylabel("ρ = $(ρ)")
    end
    =#
    if ρ == 0.8
        ax.set_xlabel("v₀ = $(λ), ρ = $(ρ)")
    end
    i += 1
end
end
display(fig)
PyPlot.savefig("hist_phase_sep_parameter_range_L=$(L).png",dpi = 300, format = "png")

h = randn(100)
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.hist(h; bins = 20)
display(fig)


display(fig)
PyPlot.savefig("time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).png",dpi = 300, format = "png")



λ = 14
ρ = 0.8