cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("booted")
##
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random

function uniform_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, d=2, Δt = 0.001)
    param = Dict{String,Any}()  
    #this is the only dimension dependent part:
    Ω = [[i,j] for i in 1:L for j in 1:L] 
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))
    γ = 0
    function angles(x,n) 
        if n == 1
            return 2*π*rand()
        else
            return -1
        end
    end
    function rates(n,m,i)
        E = [[1,0],[0,1],[0,-1],[-1,0],]
        if m[1]>0
            return 0
        elseif n[1]==0
            return 0.
        elseif n[2]==2
            return L^2*D
        else
            return L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
        end
    end
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, Ω, E, site_distribution, angles, rates
    return param
end

function extra_mixing_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, d=2, Δt = 0.001, γ = 0.01)
    param = Dict{String,Any}()  
    #this is the only dimension dependent part:
    Ω = [[i,j] for i in 1:L for j in 1:L] 
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))
    function angles(x,n) 
        if n == 1
            return 2*π*rand()
        else
            return -1
        end
    end
    function rates(n,m,i)
        E = [[1,0],[0,1],[0,-1],[-1,0],]
        if m[1]>0
            return L^2*γ
        elseif n[1]==0
            return 0.
        elseif n[2]==2
            return L^2*D + L^2*γ
        else
            return L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] + L^2*γ
        end
    end
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, Ω, E, site_distribution, angles, rates
    return param
end

function initialize_model(param::Dict{String,Any})
    @unpack name, L, D, λ, γ, ρa, ρp, Δt, Ω, E, site_distribution, angles, rates = param
    # create configuration, rates and jumps
    η = fill([],(L,L))
    c = fill(0.,(L,L,4))
    j = fill([],(L,L,4))
    #fill configuration
    for x ∈ Ω
        local w, n
        w = Weights(site_distribution[x...])
        #fill model
        n = sample([0,1,2],w)
        η[x...] = [n, angles(x,n)]
    end
    #fill rates and jumps
    for x ∈ Ω 
        for i in 1:4
            local y
            #find adjacent site
            y  = (x + E[i] +[L-1,L-1]) .% L + [1,1]
            #fill rates 
            c[x...,i] = rates(η[x...],η[y...],i)
            #fill jump vectors
            j[x...,i] = [x,y]
        end
    end
    #pack into model
    α = sum(c)
    Δτ = 0.01/α
    t = 0.
    model = Dict{String,Any}()
    @pack! model = η, c, j, t, α, Δτ
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, Ω, E, site_distribution, angles, rates
    return model
end

function model_step!(param::Dict{String,Any},model::Dict{String,Any})
    @unpack L, Ω, E, rates = param
    @unpack η, c, j, t, α, Δτ, = model
    #increase time 
    Δt = round(α)*Δτ #roughly Δt
    t += Δt  
    #see if jump occurs
    for i in 1:round(α)
        if rand() < α*Δτ 
        #select jump
            w     = Weights( [(c...)...])
            jump  = sample(j, w)
        #execute jumps
            η[jump[2]...], η[jump[1]...] = η[jump[1]...], η[jump[2]...]
        #correct propensity
            for x in jump
                for i in 1:4
                    local y
            #find adjacent site
                    y  = (x + E[i] +[L-1,L-1]) .% L + [1,1]
            #correct new rates 
                    c[x...,i]   = rates(η[x...],η[y...],i  )
                    c[y...,5-i] = rates(η[y...],η[x...],5-i)
                end
            end
            α = sum(c)
        end
    end
    #diffuse angles
    for x ∈ Ω
        if η[x...][1] == 1
            η[x...][2] = (η[x...][2] + sqrt(Δt)*randn() + 2*π) % (2*π)
        end
    end
    @pack! model = η, c, t
end

function run_model_until!(param::Dict{String,Any},model::Dict{String,Any},T; save_on =false)
    while model["t"] < T 
        model_step!(param,model)
        if save_on
            @unpack name, L, λ, γ, ρa, ρp, Δt = param
            @unpack η, t = model
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).jld2";
            data = Dict{String,Any}();
            @pack! data = param, η, t
            safesave(filename,data)
        end
    end
end

function run_and_dump_sim(param::Dict{String,Any},model::Dict{String,Any},T; dump_interval = 0.01, save_on =false)
    while model["t"] < T
        local t
        t = deepcopy(model["t"]+dump_interval)
        run_model_until!(param,model,t; save_on =false)
        if save_on
            @unpack name, L, λ, γ, ρa, ρp, Δt = param
            @unpack η, t = model
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).jld2";
            data = Dict{String,Any}();
            @pack! data = param, η, t
            safesave(filename,data)
        end
    end
end


function load_etas(param::Dict{String,Any},T; dump_interval = 0.01)
    s =0.
    t_saves = []
    η_saves = []
    while s < T
        try 
            @unpack name, L, λ, γ, ρa, ρp, Δt = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).jld2";
            data = wload(filename)
            @unpack  η, t =data
            push!(t_saves,t)
            push!(η_saves,η)
        catch
        end
        s += dump_interval
    end
    return t_saves, η_saves
end
##
#visualise data
using PyPlot, PyCall
@pyimport matplotlib.animation as anim

function animate_etas(param,t_saves,η_saves)
    frames = length(η_saves)-1
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    makeframe(i) = plot_eta(fig,ax,param, t_saves[i+1], η_saves[i+1])
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=20)
    # Convert it to an MP4 movie file and saved on disk in this format.
    @unpack name, L, λ, γ, ρa, ρp = param
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/$(name)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)"
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/$(name)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).mp4"
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function plot_eta(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Array{Any,1},2})
    @unpack name, L, λ, γ, ρa, ρp, Ω, E, Δt, site_distribution, angles, rates = param
    ax.clear()
    #collect data
    passive, active, directions = extract_points(Ω,η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t), Φ = $(Φ)")
    # Plot points
    ax.quiver(active[1,:],active[2,:], dx,dy,
        directions,
        scale_units = "x",
        pivot = "mid",
        minlength = 0.1,
        minshaft = 1,
        width =1/L,
        headlength = 5,
        headaxislength=5,
        scale = L,
    )
    ax.errorbar(passive[1,:],passive[2,:], 
        markersize = 400/L, 
        fmt= "o", 
        color = "black",
        alpha=0.8,
    )
    return fig
end 

function extract_points(Ω::Array{Array{Int64,1},1}, η::Array{Array{Any,1},2},L)
    passive = []
    active  = []
    directions  = []
    correction = [0.5,0.5]
    for x ∈ Ω
        if η[x...][1]== 2
            push!(passive,x-correction)
        elseif η[x...][1]== 1
            push!(active,x-correction)
            push!(directions,η[x...][2])
        end
    end
    passive = reshape([(passive...)...], (2,:))/L
    active  = reshape([(active...)...], (2,:))/L
    return passive, active, directions
end

##
#phase seperation metircs

function site_ρ(ηx::Array{Any,1})
    if ηx[1] == 1
        return 1
    elseif ηx[1] == 2
        return 1
    else
        return 0
    end
end

function fourier_config(η::Array{Array{Any,1},2}; Ω::Array{Array{Int64,1},1}= [],L::Int64 = 1, k::Vector{Float64}=[0,0])
    ϕLL = 0.
    for x ∈ Ω
        ϕLL += ( 1-site_ρ(η[x...]))*exp(- im* k⋅ x)
    end
    return ϕLL/L^2
end

function translation_invariance(η::Array{Array{Any,1},2}; Ω::Array{Array{Int64,1},1}= [],L::Int64 = 1)
    return norm(fourier_config(η; Ω = Ω, L= L, k =[2*π/L,0.]))+norm(fourier_config(η; Ω = Ω, L= L, k =[0.,2*π/L]))
end



#Example 
#=
#Parameters
param = uniform_initial_param(L=32, λ = 16, ρa = 0.9, ρp = 0.0, Δt = 0.01)
param = extra_mixing_initial_param(L=32, λ = 16, ρa = 0.9, ρp = 0.0, Δt = 0.01, γ = 0.01)
model = initialize_model(param)
#expand variables
@unpack name, L, λ, γ, ρa, ρp, Ω, E, site_distribution, angles, rates = param
@unpack η, c, j, t = model
#Run and save
using BenchmarkTools
T = 0.05
@time model_step!(param, model);
run_and_dump_sim(param,model, T; dump_interval = 0.01, save_on =true)
#Loading
t_saves, η_saves = load_etas(param, T; dump_interval = 0.01)
#plotting
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = length(t_saves)
plot_eta(fig,ax,param, t_saves[n], η_saves[n])
display(fig)
#video
animate_etas(param,t_saves,η_saves)
#symmetry
y = translation_invariance.(η_saves;Ω = Ω,L= L)
n = length(t_saves)
translation_invariance(η_saves[n];Ω = Ω,L= L)
clf()
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(y)
display(fig)
clf()
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = 6000
plot_eta(fig,ax,param, t_saves[n], η_saves[n])
display(fig)
