cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

#Runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random

function uniform_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, d=2)
    param = Dict{String,Any}()  
    #this is the only dimension dependent part:
    Ω = [[i,j] for i in 1:L for j in 1:L] 
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))
    Δt = 0.001/L^2
    function angles(x,n) 
        if n == 1
            return 2*π*rand()
        else
            return -1
        end
    end
    function rates(n,m,i)
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
    @pack! param = name, L, D, λ, ρa, ρp, Δt, Ω, E, site_distribution, angles, rates
    return param
end

function initialize(param::Dict{String,Any})
    @unpack name, L, λ, ρa, ρp, Ω, E, site_distribution, angles, rates = param
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
    t = 0.
    model = Dict{String,Any}() 
    @pack! model = η, c, j, t
    return model
end

function model_step!(param::Dict{String,Any},model::Dict{String,Any})
    @unpack name, L, λ, ρa, ρp, Ω, E, Δt, site_distribution, angles, rates = param
    @unpack η, c, j, t = model
    # increae timestep
    t += Δt
    # see if jump occurs
    if rand() < Δt*sum(c)
        #select jump
        w     = Weights( [(c...)...])
        jump  = sample(j, w)
        #execute jumps
        η[jump[2]...] = deepcopy(η[jump[1]...])
        η[jump[1]...] = [0,-1]
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
    end
    #diffuse angles
    for x ∈ Ω
        if η[x...][1] == 1
            η[x...][2] = (η[x...][2] + sqrt(Δt)*randn() + 2*π) % (2*π)
        end
    end
    @pack! model = η, c, j, t
    return t
end

function run_model_until!(param::Dict{String,Any},model::Dict{String,Any},T; return_all = false, save_on=true)
    if return_all
        η_saves = []
        t_saves = []
        start_time = deepcopy(model["t"])
        while model["t"] < T 
            model_step!(param,model)
            push!(η_saves,model["η"])
            push!(t_saves,model["t"])
        end
        if save_on
            @unpack name, L, λ, ρa, ρp = param
            filename = "/home/jm2386/Active_Lattice/data/sims_raw/$(name)/start_time=$(start_time)_end_time=$(T)_interval=$(param["Δt"])_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
            data = Dict{String,Any}();
            @pack! data = param, η_saves, t_saves
            safesave(filename,data)
        end
        return t_saves, η_saves
    else
        while model["t"] < T 
            model_step!(param,model)
        end
    end
end

function run_model_intervals!(param::Dict{String,Any},model::Dict{String,Any},T; interval = 0.001, save_on =false)
    η_saves = []
    t_saves = []
    start_time = deepcopy(model["t"])
    while model["t"] < T
        local t
        t = deepcopy(model["t"]+interval)
        run_model_until!(param,model,t; return_all = false)
        push!(η_saves, model["η"])
        push!(t_saves, model["t"])
    end 
    if save_on
        @unpack name, L, λ, ρa, ρp = param
        filename = "/home/jm2386/Active_Lattice/data/sims_raw/$(name)/start_time=$(start_time)_end_time=$(T)_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
        data = Dict{String,Any}();
        @pack! data = param, η_saves, t_saves
        safesave(filename,data)
    end
    return t_saves, η_saves
end

#visualise data
using PyPlot

function animate_etas(param,η_saves,t_saves)
    frames = length(η_saves)
    makeframe(i) = plot_eta(param, η_saves[i],t_saves[i])
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=20)

    # Convert it to an MP4 movie file and saved on disk in this format.
    filename = "testvid.mp4"
    myanim[:save](filename, bitrate=-1, dpi= 500, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function plot_eta(param::Dict{String,Any}, η::Array{Array{Any,1},2},t::Float64)
    @unpack name, L, λ, ρa, ρp, Ω, E, Δt, site_distribution, angles, rates = param
    #collect data
    passive, active, directions = extract_points(Ω,η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    #figure configuration
    fig, ax = PyPlot.subplots(figsize =(10, 10))
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
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
        markersize = 3, 
        fmt= "o", 
        color = "black",
        alpha=0.8,
    )
    fig
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

#Example 
#=
#Parameters
param = uniform_initial_param(L=100)
model = initialize(param)
#expand variables
@unpack name, L, λ, ρa, ρp, Ω, E, site_distribution, angles, rates = param
@unpack η, c, j, t = model
#Run and save
T = 0.001
#model_step!(param,model)
#η_saves,t_saves = run_model_until!(param,model,T; return_all = true, save_on =true)
η_saves,t_saves = run_model_intervals!(param,model,T; interval = T/100, save_on =true)
#Loading
@unpack name, L, λ, ρa, ρp = param
filename = "/home/jm2386/Active_Lattice/data/sims_raw/$(name)/time=$(T)_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
data = wload(filename)
#plotting
fig = plot_eta(param, η, 0.)
display(fig)
=#