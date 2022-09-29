cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("booted")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random

function uniform_pde_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, Nx = 1000, Nθ = 1000, δt = 0.001)
    param = Dict{String,Any}()  
    #this is the only dimension dependent part:
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ [θ] for Nθ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    Eθ = [[1] , [-1]]
    @pack! param = name, D, λ, ρa, ρp, Δt, Ω, S  E, Eθ
    return param
end

function initialize_pde(param::Dict{String,Any})
    @unpack name, D, λ, ρa, ρp, Δt, Ω, S  E, Eθ = param
    model = Dict{String,Any}()
    f = [fill(ρa,(Nx,Nx,θ)),fill( ρp,(Nx,Nx))] 


    return model
end

function U_xyθ(param::Dict{String,Any}, pde)
end

function mobility(param::Dict{String,Any}, pde)
end

function F_xyθ(param::Dict{String,Any}, pde)
end

function pde_step(param::Dict{String,Any}, pde)
end

#visualise data
using PyPlot, PyCall
@pyimport matplotlib.animation as anim

function load_pdes(param::Dict{String,Any},T)
    @unpack name, L, λ, ρa, ρp, Δt = param
    t = Δt
    η_saves = []
    t_saves = []
    while t<T
        try 
            filename = "/home/jm2386/Active_Lattice/data/pde_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
            data = wload(filename)
            push!(t_saves,data["t"])
            push!(η_saves,data["η"])
        catch
        end
        t += Δt
    end
    return t_saves, pde_saves
end

function animate_pdes(param,t_saves,pde_saves)
    @unpack name, L, λ, ρa, ρp, Δt = param
    frames = length(η_saves)-1
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    makeframe(i) = plot_eta(fig,ax,param, t_saves[i+1], η_saves[i+1])
    interval = Int64(round(5/Δt))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    pathname = "/home/jm2386/Active_Lattice/plots/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)";
    mkpath(pathname)
    filename = "/home/jm2386/Active_Lattice/plots/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(T; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function plot_pde(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Array{Any,1},2})
    @unpack name, L, λ, ρa, ρp, Ω, E, Δt, site_distribution, angles, rates = param
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

function fourier_pde(η::Array{Array{Any,1},2}; Ω::Array{Array{Int64,1},1}= [],L::Int64 = 1, k::Vector{Float64}=[0,0])
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
param = param = uniform_initial_param(; name = name, λ = 16 ,ρa = 0.95, ρp = 0., L=64, Δt = 0.01)
model = initialize_model(param)
#expand variables
@unpack name, L, λ, ρa, ρp, Ω, E, Δt, site_distribution, angles, rates = param
@unpack η, c, j, t = model
#Run and save
using BenchmarkTools
T = 0.08
@time model_step!(param,model);
@time run_model_until!(param, model, T; save_on = true);
#Loading
T = 1.0
t_saves, η_saves = load_etas(param, T)
#plotting
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_eta(fig,ax,param, t, η)
display(fig)
# plot individual frame
n = 6000
n = length(t_saves)
clf()
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_eta(fig,ax,param, t_saves[n], η_saves[n])
display(fig)
#create video
animate_etas(param,t_saves,η_saves)
#plot symmetry over time
y = translation_invariance.(η_saves;Ω = Ω,L= L)
n = length(t_saves)
translation_invariance(η_saves[n];Ω = Ω,L= L)
clf()
y = translation_invariance.(η_saves;Ω = Ω,L= L)
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(y)
display(fig)
=#
