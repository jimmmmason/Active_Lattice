cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("booted")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random, TensorOperations

function pde_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, Nx = 100, Nθ = 100, δt = 0.0001)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    Eθ = [[1] , [-1]]
    Y = fill([], (Nx,Nx,2))
    for x ∈ Ω 
        for i in 1:2
            local y
            #find adjacent site
            y  = (x + E[i] +[Nx-1,Nx-1]) .% Nx + [1,1]
            #fill jump vectors
            Y[x...,i] = y
        end
    end
    param = Dict{String,Any}()  
    @pack! param = name, D, λ, ρa, ρp, δt, Nx, Nθ, Ω, S,  E, Eθ, Y
    return param
end

function initialize_density(param::Dict{String,Any})
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, Ω, S,  E, Eθ = param
    density = Dict{String,Any}()
    fa = fill(ρa/(2π),(Nx,Nx,Nθ))
    fp = fill(ρp,(Nx,Nx))
    ρ  = fill(ρp+ρa,(Nx,Nx))
    m  = fill(0.,(Nx,Nx,2))
    t = 0.
    @pack! density = fa , fp, ρ, m, t
    return density
end

function U_xyθ(param::Dict{String,Any}, density::Dict{String,Any}; logtol = -1000000)
    @unpack λ, Nx, Nθ, Ω, S,  E, Y= param
    @unpack fa, fp, ρ, m, t = density

    logmρa = map(x -> (x>0 ? log(x) : logtol), fa);
    logmρp = map(x -> (x>0 ? log(x) : logtol), fp);
    logmρ = map(x -> (x>0 ? log(x) : logtol), 1 .- ρ);

    U1 = zeros(Nx,Nx,Nθ, 2);
    U2 = zeros(Nx,Nx,2);
    U3  = zeros(Nx,Nx,2);
    U4  = zeros(Nx,Nx,2);
    Uθ  = zeros(Nx,Nx,Nθ);

    α = π/2 -1
    dₛ = ( -ρ .+1).*( α*(2*α-1)/(2*α+1)*ρ.^2 + α*ρ .+1)

    for x ∈ Ω
        local magx, dsx
        for i in 1:2
            y = Y[x...,i]
            dsx  = (dₛ[x...]+dₛ[y...])/2
            ρx   = (ρ[x...] +ρ[y...])/2
            magx = (m[y...,:] +m[x...,:])/2
            for  θ ∈ S
                eθ = [cos(θ*2π/Nθ), sin(θ*2π/Nθ)]
                U1[x...,θ...,i] = dsx*(- Nx*(  logmρa[y...,θ]) + Nx*( logmρa[x...,θ]) + λ*(eθ⋅E[i]) )
            end 
            U2[x...,i] = dsx*(- Nx*( logmρp[y...] ) + Nx*( logmρp[x...] ))
            if ρx == 0.
            else
                U3[x...,i] = ((1-dsx)/ρx) * (- Nx*( -logmρ[y...]) + Nx*( -logmρ[x...]))
                U4[x...,i] = ((1-ρx-dsx)/ρx) * (λ*E[i] ⋅ magx)
            end
        end
        for θ ∈ S
            ϕ = ((θ % Nθ) +1)
            Uθ[x...,θ] = -Nθ*( logmρa[x...,ϕ] - logmρa[x...,θ] )/(2*π)
        end
    end
    return U1, U2, U3, U4, Uθ
end

function F_xyθ(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack Nx, Nθ, Ω, S,  E, Eθ, Y= param
    @unpack fa, fp, ρ, m, t = density

    Fa = zeros(Nx,Nx,Nθ, 2);
    Fp = zeros(Nx,Nx,2);
    Fθ  = zeros(Nx,Nx,Nθ);

    U1, U2, U3, U4, Uθ = U_xyθ(param, density)
    for x ∈ Ω
        for i in 1:2
            y = Y[x...,i]
            Fa[x...,:,i] = max.(0., U1[x...,:,i] .+U3[x...,i] .+U4[x...,i] ).*fa[x...,:] + min.(0.,U1[x...,:,i]).*fa[y...,:]
            Fp[x...,i]   = max.(0.,U2[x...,i]+U3[x...,i]+U4[x...,i]).*fp[x...] + min.(0.,U2[x...,i]+U3[x...,i]+U4[x...,i]).*fp[y...]
        end
        for θ ∈ S
            ϕ = ((θ % Nθ) +1)
            Fθ[x...,θ] = max.(0.,Uθ[x...,θ]).*fa[x...,θ] + min.(0.,Uθ[x...,θ]).*fa[x...,ϕ]
        end
    end
    return Fa, Fp, Fθ
end


function fast_U_xyθ(param::Dict{String,Any}, density::Dict{String,Any}; logtol = -1000000)
    @unpack λ, Nx, Nθ, Ω, S,  E, Y= param
    @unpack fa, fp, ρ, m, t = density

    logmρa = map(x -> (x>0 ? log(x) : logtol), fa);
    logmρp = map(x -> (x>0 ? log(x) : logtol), fp);
    logmρ = map(x -> (x>0 ? log(x) : logtol), 1 .- ρ);

    U1 = zeros(Nx,Nx,Nθ, 2);
    U2 = zeros(Nx,Nx,2);
    U3  = zeros(Nx,Nx,2);
    U4  = zeros(Nx,Nx,2);
    Uθ  = zeros(Nx,Nx,Nθ);

    α = π/2 -1
    dₛ = ( -ρ .+1).*( α*(2*α-1)/(2*α+1)*ρ.^2 + α*ρ .+1)

    for x₁ in 1:Nx, x₂ in 1:Nx
        local dsx, ρx, magx, y₁,y₂
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        dsx  = (dₛ[x₁,x₂]+dₛ[y₁,y₂])/2
        ρx   = (ρ[x₁,x₂] +ρ[y₁,y₂])/2
        magx = (m[y₁,y₂,:] +m[x₁,x₂,:])/2
        for  θ ∈ S
            U1[x₁,x₂,θ...,2] = dsx*(- Nx*(  logmρa[y₁,y₂,θ]) + Nx*( logmρa[x₁,x₂,θ]) + λ*(sin(θ*2π/Nθ)) )
        end 
            U2[x₁,x₂,2] = dsx*(- Nx*( logmρp[y₁,y₂] ) + Nx*( logmρp[x₁,x₂] ))
        if ρx == 0.
        else
                U3[x₁,x₂,2] = ((1-dsx)/ρx) * (- Nx*( -logmρ[y₁,y₂]) + Nx*( -logmρ[x₁,x₂]))
                U4[x₁,x₂,2] = ((1-ρx-dsx)/ρx) * (λ*E[2] ⋅ magx)
        end
        ## 1 direction 
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        dsx  = (dₛ[x₁,x₂]+dₛ[y₁,y₂])/2
        ρx   = (ρ[x₁,x₂] +ρ[y₁,y₂])/2
        magx = (m[y₁,y₂,:] +m[x₁,x₂,:])/2
        for  θ ∈ S
            U1[x₁,x₂,θ...,1] = dsx*(- Nx*(  logmρa[y₁,y₂,θ]) + Nx*( logmρa[x₁,x₂,θ]) + λ*(cos(θ*2π/Nθ)) )
        end 
            U2[x₁,x₂,1] = dsx*(- Nx*( logmρp[y₁,y₂] ) + Nx*( logmρp[x₁,x₂] ))
        if ρx == 0.
        else
                U3[x₁,x₂,1] = ((1-dsx)/ρx) * (- Nx*( -logmρ[y₁,y₂]) + Nx*( -logmρ[x₁,x₂]))
                U4[x₁,x₂,1] = ((1-ρx-dsx)/ρx) * (λ*E[1] ⋅ magx)
        end
        ## theta diection
        for θ ∈ S
            ϕ = ((θ % Nθ) +1)
            Uθ[x₁,x₂,θ] = -Nθ*( logmρa[x₁,x₂,ϕ] - logmρa[x₁,x₂,θ] )/(2*π)
        end
    end
    return U1, U2, U3, U4, Uθ
end

function fast_F_xyθ(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack Nx, Nθ, Ω, S,  E, Eθ, Y= param
    @unpack fa, fp, ρ, m, t = density

    Fa = zeros(Nx,Nx,Nθ, 2);
    Fp = zeros(Nx,Nx,2);
    Fθ  = zeros(Nx,Nx,Nθ);

    U1, U2, U3, U4, Uθ = fast_U_xyθ(param, density)
    for x₁ in 1:Nx, x₂ in 1:Nx
        local y₁,y₂
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        Fa[x₁,x₂,:,2] = max.(0., U1[x₁,x₂,:,2] .+U3[x₁,x₂,2] .+U4[x₁,x₂,2] ).*fa[x₁,x₂,:] + min.(0.,U1[x₁,x₂,:,2]).*fa[y₁,y₂,:]
        Fp[x₁,x₂,2]   = max.(0.,U2[x₁,x₂,2]+U3[x₁,x₂,2]+U4[x₁,x₂,2]).*fp[x₁,x₂] + min.(0.,U2[x₁,x₂,2]+U3[x₁,x₂,2]+U4[x₁,x₂,2]).*fp[y₁,y₂]
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        Fa[x₁,x₂,:,1] = max.(0., U1[x₁,x₂,:,1] .+U3[x₁,x₂,1] .+U4[x₁,x₂,1] ).*fa[x₁,x₂,:] + min.(0.,U1[x₁,x₂,:,1]).*fa[y₁,y₂,:]
        Fp[x₁,x₂,1]   = max.(0.,U2[x₁,x₂,1]+U3[x₁,x₂,1]+U4[x₁,x₂,1]).*fp[x₁,x₂] + min.(0.,U2[x₁,x₂,1]+U3[x₁,x₂,1]+U4[x₁,x₂,1]).*fp[y₁,y₂]
        for θ ∈ S
            ϕ = ((θ % Nθ) +1)
            Fθ[x₁,x₂,θ] = max.(0.,Uθ[x₁,x₂,θ]).*fa[x₁,x₂,θ] + min.(0.,Uθ[x₁,x₂,θ]).*fa[x₁,x₂,ϕ]
        end
    end
    return Fa, Fp, Fθ
end

function pde_step!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ, Ω, S, Y = param
    @unpack fa, fp, ρ, m, t = density
    Fa, Fp, Fθ = fast_F_xyθ(param, density);
    for x₁ in 1:Nx, x₂ in 1:Nx
        local y₁,y₂
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        fa[y₁,y₂,:] += -δt*Nx*(Fa[y₁,y₂,:,2]-Fa[x₁,x₂,:,2])
        fp[y₁,y₂]   += -δt*Nx*(Fp[y₁,y₂,2]  -Fp[x₁,x₂,2]  )
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        fa[y₁,y₂,:] += -δt*Nx*(Fa[y₁,y₂,:,1]-Fa[x₁,x₂,:,1])
        fp[y₁,y₂]   += -δt*Nx*(Fp[y₁,y₂,1]  -Fp[x₁,x₂,1]  )
        for θ ∈ S
            ϕ = ((θ % Nθ) +1)
            fa[x₁,x₂,ϕ] += -δt*Nθ*(Fθ[x₁,x₂,ϕ]-Fθ[x₁,x₂,θ])/(2*π)
        end
    end
    t += δt
    ρ  = fp + sum(fa; dims = 3)[:,:,1]*2*π/Nθ
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := fa[a,b,θ]*eθ[θ,d]
    end
    @pack! density = fa, fp, ρ, m, t
end

function run_pde_until!(param::Dict{String,Any},density::Dict{String,Any},T; save_on =false, max_steps = 100)
    if max_steps == false
        while density["t"] < T 
            pde_step!(param,density)
            if save_on
                @unpack name, L, λ, ρa, ρp, Δt = param
                @unpack η, t = model
                filename = "/home/jm2386/Active_Lattice/data/pde_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
                data = Dict{String,Any}();
                @pack! data = param, η, t
                safesave(filename,data)
            end
        end
    else
        for i in 1:max_steps
            if density["t"] < T
                pde_step!(param,density)
            end
            if save_on
                @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
                @unpack fa, fp, ρ, m, t = density
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
                data = Dict{String,Any}();
                @pack! data = param, fa, fp, ρ, m, t
                safesave(filename,data)
            end
        end
    end
end

function run_and_dump_pde(param::Dict{String,Any},model::Dict{String,Any},T; dump_interval = 0.001, save_on =false)
    while model["t"] < T
        local t
        t = deepcopy(model["t"]+dump_interval)
        run_pde_until!(param,density,t; save_on =false, max_steps = false)
        if save_on
            @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
            @unpack fa, fp, ρ, m, t = density
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
            data = Dict{String,Any}();
            @pack! data = param, fa, fp, ρ, m, t
            safesave(filename,data)
        end
    end
end

#=
param = pde_param()
density = initialize_density(param)
@time pde_step!(param,density);
@time run_pde_until!(param,density,1.; max_steps = 100);
density
run_and_dump_pde(param,density,0.0006; dump_interval = 0.0002, save_on =true)

@time for x∈Ω
    y = Y[x...,1]
    ρ[x...], ρ[y...] = ρ[y...], ρ[x...]
end

@time for i in 1:100, j in 1:100
    k = (i +100)%100 +1
    ρ[k,j], ρ[i,j] = ρ[i,j], ρ[k,j]
    #println(j)
end

=#

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
param = pde_param()
density = initialize_density(param)
#expand variables
@unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, Ω, S,  E, Eθ = param
@unpack fa, fp, ρ, t = density
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
