cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("booted")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random, TensorOperations, StaticArrays

function pde_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, Nx = 100, Nθ = 100, δt = 0.0001)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    param = Dict{String,Any}()
    @pack! param = name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E
    return param
end

function initialize_density(param::Dict{String,Any})
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    density = Dict{String,Any}()
    fa = fill(ρa/(2π),(Nx,Nx,Nθ))
    fp = fill(ρp,(Nx,Nx))
    t = 0.
    @pack! density = fa , fp, t
    return density
end

function fast_U_xyθ(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, E = param
    logtol::Float64 = log(1e-10);
    Ua = Array{Float64,4}(undef, Nx, Nx, Nθ,2);
    Up = Array{Float64,3}(undef, Nx, Nx,2);
    Uθ = Array{Float64,3}(undef, Nx, Nx, Nθ);
    logmρa = Array{Float64,3}(undef, Nx, Nx, Nθ);
    logmρp = Array{Float64,2}(undef, Nx, Nx);
    logmρp = Array{Float64,2}(undef, Nx, Nx);
    ρ = Array{Float64,2}(undef, Nx, Nx);
    m = Array{Float64,3}(undef, Nx, Nx, 2);
    eθ = Array{Float64,2}(undef, Nθ, 2);
    @unpack fa, fp = density

    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := fa[a,b,θ]*eθ[θ,d]
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)

    logmρa = map(x -> (x>0 ? log(x) : logtol), fa);
    logmρp = map(x -> (x>0 ? log(x) : logtol), fp);
    logmρ = map(x -> (x>0 ? log(x) : logtol), 1 .- ρ);

    Ua  = zeros(Nx,Nx,Nθ, 2);
    Up  = zeros(Nx,Nx,2);
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
        if ρx == 0.
        else
            Ua[x₁,x₂,:,2] = dsx*(- Nx*(  logmρa[y₁,y₂,:])+Nx*(logmρa[x₁,x₂,:]) + λ*(eθ[:,1]) ) .+ (((1-dsx)/ρx) * (- Nx*( -logmρ[y₁,y₂]) + Nx*( -logmρ[x₁,x₂]))+((1-ρx-dsx)/ρx) * (λ*E[2] ⋅ magx))
            Up[x₁,x₂,2] = dsx*(- Nx*( logmρp[y₁,y₂] ) + Nx*( logmρp[x₁,x₂] )) + (((1-dsx)/ρx) * (- Nx*( -logmρ[y₁,y₂]) + Nx*( -logmρ[x₁,x₂]))+((1-ρx-dsx)/ρx) * (λ*E[2] ⋅ magx))
        end
        ## 1 direction 
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        dsx  = (dₛ[x₁,x₂]+dₛ[y₁,y₂])/2
        ρx   = (ρ[x₁,x₂] +ρ[y₁,y₂])/2
        magx = (m[y₁,y₂,:] + m[x₁,x₂,:])/2
        if ρx == 0.
        else
            Ua[x₁,x₂,:,1] = dsx*(- Nx*(  logmρa[y₁,y₂,:]) + Nx*( logmρa[x₁,x₂,:]) + λ*(eθ[:,2]) ) .+ ( ((1-ρx-dsx)/ρx) * (λ*E[1] ⋅ magx) +((1-dsx)/ρx) * (- Nx*( -logmρ[y₁,y₂]) + Nx*( -logmρ[x₁,x₂])) )
        #end 
            Up[x₁,x₂,1] = dsx*(- Nx*( logmρp[y₁,y₂] ) + Nx*( logmρp[x₁,x₂] ))
        end
    end
    ## theta diection
    for θ ∈ S
        ϕ = ((θ % Nθ) +1)
        Uθ[:,:,θ] = -Nθ*( logmρa[:,:,ϕ] - logmρa[:,:,θ] )/(2*π)
    end
    return Ua, Up, Uθ
end

function fast_F_xyθ(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack Nx, Nθ, S= param
    Ua = Array{Float64,4}(undef, Nx, Nx, Nθ,2);
    Up = Array{Float64,3}(undef, Nx, Nx,2);
    Uθ = Array{Float64,3}(undef, Nx, Nx, Nθ);
    Fa = Array{Float64,4}(undef, Nx, Nx, Nθ,2);
    Fp = Array{Float64,3}(undef, Nx, Nx,2);
    Fθ = Array{Float64,3}(undef, Nx, Nx, Nθ);
    fa = Array{Float64,3}(undef, Nx, Nx, Nθ);
    fp = Array{Float64,2}(undef, Nx, Nx);
    @unpack fa, fp = density

    Fa = zeros(Nx,Nx,Nθ, 2);
    Fp = zeros(Nx,Nx,2);
    Fθ  = zeros(Nx,Nx,Nθ);

    Ua, Up, Uθ = fast_U_xyθ(param, density)
    for x₁ in 1:Nx, x₂ in 1:Nx
        local y₁,y₂
        ## 2 direction
        y₁::Int64 = x₁
        y₂::Int64= (x₂ +Nx)%Nx +1
        Fa[x₁,x₂,:,2] = max.(0., Ua[x₁,x₂,:,2] ).*fa[x₁,x₂,:] + min.(0.,Ua[x₁,x₂,:,2]).*fa[y₁,y₂,:]
        Fp[x₁,x₂,2]   = max.(0.,Up[x₁,x₂,2]).*fp[x₁,x₂] + min.(0.,Up[x₁,x₂,2]).*fp[y₁,y₂]
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        Fa[x₁,x₂,:,1] = max.(0., Ua[x₁,x₂,:,1] ).*fa[x₁,x₂,:] + min.(0.,Ua[x₁,x₂,:,1]).*fa[y₁,y₂,:]
        Fp[x₁,x₂,1]   = max.(0.,Up[x₁,x₂,1]).*fp[x₁,x₂] + min.(0.,Up[x₁,x₂,1]).*fp[y₁,y₂]
    end
    for θ ∈ S
        local ϕ 
        ϕ::Int64 = ((θ % Nθ) +1)
        Fθ[:,:,θ] = max.(0.,Uθ[:,:,θ]).*fa[:,:,θ] + min.(0.,Uθ[:,:,θ]).*fa[:,:,ϕ]
    end
    a::Float64 = maximum(abs.(Ua));
    b::Float64= maximum(abs.(Up));
    c::Float64 = maximum(abs.(Uθ));
    tempu::Float64 = 1/(6*max(a*Nx, b*Nx, c*Nθ));
    return Fa, Fp, Fθ, tempu
end

function pde_step!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ = param
    Fa = Array{Float64,4}(undef, Nx, Nx, Nθ,2);
    Fp = Array{Float64,3}(undef, Nx, Nx,2);
    Fθ = Array{Float64,3}(undef, Nx, Nx, Nθ);
    fa = Array{Float64,3}(undef, Nx, Nx, Nθ);
    fp = Array{Float64,2}(undef, Nx, Nx);
    @unpack fa, fp, t = density
    Fa, Fp, Fθ, tempu = fast_F_xyθ(param, density);
    dt::Float64= min(δt, tempu)
    for x₁ in 1:Nx, x₂ in 1:Nx
        local y₁,y₂
        ## 2 direction
        y₁ ::Int64 = x₁
        y₂ ::Int64 = (x₂ +Nx)%Nx +1
        fa[y₁,y₂,:] += -dt*Nx*(Fa[y₁,y₂,:,2]-Fa[x₁,x₂,:,2])
        fp[y₁,y₂]   += -dt*Nx*(Fp[y₁,y₂,2]  -Fp[x₁,x₂,2]  )
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        fa[y₁,y₂,:] += -dt*Nx*(Fa[y₁,y₂,:,1]-Fa[x₁,x₂,:,1])
        fp[y₁,y₂]   += -dt*Nx*(Fp[y₁,y₂,1]  -Fp[x₁,x₂,1]  )
    end
    for θ in 1:Nθ
        ϕ ::Int64 = ((θ % Nθ) +1)
        fa[:,:,ϕ] += -dt*Nθ*(Fθ[:,:,ϕ]-Fθ[:,:,θ])/(2*π)
    end
    t::Float64 += dt
    #println(t)
    @pack! density = fa, fp, t
    #println(density["t"])
    return dt
end

##

function run_pde_until!(param::Dict{String,Any},density::Dict{String,Any},T; save_on =false, max_steps = 100, save_interval = 1.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    if max_steps == false
        while density["t"] < T 
                time_since_save = 0.
                while time_since_save < min(save_interval, T)
                    dt = pde_step!(param,density)
                    time_since_save += dt 
                end
                if save_on
                    @unpack fa, fp, t = density
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
                    data = Dict{String,Any}();
                    @pack! data = fa, fp, t
                    safesave(filename,data)
                end
            end
    else
        steps ::Int64 = 0
        while (density["t"] < T)&&(steps < max_steps)
            time_since_save = 0.
            while time_since_save < min(save_interval, T)
                dt = pde_step!(param,density)
                time_since_save += dt 
            end
            if save_on
                @unpack fa, fp, t = density
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
                data = Dict{String,Any}();
                @pack! data = fa, fp, t
                safesave(filename,data)
            end
            steps += 1
        end
        return steps
    end
end

function perturb_pde!(param::Dict{String,Any}, density::Dict{String,Any}; δ = 0.01, n2 = 2)
    @unpack Nx, S, ρa, ρp, λ = param
    @unpack fa = density
    δ = min(δ, 1 - ρa - ρp);
    #from maria paper: 
    β = ρa*λ/2;
    γ = ρa;
    Λ = 1/2*(-1 - 4*γ*n2*π^2 + sqrt(1 + 8*n2*π^2*(β^2 + γ*(2*γ*n2*π^2 - 1))));
    C = 2*π*β / (1 + Λ)
    P = (x,y,θ) -> δ*(cos(2*π*x/Nx)*cos(2*π*y/Nx) - C*(sin(2*π*x/Nx)*cos(2*π*y/Nx)*cos(2π*θ/Nθ) + sin(2*π*y/Nx)*cos(2*π*x/Nx)*sin(2π*θ/Nθ) ));
    # 
    for x₁ in 1:Nx, x₂ in 1:Nx, θ in S
        fa[x₁, x₂, θ] = P(x₁, x₂, θ);
    end
    @pack! density = fa;
end

#=
param = pde_param()
density = initialize_density(param)
perturb_pde!(param,density);
run_pde_until!(param,density,0.001; save_on = true, max_steps = 20, save_interval = 0.0001)
pde_step!(param, density);
density["t"]

P = (x,y,θ) -> cos(2*pi*x)*cos(2*pi*y) - (sin(2*pi*x)*cos(2*pi*y)*cos(θ) + sin(2*pi*y)*cos(2*pi*x)*sin(θ) ));
=#
#visualise data
using PyPlot, PyCall
@pyimport matplotlib.animation as anim

function load_pdes(param::Dict{String,Any},T; save_interval = 1.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    t = save_interval
    fa_saves = []
    fp_saves = []
    t_saves = []
    while t<T
        try 
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).jld2";
            data = wload(filename)
            push!(t_saves,data["t"])
            push!(fa_saves,data["fa"])
            push!(fp_saves,data["fp"])
        catch
        end
        t += save_interval
    end
    return t_saves, fa_saves, fp_saves
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

function plot_pde(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S= param
    ax.clear()
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := 2π *fa[a,b,θ]*eθ[θ,d]/Nθ
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    cp = ax.contourf(absmag;levels = 100, set_cmap = winter() )
    fig.colorbar(cp)
    display(fig)
    return fig


    directions[2]
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
#testing
param = pde_param()
density = initialize_density(param)

@profview fast_U_xyθ(param, density);
@profview fast_F_xyθ(param, density);
@profview pde_step!(param, density);

@time fast_U_xyθ(param, density);
@time fast_F_xyθ(param, density);
@time pde_step!(param, density)

@unpack λ, Nx, Nθ, Ω, S,  E, Y= param
@unpack fa, fp, ρ, m, t = density
@time afast_U(param, density);
@profview tfast_U(fa, fp, ρ, m, t, λ, Nx, Nθ);
=#
