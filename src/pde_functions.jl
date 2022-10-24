cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("booted")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random, TensorOperations, StaticArrays

function pde_param(; name = "test", D =1. , λ =1. ,ρa = 0.5, ρp = 0.0, Nx = 100, Nθ = 100, δt = 1e-5, Dθ = 10, T= 0.001)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    param = Dict{String,Any}()
    @pack! param = name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E, Dθ, T
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

## 

function midpoint_bond_diff(f::Array{Float64,2}; Nx::Int64 = 100) 

    grad_f::Array{Float64,3} = zeros(2, Nx, Nx)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        grad_f[1,x₁,x₂] = Nx*( f[y₁,y₂] - f[x₁,x₂] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        grad_f[2,x₁,x₂] = Nx*( f[y₁,y₂] - f[x₁,x₂] ) 
    end
    return grad_f
end

function midpoint_bond_diff_θ(f::Array{Float64,3}; Nx::Int64 = 100,  Nθ::Int64 = 100) 

    grad_f::Array{Float64,4} = zeros(2, Nx, Nx, Nθ)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        grad_f[1,x₁,x₂,:] = Nx*( f[y₁,y₂,:] - f[x₁,x₂,:] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        grad_f[2,x₁,x₂,:] = Nx*( f[y₁,y₂,:] - f[x₁,x₂,:] ) 
    end
    return grad_f
end

function site_div(f::Array{Float64,3}; Nx::Int64 = 100) 

    div_f::Array{Float64,2} = zeros(Nx, Nx)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        div_f[x₁,x₂] = Nx*( f[1,y₁,y₂] - f[1,x₁,x₂] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        div_f[x₁,x₂] = Nx*( f[2,y₁,y₂] - f[2,x₁,x₂] ) 
    end
    return div_f
end

function site_div_θ(f::Array{Float64,4}; Nx::Int64 = 100,  Nθ::Int64 = 100) 

    div_f::Array{Float64,3} = zeros(Nx, Nx, Nθ)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        div_f[y₁,y₂,:] = Nx*( f[1,y₁,y₂,:] - f[1,x₁,x₂,:] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        div_f[y₁,y₂,:] = Nx*( f[2,y₁,y₂,:] - f[2,x₁,x₂,:] ) 
    end

    return div_f
end

function site_θ_diff(f::Array{Float64,3}; Nx::Int64 = 100,  Nθ::Int64 = 100) 

    div_f::Array{Float64,3} = zeros(Nx, Nx, Nθ)

    for θ in 1:Nθ
        ϕ ::Int64 = ((θ % Nθ) +1)
        div_f[:,:,ϕ] += Nθ*(f[:,:,ϕ]-f[:,:,θ])/(2*π)
    end

    return div_f
end

function midpoint_bond_av(f::Array{Float64,3}; dims =2, Nx::Int64 = 100) 
    av_f::Array{Float64,3}= zeros(2, Nx, Nx)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        av_f[1,x₁,x₂] = ( f[1,y₁,y₂] + f[1,x₁,x₂] )/2
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        av_f[2,x₁,x₂] = ( f[2,y₁,y₂] + f[2,x₁,x₂] )/2
    end
    return av_f
end

function midpoint_Θ_diff(f::Array{Float64,3}; Nx::Int64 = 100,  Nθ::Int64 = 100) 
    grad_f::Array{Float64,3} = zeros(Nx, Nx, Nθ)

    for θ in 1:Nθ
        ϕ = ((θ % Nθ) +1)
        grad_f[:,:,θ] = Nθ*( f[:,:,ϕ] - f[:,:,θ] )/(2*π)
    end
    return grad_f
end

##

function self_diff(ρ::Float64)
    α::Float64= π/2 -1;
    return ( 1-ρ).*( α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ +1)
end

function mag(f::Array{Float64,3}; Nθ = 50, Nx =100)
    eθ::Array{Float64,2} = [cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)];
    m = Array{Float64,3}(undef, 2, Nx, Nx);
    @tensor begin
        m[d,a,b] := f[a,b,θ]*eθ[θ,d]
    end
    return 2*π*m/Nθ
end

function coeff_s(rho::Float64,ds::Float64)
    return (rho>0 ? (1-rho-ds)/(rho*ds) : 0)
end

function coeff_mag_s(f::Array{Float64,3},ρ::Array{Float64,2}; Nθ::Int64 = 100,  Nx::Int64 = 100)
    m    ::Array{Float64,3} = mag(f; Nθ=Nθ, Nx=Nx );
    ds   ::Array{Float64,2} = self_diff.(ρ);
    s    ::Array{Float64,3} = reshape(coeff_s.(ρ,ds),1,Nx,Nx);
    mag_s::Array{Float64,3} = s.*m
    return mag_s
end

function p(x::Float64)
    return -(π-1)*log(1-x)   + real(  -   ( (4 -5*π +π^2)*sqrt( Complex( (π-2)/( -26 +37*π -12π^2 +π^3 )) )*atanh( sqrt( Complex( (π-2)/( -26 +37*π -12π^2 +π^3 )) )*( 1 -6*x + π*(-1+2*x) )))   +    (1/2)*(-2 + π)*log( -2 -2*x + π^2(-1+x)*x + 6*x^2 + π*(2 +3*x -5x^2))     )
end

function mob(fa::Array{Float64,3}, fp::Array{Float64,2}, ρ::Array{Float64,2})
    ds::Array{Float64,2} = self_diff.(ρ)
    return fa.*ds, fp*ds, fa
end

function upwind(U::Float64, mb_down::Float64, mb_up::Float64)
    return (U   > 0. ? U.*mb_down  : U*mb_up)
end

##

function U_apθ(fa::Array{Float64,3}, fp::Array{Float64,2}, ρ::Array{Float64,2}; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10.)
    logtol::Float64 = log(1e-10);

    eθ:: Array{Float64,4} = reshape([cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)],2,1,1,Nθ)

    logmfa::Array{Float64,3} = map(x -> (x>0 ? log(x) : logtol), fa);
    logmfp::Array{Float64,2} = map(x -> (x>0 ? log(x) : logtol), fp);
    p_rho ::Array{Float64,2} = p.(ρ)

    Ua::Array{Float64,4}  = -midpoint_bond_diff_θ(logmfa .+ p_rho; Nx=Nx, Nθ=Nθ) .+ λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ) .+ λ*eθ 
    Up::Array{Float64,3}  = -midpoint_bond_diff(  logmfp  + p_rho; Nx=Nx       )  + λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx )
    Uθ::Array{Float64,3}  = -midpoint_Θ_diff(fa; Nx=Nx, Nθ = Nθ)

    return Ua, Up, Uθ
end

function F_apθ(Ua::Array{Float64,4}, Up::Array{Float64,3}, Uθ::Array{Float64,3}, moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3}; Nx::Int64 =100, Nθ::Int64 =100 )
    Fa ::Array{Float64,4} = zeros(2,Nx,Nx,Nθ);
    Fp ::Array{Float64,3} = zeros(2,Nx,Nx);
    Fθ ::Array{Float64,3} = zeros(Nx,Nx,Nθ);
    for x₁ in 1:Nx, x₂ in 1:Nx
        local y₁,y₂
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        Fa[1,x₁,x₂,:] = upwind.(Ua[1,x₁,x₂,:], moba[x₁,x₂,:], moba[y₁,y₂,:])
        Fp[1,x₁,x₂]   = upwind( Up[1,x₁,x₂]  , mobp[x₁,x₂],   mobp[y₁,y₂]  )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        Fa[2,x₁,x₂,:] = upwind.(Ua[2,x₁,x₂,:], moba[x₁,x₂,:], moba[y₁,y₂,:])
        Fp[2,x₁,x₂]   = upwind( Up[2,x₁,x₂]  , mobp[x₁,x₂],   mobp[y₁,y₂]  )
    end
    for θ in 1:Nθ
        local ϕ 
        ϕ::Int64 = ((θ % Nθ) +1)
        Fθ[:,:,θ] = upwind.(Uθ[:,:,θ], mobθ[:,:,θ], mobθ[:,:,ϕ])
    end
    return Fa, Fp, Fθ
end

##

function time_step(fa::Array{Float64,3}, fp::Array{Float64,2}, δt::Float64; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., Dθ::Float64 = 10.)
    ρ::Array{Float64,2} = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    
    Ua::Array{Float64,4},   Up::Array{Float64,3},   Uθ::Array{Float64,3}   = U_apθ(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ)
    moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3} = mob(fa,fp,ρ)
    Fa::Array{Float64,4},   Fp::Array{Float64,3},   Fθ::Array{Float64,3}   = F_apθ(Ua, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)
    
    a::Float64 = maximum(abs.(Ua));
    b::Float64 = maximum(abs.(Up));
    c::Float64 = maximum(abs.(Uθ));
    
    tempu::Float64 = 1/(6*max(a*Nx, b*Nx, c*Nθ*Dθ/(2*π)));
    dt::Float64= min(δt, tempu)

    fa -= dt*( site_div_θ(Fa; Nx=Nx, Nθ=Nθ) + Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ))
    fp -= dt*site_div(Fp; Nx=Nx)

    return fa, fp, dt
end

function pde_step!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ, Dθ, λ = param
    @unpack fa, fp, t = density
    fa, fp , dt = time_step(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)
    t::Float64 += dt

    @pack! density = fa, fp, t
    return dt
end

##

function run_pde_until!(param::Dict{String,Any},density::Dict{String,Any},T; save_on =false, max_steps = 100, save_interval = 1.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    if save_on
        @unpack t = density
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
        safesave(filename,density)
    end
    if max_steps == false
        while density["t"] < T 
                time_since_save = 0.
                while time_since_save < min(save_interval, T)
                    dt = pde_step!(param,density)
                    time_since_save += dt 
                end
                if save_on
                    @unpack t = density
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
                    safesave(filename,density)
                end
            end
    else
        steps ::Int64 = 0
        while (density["t"] < T)&(steps < max_steps)
            time_since_save = 0.
            while (time_since_save < min(save_interval, T))&(steps < max_steps)
                dt = pde_step!(param,density)
                time_since_save += dt 
                steps += 1
            end
            if save_on
                @unpack t = density
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
                safesave(filename,density)
            end
        end
        return steps
    end
end

function perturb_pde!(param::Dict{String,Any}, density::Dict{String,Any}; δ = 0.01, n2 = 2)
    @unpack Nx, S, ρa, ρp, λ, Dθ, Nx, Nθ = param
    @unpack fa = density
    δ = min(δ, 1 - ρa - ρp);
    #from stability: 
    ρ = ρa + ρp
    α = π/2 -1
    ds = ( -ρ .+1).*( α*(2*α-1)/(2*α+1)*ρ.^2 + α*ρ .+1)
    dsp = - ( α*(2*α-1)/(2*α+1)*ρ.^2 + α*ρ .+1) + ( -ρ .+1)*(2*α*(2*α-1)/(2*α+1)*ρ + α );
    β = 1- ds;
    a = λ* ( (1-ρ-ds)/ρ - dsp) /2
    γ = -8*π^2
    Λ = 1/2*(-Dθ + γ*β + sqrt(Dθ^2 + 2*Dθ*β*γ + (β*γ)^2 - 2*a^2*γ*ρa*ρ ) )
    C = 2*π*ρa*(ρa + ρp)/ (1 + Λ)
    P = (x,y,θ) -> δ*(cos(2*π*x/Nx)*cos(2*π*y/Nx) - C*(sin(2*π*x/Nx)*cos(2*π*y/Nx)*cos(2π*θ/Nθ) + sin(2*π*y/Nx)*cos(2*π*x/Nx)*sin(2π*θ/Nθ) ));
    # 
    for x₁ in 1:Nx, x₂ in 1:Nx, θ in S
        fa[x₁, x₂, θ] += P(x₁, x₂, θ);
    end
    @pack! density = fa;
end

function perturb_pde_run(param; max_steps = 1e8, save_interval = 0.0001)
    @unpack T = param
    density = initialize_density(param)
    perturb_pde!(param,density);
    run_pde_until!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)
end

#visualise data
using PyPlot, PyCall
@pyimport matplotlib.animation as anim

function load_pdes(param::Dict{String,Any},T; save_interval = 1.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    t = 0.
    fa_saves = []
    fp_saves = []
    t_saves = []
    while t ≤ T
        try 
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt).jld2";
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
    myanim[:save](filename, bitrate=-1, dπ= 100, extra_args=["-vcodec", "libx264", "-πx_fmt", "yuv420p"])
end 

function plot_pde_mag(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #fig, ax = PyPlot.subplots(figsize =(10, 10))
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := 2π *fa[a,b,θ]*eθ[θ,d]/Nθ
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
        #ax.yaxis.set_ticks([])
        ax.axis([Δx, 1., Δx, 1.])
        ax.set_aspect("equal")
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    cp = ax.contourf(Δx:Δx:1, Δx:Δx:1,absmag; levels = 100, set_cmap = winter() )
    fig.colorbar(cp)
    #fig
    return fig
end 

function plot_pde_mass(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #fig, ax = PyPlot.subplots(figsize =(10, 10))
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := 2π *fa[a,b,θ]*eθ[θ,d]/Nθ
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    #absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
        #ax.yaxis.set_ticks([])
        ax.axis([Δx, 1., Δx, 1.])
        ax.set_aspect("equal")
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    cp = ax.contourf(Δx:Δx:1, Δx:Δx:1,ρ; levels = 100, set_cmap = winter() )
    fig.colorbar(cp)
    #fig
    return fig
end 

function plot_error(fig::Figure, ax::PyObject,param, t_saves, fa_saves, fp_saves)
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
    ax.plot(t_saves,dist_saves)
    ax.set_title("ρₐ = $(ρa), Pe = $(λ)")
    ax.axis([0, 0.3, 0., 0.15])
    ax.set_aspect("equal")
end 
##
#phase seperation metircs

function dist_from_unif(param, fa, fp)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    return 2*π*sqrt(sum( (fa .- ρa/(2*π) ).^2/(Nx*Nx*Nθ))) + sqrt(sum( (fp .- ρp).^2/(Nx*Nx)))
end

function time_dist_from_unif(param, fa_saves, fp_saves)
    n = length(fa_saves)
    dist_saves = zeros(n)
    for i in 1:n
        dist_saves[i] = dist_from_unif(param, fa_saves[i], fp_saves[i])
    end
    return dist_saves
end


#Example 
#=
#Parameters
param = pde_param(λ = 30., T = 0.003, name = "test7", Nθ =50, Dθ = 10., δt = 1e-5, ρa = 0.8)
density = initialize_density(param)
@unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S = param
@unpack fa, fp, t = density
perturb_pde!(param,density);

ρ = fp + 2*π*sum(fa;dims = 3)[:,:,1]/Nθ
p.(ρ)

midpoint_bond_diff(fp; Nx = Nx)
midpoint_bond_diff_θ(fa;  Nx = Nx,  Nθ = Nθ)
midpoint_bond_av(coeff_mag_s(fa,ρ); Nx =Nx )

@time for i in 1:10 pde_step!(param,density) end
fig, ax = PyPlot.subplots(figsize =(10, 10))
#plot_pde_mass(fig,ax,param,density)
plot_pde_mag(fig,ax,param,density)
display(fig)

@unpack fa, fp, t = density
dist_from_unif(param,fa,fp)

#expand variables
@unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S = param
@unpack fa, fp, t = density
#Run and save
run_pde_until!(param,density,0.002; save_on = true, max_steps = 20, save_interval = 0.001)
#for pmap
param = pde_param(λ = 0, T = 0.003, name = "test8")
perturb_pde_run(param; max_steps = 60, save_interval = 0.0005)
#Loading
t_saves, fa_saves, fp_saves = load_pdes(param,0.003; save_interval = 0.0005)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
#plotting
fig, ax = PyPlot.subplots(figsize =(10, 10))
y = self_diff.(x)
ax.plot(y)
display(fig)
density = initialize_density(param)
n = 1#length(t_saves)
t, fa, fp = t_saves[n], fa_saves[n], fp_saves[n]
@pack! density = fa, fp, t
#
pde_step!(param,density)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_pde_mass(fig,ax,param,density)
display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_pde_mag(fig,ax,param,density)
display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(t_saves,dist_saves)
display(fig)
##
pde_step!(param,density)
@unpack fa, fp, t = density
dist_from_unif(param,fa,fp)


for x₂ in 1:Nx
        local dsx, ρx, magx, y₁,y₂
        ## 2 direction

        y₂ = (x₂ +Nx)%Nx +1

        println(y₂)
end
=#
