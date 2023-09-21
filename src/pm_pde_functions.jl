cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
#println("Loading ...")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random, TensorOperations, StaticArrays


function pde_param_pm(; name = "test", D =1., Dx = 1., Pe =1., Dθ = 10, ρ= 0.5, χ = 1.0, Nx = 100, Nθ = 20, δt = 1e-5, T= 0.001, save_interval = 0.01, max_steps = 1e8, max_runs = 6, λ_step = 10., λmax = 100., λs = 20.:20.:100., pert = "n=1", δ = 0.01, k=20,γ = 0.0, video_length = 10000., cbar_max = 1.0, cbar_min = 0.0, frames = 1000, pert_interval = 5.)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    ρp = 0.0
    λ = Pe*sqrt(Dθ)
    ρp = (1-χ)*ρ
    ρa = χ*ρ
    param = Dict{String,Any}()
    @pack! param = k, name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E, Dθ, T, save_interval, max_steps, max_runs, λ_step, λmax, λs, pert, δ, Pe, Dx, χ, ρ, γ, video_length, cbar_max, cbar_min, frames, pert_interval
    return param
end

function initialize_density_pm(param::Dict{String,Any})
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    density = Dict{String,Any}()
    fa = fill(ρa/(2),(Nx,Nθ))
    fp = fill(ρp,(Nx))
    t = 0.
    @pack! density = fa , fp, t
    return density
end

## 

function midpoint_bond_diff_pm(f::Array{Float64,1}; Nx::Int64 = 100) 

    grad_f::Array{Float64,1} = zeros(Nx)

    for x₁ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        grad_f[x₁] = Nx*( f[y₁] - f[x₁] )
    end
    return grad_f
end

function midpoint_bond_diff_θ_pm(f::Array{Float64,2}; Nx::Int64 = 100,  Nθ::Int64 = 2) 

    grad_f::Array{Float64,2} = zeros(Nx,Nθ)

    for x₁ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        grad_f[x₁,:] = Nx*( f[y₁,:] - f[x₁,:] )
    end
    return grad_f
end

function midpoint_bond_av_pm(f::Array{Float64,1}; Nx::Int64 = 100) 
    av_f::Array{Float64,1}= zeros(Nx)

    for x₁ in 1:Nx
        ## 1 direction
        y₁::Int64 = (x₁ +Nx)%Nx +1
        av_f[x₁] = ( f[y₁] + f[x₁] )/2
    end
    return av_f
end

function midpoint_Θ_diff_pm(f::Array{Float64,2}; Nx::Int64 = 100,  Nθ::Int64 = 2) 
    grad_f::Array{Float64,2} = zeros(Nx, Nθ)

    for θ in 1:Nθ
        ϕ = ((θ % Nθ) +1)
        grad_f[:,θ] += Nθ*( f[:,ϕ] - f[:,θ] )/(2*π)
    end
    return grad_f
end

function site_div_pm(f::Array{Float64,1}; Nx::Int64 = 100) 

    div_f::Array{Float64,1} = zeros(Nx)

    for x₁ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        div_f[y₁] += Nx*( f[y₁] - f[x₁] )
    end
    return div_f
end

function site_div_θ_pm(f::Array{Float64,2}; Nx::Int64 = 100,  Nθ::Int64 = 2) 

    div_f::Array{Float64,2} = zeros(Nx, Nθ)

    for x₁ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        div_f[y₁,:] += Nx*( f[y₁,:] - f[x₁,:] )
    end

    return div_f
end

function θ_diff_pm(f::Array{Float64,2}; Nx::Int64 = 100,  Nθ::Int64 = 2) 

    div_f::Array{Float64,2} = zeros(Nx, Nθ)

    for θ in 1:Nθ
        ϕ ::Int64 = ((θ % Nθ) +1)
        div_f[:,θ] -= (f[:,ϕ]-f[:,θ])
    end

    return div_f
end

##


function self_diff(ρ::Float64;logtol = 1e-10, γ = 0.0)
    α::Float64= π/2 -1;
    if ρ ≤  0.
        ρ = 0.
    elseif ρ>1.
        ρ = 1.
    end
    return ( 1-ρ).*( α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ +1)+γ
end

function self_diff_prime(ρ::Float64)
    α::Float64= π/2 -1;
    return - ( α*(2*α-1)/(2*α+1)*ρ.^2 - α*ρ .+1) + ( -ρ .+1)*(2*α*(2*α-1)/(2*α+1)*ρ - α );
end

function mag_pm(f::Array{Float64,2}; Nθ = 2, Nx =100)
    return f[:,2] - f[:,1]
end

function coeff_mag_s_pm(f::Array{Float64,2},ρ::Array{Float64,1}; Nθ::Int64 = 2,  Nx::Int64 = 100,γ::Float64 = 0.0)
    m    ::Array{Float64,1} = mag_pm(f; Nθ=Nθ, Nx=Nx );
    ds   ::Array{Float64,1} = self_diff.(ρ; γ=γ);
    s    ::Array{Float64,1} = coeff_s.(ρ,ds);
    mag_s::Array{Float64,1} = s.*m
    return mag_s
end

function coeff_s(rho::Float64,ds::Float64)
    return ((rho*ds)>0 ? (1-rho-ds)/(rho*ds) : 0)
end

#functon p is labelled W in the pdf

function p(x::Float64;logtol = 1e-10, γ =0.)
    if γ ==0.
        if x <0.
            x = logtol
        elseif x ≥ 1.
            x = 1. -logtol
        end
        c2::Float64 = sqrt( (π-2)/( (π-1)*(26+(-11+π)*π))  )
        c3::Float64 = -(-4+π)sqrt( 1+8*(π-3)/( 26+(-11+π)*π)  )/2
        p1::Float64 = (1-π+2*(-3+π)*x)
        p2::Float64 = -2  +2*π  -(-2+π)*(-1+π)*x  +(-3+π)*(-2+π)*x^2
        return c3*log(-1-c2*p1)  -c3*log(1-c2*p1)  +(1 -π)*log(1-x)   + 0.5*(-2+π)*log(p2)
    else
        # if x <0.
        #     x = logtol
        # elseif x ≥ 1.
        #     x = 1. -logtol
        # end
        # a::Float64 = π/2 -1
        # coeff =[-1-2*a-γ-2*a*γ, 1+3*a+2*a^2,-4*a^2 ,-a+2*a^2];
        # rts::Vector{ComplexF64}= roots(Polynomial(coeff))
        # w = -(1/(1 + γ))*γ*log(x)
        # for r in rts
        #     denom = (1 + 3*a + 2 *a^2)+ (- 8* a^2)*r +(- 3 *a  + 6 *a^2 )*r^2
        #     neum = (1 + 3*a + 2a^2 )+(-4*a^2 )*r +(-a + 2*a^2)*r^2
        #     w += -(1/(1 + γ))*log(complex(x-r))*neum/denom
        # end
        # return real(w)
    end
end

function mob_pm(fa::Array{Float64,2}, fp::Array{Float64,1}, ρ::Array{Float64,1}; γ::Float64 =0.)
    ds::Array{Float64,1} = self_diff.(ρ; γ = γ )
    return fa.*ds, fp.*ds
end

function upwind(U::Float64, mb_down::Float64, mb_up::Float64)
    return (U   > 0. ? U.*mb_down  : U.*mb_up)
end

##

function U_velocities_pm(fa::Array{Float64,2}, fp::Array{Float64,1}, ρ::Array{Float64,1}; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., γ::Float64=0.)
    logtol::Float64 = log(1e-10);

    eθ:: Array{Float64,2} = reshape([-1 1],1,Nθ) 

    logmfa::Array{Float64,2} = map(x -> (x>0 ? log(x) : logtol), fa);
    logmfp::Array{Float64,1} = map(x -> (x>0 ? log(x) : logtol), fp);
    p_rho ::Array{Float64,1} = p.(ρ;γ=γ) #functon p is labelled W in the pdf

    Ua::Array{Float64,2}  = -midpoint_bond_diff_θ_pm(logmfa .+ p_rho; Nx=Nx, Nθ=Nθ).+ λ*midpoint_bond_av_pm(coeff_mag_s_pm(fa,ρ; Nθ=Nθ, Nx=Nx,γ=γ ); Nx =Nx ) .+ λ*eθ 
    Up::Array{Float64,1}  = -midpoint_bond_diff_pm(  logmfp  + p_rho; Nx=Nx       ) + λ*midpoint_bond_av_pm(coeff_mag_s_pm(fa,ρ; Nθ=Nθ, Nx=Nx,γ=γ ); Nx =Nx )
    #Uθ::Array{Float64,2}  = -midpoint_Θ_diff_pm(fa; Nx=Nx, Nθ = Nθ)

    return Ua, Up #, Uθ
end


function F_fluxes_pm(Ua::Array{Float64,2}, Up::Array{Float64,1}, moba::Array{Float64,2}, mobp::Array{Float64,1}; Nx::Int64 =100, Nθ::Int64 =100 )
    Fa ::Array{Float64,2} = zeros(Nx,Nθ);
    Fp ::Array{Float64,1} = zeros(Nx);
    # Fθ ::Array{Float64,2} = zeros(Nx,Nθ);
    for x₁ in 1:Nx
        local y₁
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        Fa[x₁,:] = upwind.(Ua[x₁,:], moba[x₁,:], moba[y₁,:])
        Fp[x₁]   = upwind( Up[x₁]  , mobp[x₁],   mobp[y₁]  )
    end
    # for θ in 1:Nθ
    #     local ϕ 
    #     ϕ::Int64 = ((θ % Nθ) +1)
    #     Fθ[:,θ] = upwind.(Uθ[:,θ], mobθ[:,θ], mobθ[:,ϕ])
    # end
    return Fa, Fp #, Fθ
end

##

function time_stepper_pm(fa::Array{Float64,2}, fp::Array{Float64,1}, δt::Float64; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., Dθ::Float64 = 10.,γ::Float64=0.)
    ρ::Array{Float64,1} = fp + sum(fa; dims =2)[:,1]
    
    Ua::Array{Float64,2},   Up::Array{Float64,1}   = U_velocities_pm(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ,γ=γ)
    moba::Array{Float64,2}, mobp::Array{Float64,1} = mob_pm(fa,fp,ρ;γ=γ)
    Fa::Array{Float64,2},   Fp::Array{Float64,1}  = F_fluxes_pm(Ua, Up, moba, mobp; Nx=Nx, Nθ=Nθ)
    
    a::Float64 = maximum(abs.(Ua));
    b::Float64 = maximum(abs.(Up));
    
    tempu::Float64 = 1/(6*max(a*Nx, b*Nx));
    dt::Float64= min(δt, tempu)

    fa -= dt*( site_div_θ_pm(Fa; Nx=Nx, Nθ=Nθ) + Dθ*θ_diff_pm(fa; Nx=Nx, Nθ=Nθ))
    fp -= dt*site_div_pm(Fp; Nx=Nx)

    return fa, fp, dt
end

function pde_stepper_pm!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ, Dθ, λ,γ = param
    @unpack fa, fp, t = density
    fa, fp , dt = time_stepper_pm(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ,γ=γ)
    t::Float64 += dt

    @pack! density = fa, fp, t
    return dt
end

##

function run_pde_until_pm!(param::Dict{String,Any},density::Dict{String,Any},T; save_on =false, max_steps = 100, save_interval = 1.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt, Dθ = param
    if save_on
        @unpack t = density
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
        safesave(filename,density)
    end
    if max_steps == false
        while density["t"] < T 
                time_since_save = 0.
                while time_since_save < min(save_interval, T)
                    dt = pde_stepper_pm!(param,density)
                    time_since_save += dt 
                end
                if save_on
                    @unpack t = density
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                    safesave(filename,density)
                end
            end
    else
        steps ::Int64 = 0
        while (density["t"] < T)&(steps < max_steps)
            time_since_save = 0.
            while (time_since_save < min(save_interval, T))&(steps < max_steps)
                dt = pde_stepper_pm!(param,density)
                time_since_save += dt 
                steps += 1
            end
            if save_on
                @unpack t = density
                if save_interval < 0.01
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                    safesave(filename,density)
                else
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                    safesave(filename,density)
                end
            end
        end
        return steps
    end
end

function perturb_pde_pm!(param::Dict{String,Any}, density::Dict{String,Any}; δ = 0.01, pert = "n=2")
    @unpack Nx, S, ρa, ρp, λ, Dθ, Nx, Nθ,Dx,Pe,Dθ,k, γ = param
    @unpack fa, fp = density
    ρ = ρa + ρp
    if ρ >0.99
        δ = min(δ, (1 - ρ)/(2*π+0.01));
    end
    #from stability: 
    if pert == "n=1"
        if ρp >0.
            K = collect(0:1:(k-1))
            matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k=k, γ= γ)
            ω = 2*π
            a, A = ap_MathieuEigen(matrix)
            Pa = (x,θ) -> real.( dot(A[2:1:(k+1),k+1],cos.(θ*K*(2*π/Nθ)))*exp(-im*x*ω/Nx) )
            Pp = (x) -> real.(A[1,k+1]*exp(-im*x*ω/Nx));
        else
            K = collect(0:1:(k-1))
            matrix = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k=k, γ= γ)
            ω = 2*π
            a, A = a_MathieuEigen(matrix)
            Pa = (x,θ) -> real.( dot(A[:,k],cos.(θ*K*(2*π/Nθ)))*exp(-im*x*ω/Nx) )
            Pp = (x) -> 0.;
        end
    end
    if pert == "pm_lin"
        ω, value, vector = pm_lin_pert(param)
        Prho = (x) -> real.( vector[1]*exp(im*x*2*π/Nx) );
        Pact = (x) -> real.( vector[2]*exp(im*x*2*π/Nx) );
        Pmag = (x) -> real.( vector[3]*exp(im*x*2*π/Nx) );
    end
    if pert == "rand"
        Pa = (x,θ) -> δ*ρa*(( rand() - 0.5 )/(ρa+0.01))/(2*π);
        Pp = (x) -> δ*ρp*( rand() - 0.5 )/(ρp+0.01);
    end
    #
    if pert == "safe_unif"
        Pa = (x,θ) -> ( rand() - 0.5 )/(2*π);
        Pp = (x) -> ( rand() - 0.5 );
    end
    #
    perta = zeros(Nx,Nθ)
    pertp = zeros(Nx)
    if pert == "safe_unif"
        for x₁ in 1:Nx, θ in S
            ρ = fp[x₁] + sum(fa[x₁,:])*2*π/Nθ
            perta[x₁, θ] += min(10*δ, ρ*(1-ρ)/2)*Pa(x₁, θ)/2;
        end
        for x₁ in 1:Nx
            ρ = fp[x₁] + sum(fa[x₁,:])*2*π/Nθ
            pertp[x₁] += min(10*δ, ρ*(1-ρ)/2)*Pp(x₁)/2;
        end
    elseif pert == "pm_lin"
        for x₁ in 1:Nx
            perta[x₁, 1] += (Pact(x₁)-Pmag(x₁))/2
            perta[x₁, 2] += (Pact(x₁)+Pmag(x₁))/2
            pertp[x₁]    += Prho(x₁)-Pact(x₁)
        end
    else
        for x₁ in 1:Nx, θ in S
            perta[x₁, θ] += Pa(x₁, θ);
        end
        for x₁ in 1:Nx
            pertp[x₁] += Pp(x₁);
        end
        if pert == "rand"
            perta[:, 1:(Nθ-1)]  = 0.5*perta[:, 1:(Nθ-1)] + 0.5*perta[:, (Nθ-1):(-1):1] 
        end
    end

    if pert == "safe_unif"
        fa += perta
        fp += pertp
    else
        c = dist_from_unif_pm(param, perta.+ρa/(2*π), pertp.+ρp)
        fa += δ*perta/c
        fp += δ*pertp/c
    end

    @pack! density = fa, fp;
end

function pm_lin_pert(param)
    @unpack S, ρa, ρp, λ, Nx, Nθ, Dx, Pe, Dθ, γ = param
    ω = 2*π/sqrt(Dθ);
    ϕa = ρa;
    ϕp = ρp;
    ϕ  = ϕa + ϕp;
    ϕ0 = 1- ϕ;
    ds = self_diff(ϕ);
    dsp = self_diff_prime(ϕ);
    DD = (1-ds)/ϕ
    s = DD - 1
    W = [-ω^2             0          -im*ω*Pe*ϕ0; 
        -ω^2*ϕa*DD      -ω^2*ds     -im*ω*Pe*(ϕa*s+ds); 
        -im*ω*Pe*ϕa*dsp -im*ω*Pe*ds -ω^2*ds-2         ]
    values,vectors = eigen(W)
    return ω, values[3], vectors[:,3]
end

function nudge_pde_pm!(param::Dict{String,Any}, density::Dict{String,Any}; δ = 0.01)
    @unpack Nx, S, ρa, ρp, λ, Dθ, Nx, Nθ,Dx,Pe,Dθ,k, γ = param
    @unpack fa, fp = density

    ρ = fp + sum(fa; dims = 2)[:,1]*2*π/Nθ
    x0 = argmax(ρ)

    K = collect(0:1:(k-1))
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k=k, γ= γ)
    ω = 2*π
    a, A = ap_MathieuEigen(matrix)
    Pa = (x,θ) -> real.( dot(A[2:1:(k+1),k+1],cos.(θ*K*(2*π/Nθ)))*exp(-im*x*ω/Nx) )
    Pp = (x) -> real.(A[1,k+1]*exp(im*x*ω/Nx));
    #A[:,k+1]
    #
    perta = zeros(Nx,Nθ)
    pertp = zeros(Nx)

    for x₁ in 1:Nx, θ in S
        ρ = fp[x₁] + sum(fa[x₁,:])*2*π/Nθ
        perta[x₁, θ] += min(δ, ρ*(1-ρ)/2)*Pa(x₁, θ)/5;
    end
    for x₁ in 1:Nx
        ρ = fp[x₁] + sum(fa[x₁,:])*2*π/Nθ
        pertp[x₁] += min(δ, ρ*(1-ρ)/2)*Pp(x₁)/5;
    end

    pert_ρ = pertp + sum(perta; dims = 2)[:,1]*2*π/Nθ
    x1 = argmin(pert_ρ)

    for x₁ in 1:Nx
        x2 = x₁+(x1-x0)
        if  x2<1
            x = x2+Nx
        elseif x2 > Nx
            x = x2-Nx
        else
            x = x2
        end

        fp[x₁] += pertp[x]
        for θ in S
            fa[x₁, θ] += perta[x, θ]
        end
    end

    @pack! density = fa, fp;
end

function perturb_pde_run_pm(param)
    @unpack T, save_interval, max_steps, pert, δ = param
    density = initialize_density_pm(param)
    perturb_pde_pm!(param,density; pert = pert, δ = δ);
    run_pde_until_pm!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)
end

function binodal_pde_run_pm(param)
    @unpack T, save_interval, max_steps = param
    density = initialize_sol_pm_full(param)
    run_pde_until_pm!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)
end

function repeating_perturb_pde_run_pm(param)
    @unpack T, save_interval, max_steps, pert, δ, pert_interval = param
    density = initialize_density_pm(param)
    t = 0.
    perturb_pde_pm!(param,density; pert = pert, δ = δ);
    while t<T
        run_pde_until_pm!(param,density,t+pert_interval; save_on = true, max_steps = max_steps, save_interval = save_interval)
        nudge_pde_pm!(param,density; δ = δ);
        t += pert_interval
    end
end

function load_pde_run_pm(param)
    @unpack T, save_interval, max_steps, pert, δ = param
    density = initialize_density_pm(param)
    try
        t_saves, fa_saves, fp_saves = load_pdes_pm(param,T; save_interval = 10*save_interval, start_time = 10*save_interval)
        k = length(t_saves)
        t, fa, fp = t_saves[k], fa_saves[k], fp_saves[k]
        @pack! density = t, fa, fp
        println("loaded t = $(t)")
    catch
        println("load failed")
        perturb_pde_pm!(param,density; pert = pert, δ = δ);
    end
    run_pde_until_pm!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)
end

function load_pde_nudge_run_pm(param)
    @unpack T, save_interval, max_steps, pert, δ, pert_interval = param
    density = initialize_density_pm(param)
    t = 0.
    try
        t_saves, fa_saves, fp_saves = load_pdes_pm(param,T; save_interval = 10*save_interval, start_time = 10*save_interval)
        k = length(t_saves)
        t, fa, fp = t_saves[k], fa_saves[k], fp_saves[k]
        @pack! density = t, fa, fp
        println("loaded t = $(t)")
    catch
        println("load failed")
        perturb_pde_pm!(param,density; pert = pert, δ = δ);
    end
    while t<T
        run_pde_until_pm!(param,density,t+pert_interval; save_on = true, max_steps = max_steps, save_interval = save_interval)
        nudge_pde_pm!(param,density; δ = δ);
        t += pert_interval
    end
end


#=
for i in 1:100
pde_stepper_pm!(param,density)
end
fig, ax = PyPlot.subplots(figsize =(10, 10))
@unpack fa, fp, t = density
ax.plot(fa)
display(fig)
t
=#
##
#phase seperation metircs

function dist_from_unif_pm(param, fa, fp)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param

    #=
    grad_fa = midpoint_bond_diff_θ_pm(fa; Nx = Nx,  Nθ = Nθ) 
    grad_fp =  midpoint_bond_diff_pm(fp; Nx = Nx)
    grad_grad_fa = midpoint_bond_diff_θ_pm(grad_fa; Nx = Nx,  Nθ = Nθ) 
    grad_grad_fp =  midpoint_bond_diff_pm(grad_fp; Nx = Nx) 

    partialθ_fa = midpoint_Θ_diff_pm(fa; Nx = Nx,  Nθ = Nθ)
    mixD_fa = midpoint_Θ_diff_pm(grad_fa; Nx = Nx,  Nθ = Nθ)
    partial_partialθ_fa = midpoint_Θ_diff_pm(partialθ_fa; Nx = Nx,  Nθ = Nθ)
    =#
    L2 = 2*π*sum( (fa .- ρa/2 ).^2)/(Nx*Nθ) + sum( (fp .- ρp).^2)/(Nx)
    #=
    H1x = 2*π*sum( (grad_fa ).^2)/(Nx*Nθ) + sum( (grad_fp).^2)/(Nx)

    H1θ = 2*π*sum( (partialθ_fa ).^2)/(Nx*Nθ)

    H2x = 2*π*sum( (grad_grad_fa ).^2)/(Nx*Nθ) + sum( (grad_grad_fp).^2)/(Nx)

    H2θ = 2*π*sum( (partial_partialθ_fa ).^2)/(Nx*Nθ)

    H2_mix = 2*π*sum( (mixD_fa ).^2)/(Nx*Nθ)
    =#
    return sqrt(L2) 
end

function dist_between_pm(param, fa1, fp1, fa2, fp2)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param

    L2 = 2*π*sum( (fa1- fa2 ).^2)/(Nx*Nθ) + sum( (fp1 - fp2).^2)/(Nx)
  
    return sqrt(L2) 
end

function time_dist_from_past_pm(param, fa_saves, fp_saves, skip)
    N = length(fa_saves)
    dist_saves = zeros(N-skip)
    for n in 1:(N-skip)
        fa1 = fa_saves[n]
        fa2 = fa_saves[n+skip]
        fp1 = fp_saves[n]
        fp2 = fp_saves[n+skip]
        dist_saves[n] = dist_between_pm(param, fa1, fp1, fa2, fp2)
    end
    return dist_saves
end

function time_dist_from_unif_pm(param, fa_saves, fp_saves)
    n = length(fa_saves)
    dist_saves = zeros(n)
    for i in 1:n
        dist_saves[i] = dist_from_unif_pm(param, fa_saves[i], fp_saves[i])
    end
    return dist_saves
end

# # loading
# using Pkg
# Pkg.add("JLD2")
# using JLD2


function load_pdes_pm(param::Dict{String,Any},T; save_interval = 1., start_time = 0.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt, Dθ = param
    t = start_time
    fa_saves = []
    fp_saves = []
    t_saves = []
    while t ≤ T
        try
            if save_interval < 0.01
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
            else
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
            end
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

# stability

function λsym(ϕ::Float64; Dθ::Float64 = 10., Dx = 1.)
    α = π/2 -1
    expr1 = Dθ*(1+2*α) + 4*Dx*π^2*(1-ϕ)*(1- 2*α^2*(1-ϕ)*ϕ - ϕ^2 + α*(2-ϕ+ϕ^2)  )
    expr2 = ϕ*(1 + 3*ϕ - 4*ϕ ^2 + 4*α^2*(1 - 3 *ϕ  + 2 *ϕ^2) + α* (4 - 6 *ϕ  + 4 *ϕ ^2)  )
    return  sqrt(8*Dx)*sqrt(1+2*α)*sqrt(expr1)/expr2
end

function refresh_stab_data_pm(;stabdata = Dict{String,Any}(), ρs = 0.05:0.05:0.95,   Dθ = 10., Nx = 50, Nθ = 20, λs = 5.:5.:100., name = "high_density_stability_v4", save_on = true, t_end = 1.0, ρp = 0.)
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_ρp=$(ρp).jld2"
    if save_on
        try 
            stabdata = wload(filename)
        catch
        end
    end

    # for ρ in 0.90:0.01:0.99
    for ρ in ρs
        stable = []
        unstable = []
        unsure = []
        lin_stab = []
        data = Dict{String,Any}()
        #load ans
        for λ ∈ λs
            try
                param = pde_param_pm(;name = name, Pe = λ/sqrt(Dθ) , ρa = ρ-ρp, ρp = ρp, T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, save_interval = 0.01, max_steps = 1e8)
                t_saves, fa_saves, fp_saves = load_pdes_pm(param,0.02; save_interval = 0.001, start_time = 0.0)

                stab_dsit0 = dist_from_unif_pm(param, fa_saves[1], fp_saves[1])

                stab_dsit1 = dist_from_unif_pm(param, fa_saves[2], fp_saves[2])

                #=
                if (stab_dsit0>stab_dsit1)&(λ ∉ lin_stab)
                    push!(lin_stab, λ)
                end
                =#
            catch
                    #break
            end
            try
                    param = pde_param_pm(;name = name, Pe = λ/sqrt(Dθ) , ρa = ρ-ρp, ρp = ρp, T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, save_interval = 0.01, max_steps = 1e8)
                    t_saves, fa_saves, fp_saves = load_pdes_pm(param,t_end; save_interval = 0.001, start_time = (t_end-0.02))

                    stab_dsit = dist_from_unif_pm(param, fa_saves[1], fp_saves[1])

                    if (stab_dsit>0.02)&(λ ∉ unstable)
                        push!(unstable,  λ)
                    elseif (stab_dsit<0.01)&(λ ∉ stable)&(λ ∉ lin_stab)
                        push!(stable, λ)
                    elseif (λ ∉ unsure)
                        push!(unsure, λ)
                    end
            catch
                        #break
            end
        end
        #calc ans
        λcrit = 0 #λsym(ρ; Dθ = Dθ)
        #=
        for λ ∈ λs
                if (λ ≤ λcrit)&(λ ∉ lin_stab)
                    push!(lin_stab,  λ)
                end
        end
        =#
        #fillout ans
            try
                    λmin = maximum(stable)
                    for λ ∈ λs
                        if (λ ≤ λmin)&(λ ∉ stable)&(λ ∉ lin_stab)
                            push!(stable,  λ)
                        end
                    end
            catch
            end
            #=
            try 
                λminlin = maximum(lin_stab)
                for λ ∈ λs
                    if (λ ≤ λminlin)&(λ ∉ lin_stab)
                        push!(lin_stab,  λ)
                    end
                end
            
            catch
            end
            =#
            try 
                λmax = minimum(unstable)
                for λ ∈ λs
                    if (λ ≥ λmax)&(λ ∉ unstable)
                        push!(unstable,  λ)
                    end
                end
            catch
            end
            @pack! data = λcrit, stable, unstable, unsure, lin_stab
            stabdata["ρ = $(ρ)"] = data
    end
    if save_on
        wsave(filename,stabdata)
    end
    #println(stabdata["ρ = $(ρs[1])"])
    return stabdata
end

function next_param_pm(stabdata, param; λmax = 1e8, λ_step = 10.)
    @unpack ρa, Dθ, = param
    stable = []
    unstable = [] 
    unsure = []
    try
        @unpack λcrit, stable, unstable, unsure = stabdata["ρ = $(ρa)"]
    catch
    end
    λmin = λsym(ρa; Dθ=Dθ, Dx=1. )
    try 
        λmin = maximum(stable)
    catch
    end
    λmax = λmax
    try
        λmax = minimum(unstable)
    catch
    end
    λmid = λmax
    try
        λmid = minimum(unsure)
    catch
    end

    if (λmid+λ_step < λmax)
        next_λ = λmid
        unfinished = true
    elseif (λmin+λ_step < λmid)
        next_λ = λmin
        unfinished = true
    else
        next_λ = 0.
        unfinished = false
    end

    approx_next = 0.
    while next_λ ≥ approx_next
        approx_next += λ_step
    end

    λ = approx_next
    @pack! param = λ
    return unfinished, approx_next, param
end

function run_stab_search_pm(param)
    @unpack name, Nx, Nθ, max_runs, λ_step, ρa, Dθ, Nx, Nθ, λmax, λs = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ).jld2"
    stabdata = Dict{String,Any}()
    try 
        stabdata = wload(filename)
    catch
        println(" error for ρ = $(ρa)")
        host = gethostname()
        println(host)
    end

    i = 0
    unfinished = true
    while  (i < max_runs)&unfinished
        unfinished, λ_next, param = next_param_pm(stabdata, param; λ_step = λ_step, λmax=λmax)
        if unfinished
                println("running λ = $(λ_next) ρ = $(ρa)")
                perturb_pde_run_pm(param)
                stabdata = refresh_stab_data_pm(;stabdata = stabdata,  ρs = [ρa],   Dθ = Dθ, Nx = Nx, Nθ = Nθ, λs = λs , name = name, save_on = false)
                println("finished λ = $(λ_next) ρ = $(ρa)")
        end
        i += 1 
        println("iteration = $(i), unfinished = $(unfinished)")
        println("while condition = $((i < max_runs)&unfinished))")
    end
    println("i am done")

    return ρa, stabdata["ρ = $(ρa)"]
end

function run_stab_search_stupid_mode_pm(param)
    @unpack name, Nx, Nθ, max_runs, λ_step, ρa, Dθ, Nx, Nθ, λmax, λs = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ).jld2"
    stabdata = Dict{String,Any}()
    try 
        stabdata = wload(filename)
    catch
        println(" no save for ρ = $(ρa)")
        host = gethostname()
        println(host)
    end

    i = 0
    unfinished, λ_next, param = next_param_pm(stabdata, param; λ_step = λ_step, λmax=λmax)
    while  (i < max_runs)&unfinished
        @unpack λ = param
        println("running λ = $(λ) ρ = $(ρa)")
        perturb_pde_run_pm(param)
        #stabdata = refresh_stab_data(;stabdata = stabdata,  ρs = [ρa],   Dθ = Dθ, Nx = Nx, Nθ = Nθ, λs = λs , name = name, save_on = false)
        println("finished λ = $(λ) ρ = $(ρa)")
        λ += λ_step
        @pack! param = λ
        i += 1 
        println("iteration = $(i)")
    end
    println("i am done")
end

#

function rho_eθ_av(f, eθ; Nx =Nx, Nθ=Nθ )
    av_f::Array{Float64,3}= zeros(2, Nx, Nx)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ = (x₁ +Nx)%Nx +1
        y₂ = x₂
        av_f[1,x₁,x₂] = ( f[y₁,y₂] + f[x₁,x₂] )/2
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        av_f[2,x₁,x₂] = ( f[y₁,y₂] + f[x₁,x₂] )/2
    end

    av_f2::Array{Float64,4}= zeros(2, Nx, Nx, Nθ)

    for x₁ in 1:Nx, x₂ in 1:Nx, θ in 1:Nθ
        av_f2[1,x₁,x₂,θ] = av_f[1,x₁,x₂]*eθ[1,1,1,θ]
        av_f2[2,x₁,x₂,θ] = av_f[2,x₁,x₂]*eθ[2,1,1,θ]
    end

    return av_f2
end

function U_velocities_sym(fa::Array{Float64,3}, fp::Array{Float64,2}, ρ::Array{Float64,2}; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., ρa = 0.6)
    logtol::Float64 = log(1e-10);

    eθ:: Array{Float64,4} = reshape([cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)]',2,1,1,Nθ)

    logmfa::Array{Float64,3} = map(x -> (x>0 ? log(x) : logtol), fa);
    logmfp::Array{Float64,2} = map(x -> (x>0 ? log(x) : logtol), fp);
    p_rho ::Array{Float64,2} = p.(ρ) #functon p is labelled W in the pdf

    m ::Array{Float64,3} = mag(fa; Nθ=Nθ, Nx=Nx );
    alpha_ds ::Float64 = ((1-ρa-self_diff(ρa))/ρa - self_diff_prime(ρa))/(2*self_diff(ρa))


    Ua::Array{Float64,4}  = -midpoint_bond_diff_θ(logmfa .+ p_rho; Nx=Nx, Nθ=Nθ) + λ*alpha_ds*( rho_eθ_av(ρ, eθ; Nx =Nx, Nθ=Nθ ) .+ - midpoint_bond_av(m; Nx =Nx )  )
    Up::Array{Float64,3}  = -midpoint_bond_diff(  logmfp  + p_rho; Nx=Nx       ) 
    Uθ::Array{Float64,3}  = -midpoint_Θ_diff(fa; Nx=Nx, Nθ = Nθ)

    return Ua, Up, Uθ
end


function F_fluxes_sym(Ua::Array{Float64,4}, Up::Array{Float64,3}, Uθ::Array{Float64,3}, moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3}; Nx::Int64 =100, Nθ::Int64 =100 )
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

function time_step_sym(fa::Array{Float64,3}, fp::Array{Float64,2}, δt::Float64; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., Dθ::Float64 = 10.)
    ρ::Array{Float64,2} = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    
    Ua::Array{Float64,4},   Up::Array{Float64,3},   Uθ::Array{Float64,3}   = U_velocities_sym(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ)
    moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3} = mob(fa,fp,ρ)
    Fa::Array{Float64,4},   Fp::Array{Float64,3},   Fθ::Array{Float64,3}   = F_fluxes_sym(Ua, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)
    
    a::Float64 = maximum(abs.(Ua));
    b::Float64 = maximum(abs.(Up));
    c::Float64 = maximum(abs.(Uθ));
    
    tempu::Float64 = 1/(6*max(a*Nx, b*Nx, c*Nθ*Dθ/(2*π)));
    dt::Float64= min(δt, tempu)

    fa -= dt*( site_div_θ(Fa; Nx=Nx, Nθ=Nθ) + Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ))
    fp -= dt*site_div(Fp; Nx=Nx)

    return fa, fp, dt
end

function pde_step_sym!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ, Dθ, λ = param
    @unpack fa, fp, t = density
    fa, fp , dt = time_step_sym(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)
    t::Float64 += dt

    @pack! density = fa, fp, t
    return dt
end

# plots 
using TensorOperations, PyPlot, PyCall

function make_phase_video_pm(param; frames = 100)
    @unpack T, save_interval = param
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes_pm(param,T; save_interval = save_interval)
    frames = Int64(round(length(t_saves)))-1
    animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves; frames = frames-1)
end


function vid_pde_plot_pm(fig::Figure, axs, param::Dict{String,Any}, t_saves, fa_saves, fp_saves, i)
    @unpack Nx, Nθ, ρa, ρp, χ, Dθ, Dx, k, γ,Pe = param
    ρa_saves, ρp_saves = deepcopy(spatial_density_pm.(fa_saves)), deepcopy(fp_saves)

    push!(ρa_saves[i], ρa_saves[i][1])
    push!(ρp_saves[i], ρp_saves[i][1])

    ρsum = ρp_saves[i]+ρa_saves[i]

    axs[1].plot((0:1:Nx)/Nx,ρa_saves[i].-ρa, color = "red", label = L"\rho^a - \phi^a")
    axs[1].plot((0:1:Nx)/Nx,ρsum.-ρa.-ρp, color = "black", label = L"\rho - \phi")
    axs[1].plot((0:1:Nx)/Nx,ρp_saves[i].-ρp, color = "blue", label = L"\rho^p - \phi^p")

    axs[1].xaxis.set_ticks(0.:0.2:1.0)
    axs[1].xaxis.set_tick_params(labelsize=15)
    axs[1].yaxis.set_tick_params(labelsize=15)
    rhomax = maximum(maximum(ρa_saves).-ρa)+maximum(maximum(ρp_saves).-ρp)
    axs[1].axis([0., 1., -rhomax , rhomax])
    #axs[1].axis([0., 1., min(minimum(minimum(ρa_saves)),minimum(minimum(ρp_saves))),maximum(maximum( ρa_saves+ρp_saves ))])
    axs[1].set_xlabel(L"x",fontsize=20)
    #axs[1].set_ylabel(L"\rho,",fontsize=20)
    title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ), \\phi = $(ρa+ρp), \\mathrm{Pe} = $(round(Pe; digits = 3)), t = $(round(t_saves[i]; digits = 3))\$")
    axs[1].set_title(title,fontsize=20)

    mat1 = zeros(1, Nx+1)
    mat2= zeros(1, Nx+1)
    mags = mag_pm(fa_saves[i]; Nθ = Nθ)
    push!(mags,mags[1])
    mat1[1, :] = mags
    mat2[1, :] = mags.*(-ρsum.+1)

    #colmap = PyPlot.plt.cm.seismic
    colmap = PyPlot.plt.cm.PRGn
    norm1 = matplotlib.colors.Normalize(vmin= -rhomax*0.5 , vmax= rhomax*0.5) 
    #norm1 = matplotlib.colors.Normalize(vmin= -maximum(abs.(mags)) , vmax= maximum(abs.(mags)) )
    #norm2 = matplotlib.colors.Normalize(vmin= minimum(mags/10) , vmax= maximum(mags)/10 )

    axs[2].matshow(mat1; norm = norm1,  cmap = colmap, extent = [0., 1., 0., 0.1])
    #axs[3].matshow(mat2; norm = norm2,  cmap = colmap, extent = [0., 1., 0., 0.1])

    axs[2].set_aspect(1.)
    #axs[3].set_aspect(1.)

    axs[2].xaxis.set_ticks(0.:0.2:1.0)
    axs[2].yaxis.set_ticks([])
    axs[2].xaxis.set_tick_params(labelsize=15)
    axs[2].xaxis.tick_bottom()
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    #ax.set_xlabel(L"x",fontsize=20)

    axs[2].set_ylabel(L"\mathbf{p}", fontsize=20, rotation=0)
    axs[2].yaxis.set_label_coords(-.05, .5)

    lines, labels = axs[1].get_legend_handles_labels()
    fig.tight_layout()
    ldg = fig.legend(lines, labels, loc = "center right", fontsize=20, bbox_to_anchor = (0.25, 0.25, 1, 1),
    bbox_transform = plt.gcf().transFigure)

    return fig
end

function animate_phase_pdes_pm(param,t_saves,fa_saves,fp_saves; frames = 99)
    @unpack name, λ, ρa, ρp, Nx, Nθ, δt, Dθ, χ, γ = param
    fig, axs = plt.subplots(2, 1, figsize=(10,10))
    function makeframe(i)
        clf()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        axs = ax1, ax2
        vid_pde_plot_pm(fig, axs, param, t_saves, fa_saves, fp_saves, i+1)
        return fig
    end
    interval = 5*Int64(round(20000/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    T = t_saves[Int64(round((frames+1)))]
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(T; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function cal_rho_saves(fa, fp; Nθ = Nθ)
    return fp + sum(fa; dims =2)[:,1]
end

function spatial_density_pm(fa)
    return sum(fa; dims =2)[:,1]
end

function bump_funciton(x; ϵ = 1, x0 = 0.)
    if (-1<((x-x0)/ϵ)<1)
        return exp(-1/(1-((x-x0)/ϵ)^2))/ϵ
    else
        return 0
    end
end

function initialize_bin_pm(param::Dict{String,Any}; )
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    ϕa = ρa
    ϕp = ρp
    initial_Δ = 1e-4;
    max_iter = 40;
    tol = 1e-3;
    atol = 1e-12;
    rho_max = (1-10e-20);
    γ = (1-ϕa)/(1-ρ);

    find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ, max_iter = max_iter, tol = tol, atol = atol);
    ϕg = lower_limits[1]
    ϕl = upper_limits[1]

    ϕag, ϕpg = gamma_converter(γ, ϕg)
    ϕal, ϕpl = gamma_converter(γ, ϕl)

    density = Dict{String,Any}()
    fa = fill(ϕag/(2),(Nx,Nθ))
    fp = fill(ϕpg,(Nx))
    Nx2 = Int64(round(Nx/2))
    fa[1:Nx2,:] = fill(ϕal/(2),(Nx2,Nθ))
    fp[1:Nx2] = fill(ϕpl,(Nx2))
    t = 0.
    @pack! density = fa , fp, t
    return density
end

function smooth_density(param,density)
    @unpack name, Dx, Dθ, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    @unpack t, fa, fp = density 
    
    mollifier = bump_funciton.((1:Nx)/Nx; ϵ = 5*sqrt(Dx/Dθ), x0 = 1/2)
    mollifier = mollifier/sum(mollifier)

    Nx2 = Int64(round(Nx/2))

    fa2 = zeros(Nx,2)
    fp2 = zeros(Nx)
    for i in 1:Nx, j in 1:Nx
        ## k = i - j 
        k = (i-j +Nx2 +Nx -1 )%Nx +1
        fa2[i,1] += mollifier[k].*fa[j,1]
        fa2[i,2] += mollifier[k].*fa[j,2]
        fp2[i]   += mollifier[k].*fp[j]
    end

    density = Dict{String,Any}()
    fa = fa2
    fp = fp2
    t = 0.
    @pack! density = fa , fp, t
    return density
end


function initialize_arctan_pm(param::Dict{String,Any}; )
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    ϕa = ρa
    ϕp = ρp
    ϕ = ϕa+ϕp
    initial_Δ = 1e-4;
    max_iter = 40;
    tol = 1e-3;
    atol = 1e-12;
    rho_max = (1-10e-20);
    γ = (1-ϕa)/(1-ϕ);

    find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ, max_iter = max_iter, tol = tol, atol = atol);
    ϕg = lower_limits[1]
    ϕl = upper_limits[1]

    ϕag, ϕpg = gamma_converter(γ, ϕg)
    ϕal, ϕpl = gamma_converter(γ, ϕl)


    l = sqrt(Dx/Dθ)
    ϵ = 2.0
    
    density = Dict{String,Any}()

    ρa  = (ϕal-ϕag)*(  tanh.( (collect(1:Nx)/Nx .-0.75)/ϵ/l )/2 .+0.5) .+ ϕag
    fp  = (ϕpl-ϕpg)*(tanh.( (collect(1:Nx)/Nx .-0.75)/ϵ/l )/2 .+0.5) .+ ϕpg
    ρ = ρa+fp

    m   = (ϕl-ϕg)*( -tanh.( (collect(1:Nx)/Nx .-0.75)/ϵ/l ).^2 .+1 )/2 ./(-ρ.+1)/λ/l/ϵ

    Nx2 = Int64(round(Nx/2))
    ρa[1:Nx2]  = -(ϕal-ϕag)*(  tanh.( (collect(1:Nx2)/Nx .-0.25)/ϵ/l )/2 .+0.5  ) .+ ϕal
    fp[1:Nx2] = -(ϕpl-ϕpg)*(tanh.( (collect(1:Nx2)/Nx .-0.25)/ϵ/l )/2 .+0.5  ) .+ ϕpl
    ρ = ρa+fp
    m[1:Nx2]  = -(ϕl-ϕg)*( -tanh.( (collect(1:Nx2)/Nx .-0.25)/ϵ/l ).^2 .+1 )/2 ./(-ρ[1:Nx2].+1)/λ/l/ϵ
    

    fa = zeros(Nx,2)
    fa[:,1] = (ρa - m)/2
    fa[:,2] = circshift((ρa + m)/2,0)
    t = 0.

    @pack! density = fa , fp, t
    return density
end


function initialize_sol_pm(param::Dict{String,Any}, sol, ϕg, ϕl, γ;)
   @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E= param
    γ = (1-ρa)/(1-ρa-ρp);
    #find central time
    t_mid_arg = argmax(sol[2,:])
    t_middle = sol.t[t_mid_arg]
    t_max = maximum(sol.t[:].-t_middle)
    t_min = minimum(sol.t[:].-t_middle)
    t_lim = round(min(t_max, - t_min))
    #t_lim_index = argmin(abs.(sol.t[:].-2*t_lim))
    
    l = sqrt(Dx/Dθ)
    density = Dict{String,Any}()
    
    ρ = fill(ϕg, Nx)
    m = fill(0., Nx)
    
    Nx4  = Int64(round(Nx/4))
    Nx2  = 2*Nx4
    Nx34 = 3*Nx4
    
    ρ[(Nx4+1):1:Nx34] = fill(ϕl, Nx2)
    
    n_index = min(Int64(round(t_lim*l*Nx)),Nx4)
    
    ts1 = collect((1:(2*n_index))/Nx/l) .+ sol.t[t_mid_arg] .-n_index/Nx/l
    ts2 = collect(((2*n_index):(-1):1)/Nx/l) .+ sol.t[t_mid_arg] .-n_index/Nx/l
    
    sol_index_1 = (Nx4-n_index+1):1:(Nx4+n_index)
    sol_index_2 = (Nx34-n_index+1):1:(Nx34+n_index)
    
    ρ[sol_index_1] = sol.(ts1, idxs = 1)
    m[sol_index_1] = sol.(ts1, idxs = 2)
    
    ρ[sol_index_2] =  sol.(ts2, idxs = 1)
    m[sol_index_2] = -sol.(ts2, idxs = 2)
    
    ρa = γ*(ρ .- 1) .+1
    fp = ρ - ρa
    fa = zeros(Nx,2)
    fa[:,1] = (ρa - m)/2
    fa[:,2] = (ρa + m)/2
    t = 0.
    
    @pack! density = fa , fp, t
    return density
end

function spatial_currents_pm(fa_saves, fp_saves; param = param)
    @unpack Nx, Nθ, ρa, ρp, χ, Dθ, Dx, k, γ,Pe,λ  = param
    ja_saves, jp_saves = [], []
    for i in eachindex(fa_saves)
        fa = fa_saves[i]
        fp =  fp_saves[i]

        ρ::Array{Float64,1} = fp + sum(fa; dims =2)[:,1]
    
        Ua::Array{Float64,2},   Up::Array{Float64,1} = U_velocities_pm(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ,γ=γ)
        moba::Array{Float64,2}, mobp::Array{Float64,1} = mob_pm(fa,fp,ρ;γ=γ)
        Fa::Array{Float64,2},   Fp::Array{Float64,1}  = F_fluxes_pm(Ua, Up, moba, mobp; Nx=Nx, Nθ=Nθ)
        
        ja = sum(Fa; dims =2)[:,1]
        jp = Fp
        push!(ja_saves,ja)
        push!(jp_saves,jp)
    end
    return ja_saves, jp_saves
end

using DifferentialEquations

function f(du,u,parameters,t)
    Pe = parameters[1]
    γ = parameters[2]
    ϕ1 = parameters[3]
    du[1] = Pe*(1-u[1])*u[2]
    du[2] = -Pe*u[2]^2 + Pe*( (1-γ*(1-u[1]))*self_diff(u[1]) -(1-γ*(1-ϕ1))*self_diff(ϕ1) )/self_diff(u[1]) -(2/Pe)*log( (1-u[1])/(1-ϕ1) )/self_diff(u[1])
    return du
end

function f_jac(J,u,parameters,t)
    Pe = parameters[1]
    γ = parameters[2]
    ϕ1 = parameters[3]
    J[1,1] = -Pe*u[2]
    J[1,2] =  Pe*(1-u[1])
    J[2,1] =  self_diff_prime(u[1])*            (2/Pe)*log( (1-u[1])/(1-ϕ1) )/self_diff(u[1])^2              + (2/Pe)/(1-u[1])/self_diff(u[1])
    J[2,1] += self_diff_prime(u[1])*Pe*((1-γ*(1-u[1]))*self_diff(u[1]) -(1-γ*(1-ϕ1))*self_diff(ϕ1) )/self_diff(u[1])^2  + Pe*(γ*self_diff(u[1]) +(1-γ*(1-u[1]))*self_diff_prime(u[1]) )/self_diff(u[1])
    J[2,2] = -2*Pe*u[2]
    return J
end

function initialize_sol_pm_full(param::Dict{String,Any})
    @unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, Pe, Dx, Dθ= param
    #copmute binodal values

    ϕa = ρa
    ϕp = ρp
    ϕ = ϕa+ϕp

    initial_Δ = 1e-4;
    max_iter = 40;
    tol = 1e-16;
    atol = 1e-16;
    rho_max = (1-10e-20);
    γ = (1-ϕa)/(1-ϕ);

    find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ, max_iter = max_iter, tol = tol, atol = atol);
    ϕg = lower_limits[1]
    ϕl = upper_limits[1]
    ϕg, ϕl 
    find_sol

    #compute solution
    J = zeros(2,2)
    u = zeros(2)
    parameters = (Pe, γ, ϕg, ϕl)
    t = 0.
    J = f_jac(J,u,parameters,t)
    values, vectors = eigen(J)
    evector2 = vectors[:,2]

    ff = ODEFunction(f;jac=f_jac)
    ϵ = 1e-15
    initial_position = [ϕg, 0.0] + ϵ*evector2
    time_interval = (0.0, 15.0)

    ff = ODEFunction(f;jac=f_jac)
    prob = ODEProblem(ff,initial_position,time_interval, parameters)

    sol = DifferentialEquations.solve(prob,abstol = 1e-14, reltol = 1e-14);

    #find central time
    t_mid_arg = argmax(sol[2,:])
    t_middle = sol.t[t_mid_arg]
    t_max = maximum(sol.t[:].-t_middle)
    t_min = minimum(sol.t[:].-t_middle)
    t_lim = round(min(t_max, - t_min))
    #t_lim_index = argmin(abs.(sol.t[:].-2*t_lim))
    
    l = sqrt(Dx/Dθ)
    density = Dict{String,Any}()
    
    ρ = fill(ϕg, Nx)
    m = fill(0., Nx)
    
    Nx4  = Int64(round(Nx/4))
    Nx2  = 2*Nx4
    Nx34 = 3*Nx4
    
    ρ[(Nx4+1):1:Nx34] = fill(ϕl, Nx2)
    
    n_index = min(Int64(round(t_lim*l*Nx)),Nx4)
    
    ts1 = collect((1:(2*n_index))/Nx/l) .+ sol.t[t_mid_arg] .-n_index/Nx/l
    ts2 = collect(((2*n_index):(-1):1)/Nx/l) .+ sol.t[t_mid_arg] .-n_index/Nx/l
    
    sol_index_1 = (Nx4-n_index+1):1:(Nx4+n_index)
    sol_index_2 = (Nx34-n_index+1):1:(Nx34+n_index)
    
    ρ[sol_index_1] = sol.(ts1, idxs = 1)
    m[sol_index_1] = sol.(ts1, idxs = 2)
    
    ρ[sol_index_2] =  sol.(ts2, idxs = 1)
    m[sol_index_2] = -sol.(ts2, idxs = 2)
    
    ρa = γ*(ρ .- 1) .+1
    fp = ρ - ρa
    fa = zeros(Nx,2)
    fa[:,1] = (ρa - m)/2
    fa[:,2] = (ρa + m)/2
    t = 0.
    
    @pack! density = fa , fp, t
    return density
end

println("loaded")
