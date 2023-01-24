cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")
#runnning simulation
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random, TensorOperations, StaticArrays

function pde_param_k(; name = "test", D =1. , Pe =1. ,ρa = 0.5, ρp = 0.0, Nx = 100, Nθ = 100, δt = 1e-5, Dθ = 10, T= 0.001, save_interval = 0.01, max_steps = 1e8, max_runs = 6, λ_step = 10., λmax = 100., λs = 20.:20.:100., pert = "n=1", δ = 0.01,Dx = 1.,k=20)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    λ = Pe*sqrt(Dθ)
    param = Dict{String,Any}()
    @pack! param = k, name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E, Dθ, T, save_interval, max_steps, max_runs, λ_step, λmax, λs, pert, δ, Pe, Dx
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
        grad_f[:,:,θ] += Nθ*( f[:,:,ϕ] - f[:,:,θ] )/(2*π)
    end
    return grad_f
end

function site_div(f::Array{Float64,3}; Nx::Int64 = 100) 

    div_f::Array{Float64,2} = zeros(Nx, Nx)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        div_f[y₁,y₂] += Nx*( f[1,y₁,y₂] - f[1,x₁,x₂] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        div_f[y₁,y₂] += Nx*( f[2,y₁,y₂] - f[2,x₁,x₂] ) 
    end
    return div_f
end

function site_div_θ(f::Array{Float64,4}; Nx::Int64 = 100,  Nθ::Int64 = 100) 

    div_f::Array{Float64,3} = zeros(Nx, Nx, Nθ)

    for x₁ in 1:Nx, x₂ in 1:Nx
        ## 1 direction
        y₁ ::Int64 = (x₁ +Nx)%Nx +1
        y₂ ::Int64 = x₂
        div_f[y₁,y₂,:] += Nx*( f[1,y₁,y₂,:] - f[1,x₁,x₂,:] )
        ## 2 direction
        y₁ = x₁
        y₂ = (x₂ +Nx)%Nx +1
        div_f[y₁,y₂,:] += Nx*( f[2,y₁,y₂,:] - f[2,x₁,x₂,:] ) 
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

##

function self_diff(ρ::Float64;logtol = 1e-10, γ = 0.0)
    α::Float64= π/2 -1;
    if ρ ≤  0.
        ρ = logtol
    elseif ρ>1.
        ρ = 1.
    end
    return ( 1-ρ).*( α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ +1)+γ
end

function self_diff_prime(ρ::Float64;logtol = 1e-10)
    α::Float64= π/2 -1;
    if ρ ≤  0.
        ρ = logtol
    elseif ρ>1.
        ρ = 1.
    end
    return - ( α*(2*α-1)/(2*α+1)*ρ.^2 - α*ρ .+1) + ( -ρ .+1)*(2*α*(2*α-1)/(2*α+1)*ρ - α );
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
    return ((rho*ds)>0 ? (1-rho-ds)/(rho*ds) : 0)
end

function coeff_mag_s(f::Array{Float64,3},ρ::Array{Float64,2}; Nθ::Int64 = 100,  Nx::Int64 = 100,γ::Float64 = 0.0)
    m    ::Array{Float64,3} = mag(f; Nθ=Nθ, Nx=Nx );
    ds   ::Array{Float64,2} = self_diff.(ρ; γ=γ);
    s    ::Array{Float64,3} = reshape(coeff_s.(ρ,ds),1,Nx,Nx);
    mag_s::Array{Float64,3} = s.*m
    return mag_s
end
#functon p is labelled W in the pdf
#=
function p(x::Float64; logtol = 1e-10)
    if x <0.
        x = logtol
    elseif x ≥ 1.
        x = 1. -logtol
    end
    c2::Float64 = sqrt( (π-2)/( (π-1)*(26+(-11+π)*π))  )
    c3::Float64 = -(-4+π)sqrt( 1+8*(π-3)/( 26+(-11+π)*π)  )/2
    p1::Float64 = (1-π+2*(-3+π)*x)
    p2::Float64 = -2  +2*π  -(-2+π)*(-1+π)*x  +(-3+π)*(-2+π)*x^2
    return c3*log(abs(-1-c2*p1))  -c3*log(abs(1-c2*p1))  +(1 -π)*log(1-x)   + 0.5*(-2+π)*log(abs(p2))
end
=#
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
        if x <0.
            x = logtol
        elseif x ≥ 1.
            x = 1. -logtol
        end
        a::Float64 = π/2 -1
        coeff =[-1-2*a-γ-2*a*γ, 1+3*a+2*a^2,-4*a^2 ,-a+2*a^2];
        rts::Vector{ComplexF64}= roots(Polynomial(coeff))
        w = -(1/(1 + γ))*γ*log(x)
        for r in rts
            denom = (1 + 3*a + 2 *a^2)+ (- 8* a^2)*r +(- 3 *a  + 6 *a^2 )*r^2
            neum = (1 + 3*a + 2a^2 )+(-4*a^2 )*r +(-a + 2*a^2)*r^2
            w += -(1/(1 + γ))*log(complex(x-r))*neum/denom
        end
        return real(w)
    end
end

function mob(fa::Array{Float64,3}, fp::Array{Float64,2}, ρ::Array{Float64,2}; γ=0.)
    ds::Array{Float64,2} = self_diff.(ρ;γ=γ)
    return fa.*ds, fp.*ds, fa
end

function upwind(U::Float64, mb_down::Float64, mb_up::Float64)
    return (U   > 0. ? U.*mb_down  : U.*mb_up)
end

##

function U_velocities(fa::Array{Float64,3}, fp::Array{Float64,2}, ρ::Array{Float64,2}; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10.)
    logtol::Float64 = log(1e-10);

    eθ:: Array{Float64,4} = reshape([cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)]',2,1,1,Nθ)

    logmfa::Array{Float64,3} = map(x -> (x>0 ? log(x) : logtol), fa);
    logmfp::Array{Float64,2} = map(x -> (x>0 ? log(x) : logtol), fp);
    p_rho ::Array{Float64,2} = p.(ρ) #functon p is labelled W in the pdf

    Ua::Array{Float64,4}  = -midpoint_bond_diff_θ(logmfa .+ p_rho; Nx=Nx, Nθ=Nθ).+ λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ) .+ λ*eθ 
    Up::Array{Float64,3}  = -midpoint_bond_diff(  logmfp  + p_rho; Nx=Nx       ) + λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx )
    Uθ::Array{Float64,3}  = -midpoint_Θ_diff(fa; Nx=Nx, Nθ = Nθ)

    return Ua, Up, Uθ
end


function F_fluxes(Ua::Array{Float64,4}, Up::Array{Float64,3}, Uθ::Array{Float64,3}, moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3}; Nx::Int64 =100, Nθ::Int64 =100 )
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

function time_stepper(fa::Array{Float64,3}, fp::Array{Float64,2}, δt::Float64; Nx::Int64 =100, Nθ::Int64 =100, λ::Float64 = 10., Dθ::Float64 = 10.)
    ρ::Array{Float64,2} = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    
    Ua::Array{Float64,4},   Up::Array{Float64,3},   Uθ::Array{Float64,3}   = U_velocities(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ)
    moba::Array{Float64,3}, mobp::Array{Float64,2}, mobθ::Array{Float64,3} = mob(fa,fp,ρ)
    Fa::Array{Float64,4},   Fp::Array{Float64,3},   Fθ::Array{Float64,3}   = F_fluxes(Ua, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)
    
    a::Float64 = maximum(abs.(Ua));
    b::Float64 = maximum(abs.(Up));
    c::Float64 = maximum(abs.(Uθ));
    
    tempu::Float64 = 1/(6*max(a*Nx, b*Nx, c*Nθ*Dθ/(2*π)));
    dt::Float64= min(δt, tempu)

    fa -= dt*( site_div_θ(Fa; Nx=Nx, Nθ=Nθ) + Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ))
    fp -= dt*site_div(Fp; Nx=Nx)

    return fa, fp, dt
end

function pde_stepper!(param::Dict{String,Any}, density::Dict{String,Any})
    @unpack δt, Nx, Nθ, Dθ, λ = param
    @unpack fa, fp, t = density
    fa, fp , dt = time_stepper(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)
    t::Float64 += dt

    @pack! density = fa, fp, t
    return dt
end

##

function run_pde_until!(param::Dict{String,Any},density::Dict{String,Any},T; save_on =false, max_steps = 100, save_interval = 1.)
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
                    dt = pde_stepper!(param,density)
                    time_since_save += dt 
                end
                if save_on
                    if save_interval < 1e-4
                        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                        safesave(filename,density)
                    elseif save_interval < 0.01
                        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 4))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                        safesave(filename,density)
                    else
                        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                        safesave(filename,density)
                    end
                end
            end
    else
        steps ::Int64 = 0
        while (density["t"] < T)&(steps < max_steps)
            time_since_save = 0.
            while (time_since_save < min(save_interval, T))&(steps < max_steps)
                dt = pde_stepper!(param,density)
                time_since_save += dt 
                steps += 1
            end
            if save_on
                @unpack t = density
                if save_interval < 1e-4
                    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
                    safesave(filename,density)
                elseif save_interval < 0.01
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

function perturb_pde!(param::Dict{String,Any}, density::Dict{String,Any}; δ = 0.01, pert = "n=2")
    @unpack Nx, S, ρa, ρp, λ, Dθ, Nx, Nθ, Dx, Pe, k  = param
    @unpack fa, fp = density
    ρ = ρa + ρp
    if ρ >0.9
        δ = min(δ, (1 - ρ)/(2*π+0.01));
    end
    #from stability: 
    α = π/2 -1
        ds  = self_diff(ρ)
        dsp = self_diff_prime(ρ)
        β = 1- ds;
        a = λ* ( (1-ρ-ds)/ρ - dsp) /2
    if pert == "n=2"
        γ = -8*π^2
        Λ = 1/2*(-Dθ + γ*β + sqrt(Dθ^2 + 2*Dθ*β*γ + (β*γ)^2 - 2*a^2*γ*ρa*ρ ) )
        C = 2*π*ρa*(ρa + ρp)/ (1 + Λ)
        Pa = (x,y,θ) -> δ*( ρa*cos(2*π*x/Nx)*cos(2*π*y/Nx) - C*(sin(2*π*x/Nx)*cos(2*π*y/Nx)*cos(2π*θ/Nθ) - sin(2*π*y/Nx)*cos(2*π*x/Nx)*sin(2π*θ/Nθ) ));
        Pp = (x,y) -> δ*( ρp*cos(2*π*x/Nx)*cos(2*π*y/Nx) );
    end
    if pert == "n=1"
        K = collect(0:1:(k-1))
        c,ω = ap_mstabparams_lite(ρa,ρp,Dx,Pe,Dθ)
        matrix = ap_MathieuMatrix(c, ω; k=k)
        a, A = ap_MathieuEigen(matrix)
        Pa = (x,y,θ) -> real.( dot(A[2:1:(k+1),k+1],cos.(θ*K*(2*π/Nθ)))*exp(im*x*ω/Nx) )
        Pp = (x,y) -> real.(A[1,1]*exp(im*x*ω/Nx));
    end
    if pert == "rand"
        Pa = (x,y,θ) -> ρa*(( rand() - 0.5 )/(ρa+0.01))/(2*π);
        Pp = (x,y) -> ρp*( rand() - 0.5 )/(ρp+0.01);
    end
    # 
    perta = zeros(Nx,Nx,Nθ)
    pertp = zeros(Nx,Nx)
    for x₁ in 1:Nx, x₂ in 1:Nx, θ in S
        perta[x₁, x₂, θ] += Pa(x₁, x₂, θ);
    end
    for x₁ in 1:Nx, x₂ in 1:Nx
        pertp[x₁, x₂] += Pp(x₁, x₂);
    end
    c = dist_from_unif(param, perta.+ρa/(2*π), pertp.+ρp)
    fa += δ*perta/c
    fp += δ*pertp/c
    @pack! density = fa, fp;
end

function perturb_pde_run(param)
    @unpack T, save_interval, max_steps, pert, δ = param
    density = initialize_density(param)
    perturb_pde!(param,density; pert = pert, δ = δ);
    run_pde_until!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)
end

##
#phase seperation metircs

function dist_from_unif(param, fa, fp)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
    return sqrt(2*π*sum( (fa .- ρa/(2*π) ).^2)/(Nx*Nx*Nθ)) + sqrt(sum( (fp .- ρp).^2/(Nx*Nx)))
end

function time_dist_from_unif(param, fa_saves, fp_saves)
    n = length(fa_saves)
    dist_saves = zeros(n)
    for i in 1:n
        dist_saves[i] = dist_from_unif(param, fa_saves[i], fp_saves[i])
    end
    return dist_saves
end

# loading

function load_pdes(param::Dict{String,Any},T; save_interval = 1., start_time = 0.)
    @unpack name, Nx, Nθ, λ, ρa, ρp, δt, Dθ = param
    t = start_time
    fa_saves = []
    fp_saves = []
    t_saves = []
    while t ≤ T
        try
            if save_interval < 1e-4
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";
            elseif save_interval < 0.01
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

function refresh_stab_data(;stabdata = Dict{String,Any}(), ρs = 0.05:0.05:0.95,   Dθ = 10., Nx = 50, Nθ = 20, λs = 5.:5.:100., name = "high_density_stability_v4", save_on = true, t_end = 1.0)
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ).jld2"

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
                param = pde_param(;name = name, Pe = λ/sqrt(Dθ) , ρa = ρ, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, save_interval = 0.01, max_steps = 1e8)
                t_saves, fa_saves, fp_saves = load_pdes(param,0.02; save_interval = 0.001, start_time = 0.0)

                stab_dsit0 = dist_from_unif(param, fa_saves[1], fp_saves[1])

                stab_dsit1 = dist_from_unif(param, fa_saves[2], fp_saves[2])

                if (stab_dsit0>stab_dsit1)&(λ ∉ lin_stab)
                    push!(lin_stab, λ)
                end
            catch
                    #break
            end
            try
                    param = pde_param(;name = name, Pe = λ/sqrt(Dθ) , ρa = ρ, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, save_interval = 0.01, max_steps = 1e8)
                    t_saves, fa_saves, fp_saves = load_pdes(param,t_end; save_interval = 0.001, start_time = (t_end-0.02))

                    stab_dsit = dist_from_unif(param, fa_saves[1], fp_saves[1])

                    if (stab_dsit>0.1)&(λ ∉ unstable)
                        push!(unstable,  λ)
                    elseif (stab_dsit<0.05)&(λ ∉ stable)&(λ ∉ lin_stab)
                        push!(stable, λ)
                    elseif (λ ∉ unsure)
                        push!(unsure, λ)
                    end
            catch
                        #break
            end
        end
        #calc ans
        λcrit = λsym(ρ; Dθ = Dθ)
        for λ ∈ λs
                if (λ ≤ λcrit)&(λ ∉ lin_stab)
                    push!(lin_stab,  λ)
                end
        end
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
            try 
                λminlin = maximum(lin_stab)
                for λ ∈ λs
                    if (λ ≤ λminlin)&(λ ∉ lin_stab)
                        push!(lin_stab,  λ)
                    end
                end
        catch
        end
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

function next_param(stabdata, param; λmax = 1e8, λ_step = 10.)
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

function run_stab_search(param)
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
        unfinished, λ_next, param = next_param(stabdata, param; λ_step = λ_step, λmax=λmax)
        if unfinished
                println("running λ = $(λ_next) ρ = $(ρa)")
                perturb_pde_run(param)
                stabdata = refresh_stab_data(;stabdata = stabdata,  ρs = [ρa],   Dθ = Dθ, Nx = Nx, Nθ = Nθ, λs = λs , name = name, save_on = false)
                println("finished λ = $(λ_next) ρ = $(ρa)")
        end
        i += 1 
        println("iteration = $(i), unfinished = $(unfinished)")
        println("while condition = $((i < max_runs)&unfinished))")
    end
    println("i am done")

    return ρa, stabdata["ρ = $(ρa)"]
end

function run_stab_search_stupid_mode(param)
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
    unfinished, λ_next, param = next_param(stabdata, param; λ_step = λ_step, λmax=λmax)
    while  (i < max_runs)&unfinished
        @unpack λ = param
        println("running λ = $(λ) ρ = $(ρa)")
        perturb_pde_run(param)
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





println("booted")