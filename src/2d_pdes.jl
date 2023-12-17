cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using TensorOperations, LinearAlgebra, Distributions

## utility
    function ρ_f(f::Array{Float64,3}; Nx::Int64 = 100, Nθ::Int64 = 20 )
        return sum(f[:,:,2:end],dims = 3)[:,:,1]*2*π/Nθ + f[:,:,1]
    end

    function mag_f(f::Array{Float64,3}; Nx::Int64 = 100, Nθ::Int64 = 20)
        eθ:: Array{Float64,4} = zeros(2,1,1,Nθ+1)
        eθ[1,1,1,2:end] = cos.((1:Nθ)*2π/Nθ)
        eθ[2,1,1,2:end] = sin.((1:Nθ)*2π/Nθ)
        f:: Array{Float64,4}  = reshape(f,1,Nx,Nx,Nθ+1)
        return sum(f.*eθ; dims = 4)[:,:,:,1]*2*π/Nθ
    end

#

## bond diff
    function bond_diff(f::Array{Float64,3}; Nx::Int64 = 100, Nθ::Int64 = 100, Lx::Float64 = 1.0) 

        grad_f::Array{Float64,4} = zeros(2, Nx, Nx, Nθ+1)

        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            y₂ ::Int64 = x₂
            grad_f[1,x₁,x₂,:] = Nx*( f[y₁,y₂,:] - f[x₁,x₂,:] )/Lx
            ## 2 direction
            y₁ = x₁
            y₂ = (x₂ +Nx)%Nx +1
            grad_f[2,x₁,x₂,:] = Nx*( f[y₁,y₂,:] - f[x₁,x₂,:] )/Lx
        end
        return grad_f
    end

    function bond_diff_ρ(ρ::Vector{Float64}; Nx::Int64 = 100, Nθ::Int64 = 100, Lx::Float64 = 1.0) 
        grad_ρ::Array{Float64,4} = zeros(2, Nx, Nx)

        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            y₂ ::Int64 = x₂
            grad_ρ[1,x₁,x₂] = Nx*( ρ[y₁,y₂] - ρ[x₁,x₂] )/Lx
            ## 2 direction
            y₁ = x₁
            y₂ = (x₂ +Nx)%Nx +1
            grad_ρ[2,x₁,x₂] = Nx*( ρ[y₁,y₂] - ρ[x₁,x₂] )/Lx
        end
        return grad_ρ
    end 

    function bond_av(f::Array{Float64,3}; Nx::Int64 = 100, Nθ::Int64 = 100, Lx::Float64 = 1.0) 
        av_f::Array{Float64,3}= zeros(2, Nx, Nx)

        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
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

    function bond_Θ_diff(f::Array{Float64,3}; Nx::Int64 = 100,  Nθ::Int64 = 100) 
        grad_f::Array{Float64,3} = zeros(Nx, Nx, Nθ+1)
    
        for θ in 1:Nθ
            ϕ = ((θ % Nθ) +1)
            grad_f[:,:,θ] += Nθ*( f[:,:,ϕ] - f[:,:,θ] )/(2*π)
        end
        return grad_f
    end
    #
    function site_div(f::Array{Float64,4}; Nx::Int64 = 100, Nθ::Int64 = 100, Lx::Float64 = 1.0) 

        div_f::Array{Float64,3} = zeros(Nx, Nx, Nθ+1)
    
        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            y₂ ::Int64 = x₂
            div_f[y₁,y₂,:] += Nx*( f[1,y₁,y₂,:] - f[1,x₁,x₂,:] )/Lx
            ## 2 direction
            y₁ = x₁
            y₂ = (x₂ +Nx)%Nx +1
            div_f[y₁,y₂,:] += Nx*( f[2,y₁,y₂,:] - f[2,x₁,x₂,:] )/Lx
        end
    
        return div_f
    end

    function site_div_ρ(f::Array{Float64,3}; Nx::Int64 = 100, Nθ::Int64 = 100, Lx::Float64 = 1.0) 

        div_f::Array{Float64,2} = zeros(Nx, Nx)
    
        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            y₂ ::Int64 = x₂
            div_f[y₁,y₂] += Nx*( f[1,y₁,y₂] - f[1,x₁,x₂] )/Lx
            ## 2 direction
            y₁ = x₁
            y₂ = (x₂ +Nx)%Nx +1
            div_f[y₁,y₂] += Nx*( f[2,y₁,y₂] - f[2,x₁,x₂] )/Lx
        end
    
        return div_f
    end

    function site_θ_diff(f::Array{Float64,4}; Nx::Int64 = 100,  Nθ::Int64 = 100) 

        div_f::Array{Float64,3} = zeros(Nx, Nx, Nθ+1)
    
        for θ in 1:Nθ
            ϕ ::Int64 = ((θ % Nθ) +1)
            div_f[:,:,ϕ] += Nθ*(f[3,:,:,ϕ]-f[3,:,:,θ])/(2*π)
        end
    
        return div_f
    end
    #
#

## computation
    function p(x::Float64;logtol = 1e-10, γ =0.)
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
    end
    
    function mob(f::Array{Float64,3}, ρ::Matrix{Float64}; Nx::Int64 = 100, Nθ::Int64 = 20)
        ds::Array{Float64,3}     = reshape(self_diff.(ρ),Nx,Nx,1)
        mobf::Array{Float64,4}  = zeros(2,Nx,Nx,Nθ+1)
        mobf[1,:,:,:] = f.*ds
        mobf[2,:,:,:] = f
        return mobf
    end
    
    function upwind(U::Float64, mb_down::Float64, mb_up::Float64)
        return (U   > 0. ? U.*mb_down  : U.*mb_up)
    end
#

# coefficients
    function self_diff(ρ::Float64;logtol::Float64 = 1e-10)
        α::Float64= π/2 -1;
        if ρ ≤  0.
            ρ = logtol
        elseif ρ>1.
            ρ = 1.
        end
        return ( 1-ρ).*( α*(2*α-1)/(2*α+1)*ρ^2 - α*ρ +1)
    end

    function self_diff_prime(ρ::Float64;logtol::Float64 = 1e-10)
        α::Float64= π/2 -1;
        if ρ ≤  0.
            ρ = logtol
        elseif ρ>1.
            ρ = 1.
        end
        return - ( α*(2*α-1)/(2*α+1)*ρ.^2 - α*ρ .+1) + ( -ρ .+1)*(2*α*(2*α-1)/(2*α+1)*ρ - α );
    end

    function coeff_s(rho::Float64,ds::Float64)
        return ((rho*ds)>0 ? (1-rho-ds)/(rho*ds) : 0)
    end

    function coeff_mag_s(f::Array{Float64,3},ρ::Matrix{Float64};Nx::Int64 = 100,Nθ::Int64 = 20 )
        m    ::Array{Float64,3} = mag_f(f; Nx = Nx, Nθ= Nθ);
        ds   ::Matrix{Float64} = self_diff.(ρ);
        s    ::Array{Float64,3} = reshape(coeff_s.(ρ,ds),1,Nx,Nx);
        mag_s::Array{Float64,3} = s.*m
        return mag_s
    end
#

## running pde 
    function new_pde_param(DT::Float64, v0::Float64, DR::Float64, Δx::Float64, Lx::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
        param::Dict{String, Any} = Dict{String,Any}()
        Nx::Int64 = Int64(Lx/Δx ÷ 1)
        @pack! param = DT, v0, DR, Δx, Lx, ϕa, ϕp, δt, δ, T , name, Nx, save_interval, save_on
        return param
    end

    function initiate_uniform_pde(ϕa::Float64, ϕp::Float64, Nx::Int64 = 100, Nθ::Int64 = 20)
        f::Array{Float64,3} = zeros(Nx,Nx,Nθ+1)
        f[:,:,2:end]= fill(ϕa/(2*π),(Nx,Nx,Nθ))
        f[:,:,1]  = fill(ϕp,(Nx,Nx))
        return f
    end

    function U_velocities(f::Array{Float64,3}, ρ::Matrix{Float64}; Nx::Int64 =100, Nθ::Int64 = 20, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10.,DR::Float64 = 1.)
        logtol::Float64 = log(1e-10);

        eθ:: Array{Float64,4} = zeros(2,1,1,Nθ+1)
        eθ[1,1,1,2:end] = cos.((1:Nθ)*2π/Nθ)
        eθ[2,1,1,2:end] = sin.((1:Nθ)*2π/Nθ)
    
        logmf::Array{Float64,3}  = map(x -> (x>0 ? log(x) : logtol), f);
        p_rho ::Array{Float64,2} = p.(ρ) #functon p is labelled W in the pdf
        
        U::Array{Float64,4} = zeros(3,Nx,Nx,Nθ+1)
        U[1:2,:,:,:]  = -DT*bond_diff(logmf .+ p_rho; Nx=Nx, Nθ=Nθ, Lx=Lx).+ v0*bond_av(coeff_mag_s(f,ρ; Nx=Nx, Nθ=Nθ ); Nx=Nx ) .+ v0*eθ 
        U[3,:,:,:]    = -DR*bond_Θ_diff(logmf; Nx=Nx, Nθ = Nθ)
    
        return U
    end
    
    function F_fluxes(U::Array{Float64,4}, mobf::Array{Float64,4}; Nx::Int64 =100, Nθ::Int64 =20)
        F ::Array{Float64,4} = zeros(3,Nx,Nx,Nθ+1);
        for x₁ in 1:Nx, x₂ in 1:Nx
            local y₁,y₂
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            y₂ ::Int64 = x₂
            F[1,x₁,x₂,:] = upwind.(U[1,x₁,x₂,:], mobf[1,x₁,x₂,:], mobf[1,y₁,y₂,:])
            ## 2 direction
            y₁ = x₁
            y₂ = (x₂ +Nx)%Nx +1
            F[2,x₁,x₂,:] = upwind.(U[2,x₁,x₂,:], mobf[1,x₁,x₂,:], mobf[1,y₁,y₂,:])
        end
        for θ in 1:Nθ
            local ϕ
            ϕ::Int64 = ((θ % Nθ) +1)
            F[3,:,:,θ] = upwind.(U[3,:,:,θ], mobf[2,:,:,θ], mobf[2,:,:,ϕ])
        end
        return F
    end

    function time_step!(t::Float64, f::Array{Float64,3}; δt::Float64 = δt, Nx::Int64 =100, Nθ::Int64 =20, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10., DR::Float64 = 1.0)
        ρ::Matrix{Float64} = ρ_f(f)
        
        U::Array{Float64,4}     = U_velocities(f,ρ; Nx=Nx, Nθ=Nθ, Lx=Lx, DT=DT, v0=v0, DR=DR);
        mobf::Array{Float64,4}  = mob(f,ρ; Nx=Nx, Nθ=Nθ);
        F::Array{Float64,4}     = F_fluxes(U, mobf; Nx=Nx, Nθ = Nθ);
        
        a::Float64 = maximum(abs.(U));
        
        tempu::Float64 = 1/(6*a*Nx);
        dt::Float64= min(δt, tempu);
    
        f -= dt*(  site_div(F; Nx=Nx, Nθ=Nθ, Lx=Lx) + site_θ_diff(F; Nx=Nx, Nθ=Nθ) ) 
        t += dt
        return t, f
    end
#

## saving funcitons 
    function pde_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, Δx, Nθ, Lx, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/2d_pdes_raw/$(name)/[DT,v0,DR,Δx,Nθ,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Nθ,Lx,ϕa,ϕp])/t=$(s).jld2"
    end

    function pde_time_series_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, Δx, Nθ, Lx, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/2d_pdes_pro/$(name)/[DT,v0,DR,Δx,Nθ,Lx,ϕa,ϕp]=$([DT,v0,DR,Δx,Nθ,Lx,ϕa,ϕp])/T=$(s)_Δt=$(save_interval).jld2"
    end

    function load_pde(param::Dict{String, Any},t::Float64)
        filename::String = pde_save_name(param::Dict{String, Any},t::Float64)
        data::Dict{String, Any} = load(filename)
        @unpack t, f = data 
        return t, f
    end

    function load_compress_pde(param::Dict{String, Any})
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, save_interval, save_on = param
        t_saves::Vector{Float64}            = []
        f_saves::Vector{Array{Float64,3}}  = []
        data::Dict{String, Any} = Dict()

        try
            try
                filename::String = pde_time_series_save_name(param,T)
                data = load(filename)
            catch
                filename::String = pde_time_series_save_name(param,T-save_interval)
                data= load(filename)
            end
            @unpack t_saves, f_saves = data
            println("fast load")
        catch
            println("full load")
            s = 0.
            t = 0.
            while s<T
                try
                    t, f = load_pde(param,s)
                    push!(f_saves,f)
                    push!(t_saves,t)
                    s += save_interval
                catch
                    s += save_interval
                end
            end
            if t > 0.
                filename::String = pde_time_series_save_name(param,t)
                data = Dict("f_saves" => f_saves, "t_saves" => t_saves)
                safesave(filename,data)
                println("saved")
            end
        end

        return t_saves, f_saves
    end

    function pde_vid_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Nθ, Lx, Ly, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/2d_pdes_vids/$(name)/[DT,v0,DR,N,Nθ,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Nθ,Lx,Ly,ϕa,ϕp])/t=$(s).mp4"
        pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/2d_pdes_vids/$(name)/[DT,v0,DR,N,Nθ,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Nθ,Lx,Ly,ϕa,ϕp])"
        return filename, pathname
    end
#

## pertubation

    function ap_MathieuMatrix(ρa,ρp,DT,v0,DR;k::Int64 = 10, Lx::Float64 = 2.0)
        ρ = ρa + ρp
        ω = 2*π/Lx
        #
        ds = self_diff(ρ) 
        dp = self_diff_prime(ρ)
        DD = (1-ds)/ρ
        #
        s = DD - 1
        p = -DT*ds*ω^2
        q = -v0*im*ω*ds/2 
        #
        matrix = Complex.(zeros(k+1, k+1))
        for u in 1:(k+1)
            for v in 1:(k+1)
                if abs(u - v) == 1
                    matrix[u, v] = q
                elseif u == v
                    matrix[u, v] = p - DR*(u-2)^2 
                end
            end
        end
        matrix[1, 1] = - DT*(ω^2)*(ds+ρp*DD)
        matrix[1, 2] = - 2*π*DT*(ω^2)*ρp*DD
        matrix[1, 3] = - π*v0*im*ω*ρp*s
        matrix[2, 1] = - DT*(ω^2)*(ρa/(2*π))*DD
        matrix[2, 2] = - DT*(ω^2)*(ds+ρa*DD)
        matrix[2, 3] = - v0*im*ω*(ρa*s+ds)/2
        matrix[3, 1] = - v0*im*ω*(ρa/(2*π))*dp
        matrix[3, 2] = - v0*im*ω*(ρa*dp+ds)
        return matrix
    end

    function lin_pert_values(param; k=10)
        @unpack DT, v0, DR, Lx, ϕa, ϕp = param
        W = ap_MathieuMatrix(ϕa,ϕp,DT,v0,DR; Lx = Lx, k = k)
        values,vectors = eigen(W)
        return values[k+1], vectors[:,k+1]
    end

    function dist_from_unif(f, param)
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt, δ = param
        sum = 0.
        sum += sum( (f[:,:,2:end] .- ϕa/(2*pi)).^2 ) 
        sum += sum( (f[:,:,1]     .- ϕp       ).^2 )
        return sqrt(sum)
    end

    function L2_norm(f; Nx = 100, Nθ = 20, Lx = 2)
        total = 0.
        total += Lx*Lx*2*π*sum( (f[:,:,2:end] ).^2 )/Nx/Nx/Nθ
        total += Lx*Lx*sum( (f[:,:,1]     ).^2 )/Nx/Nx
        return sqrt(total)
    end

    function perturb_pde!(f::Array{Float64,3}, param::Dict{String, Any}; type = "rand")
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, Nθ, save_interval, save_on, δt, δ = param
        pertf::Array{Float64,3}    = zeros(Nx,Nx,Nθ+1)
        if type == "rand"
            pertf    = rand(Uniform(-1,1),Nx,Nx,Nθ+1)
            if ϕp ==0.
                pertf[:,:,1] = zeros(Nx,Nx)
            end
        else
            k = 20
            value, vector = lin_pert_values(param; k = k)
            
            wave::Vector{ComplexF64}   = exp.((1:Nx)*(-im*2*π/Nx))

            K = collect(0:1:(k-1))

            for x in 1:Nx, y in 1:Nx
                pertf[x,y, 1] = real(vector[1]*wave[x])
                for θ in 1:Nθ
                    pertf[x,y,θ+1] = real( dot(vector[2:end],cos.(θ*K*(2*π/Nθ)))*wave[x] )
                end
            end
        end

        c = L2_norm(pertf; Nx=Nx, Nθ=Nθ, Lx=Lx)
        f += δ*pertf/c

        return f
    end
#

## solving functions 

    function run_new_pde(param::Dict{String, Any})
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, Nθ, save_interval, save_on, δt, pert = param
        # configuration
        f::Array{Float64,3} = initiate_uniform_pde(ϕa, ϕp, Nx, Nθ);
        f = perturb_pde!(f, param; type = pert);
        t::Float64 = 0.;
        s::Float64 = save_interval

        #inital save
        if save_on
            filename::String        = pde_save_name(param,t)
            data::Dict{String, Any} = Dict("f" => f, "t" => t)
            safesave(filename,data)
        end

        while t < T
            while t < s
                t, f = time_step!(t, f; δt=δt, Nx=Nx, Nθ, Lx=Lx, DT=DT, v0=v0, DR=DR);
            end
            #save snapshot
            if save_on
                filename    = pde_save_name(param,t)
                data        = Dict("f" => f, "t" => t)
                safesave(filename,data)
            end
            s += save_interval
        end
        return t, f
    end

    function run_current_pde(param::Dict{String, Any},dt::Float64, f::Array{Float64,3},t::Float64)
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, Nθ,save_interval, save_on, δt = param
        # configuration
        s::Float64 = t + save_interval
        t_end::Float64 = t+dt

        #inital save
        if save_on
            filename::String = pde_save_name(param,t)
            data::Dict{String, Any} = Dict("f" => f, "t" => t)
            safesave(filename,data)
        end

        while t < t_end
            while t < s
                t, f = time_step!(t, f; δt=δt, Nx=Nx, Nθ, Lx=Lx, DT=DT, v0=v0, DR=DR);
            end
            #save snapshot
            if save_on
                filename    = pde_save_name(param,t)
                data        = Dict("f" => f, "t" => t)
                safesave(filename,data)
            end
            s += save_interval
        end
        return t, f
    end

    function load_and_run_pde(param::Dict{String, Any})
        @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, Nθ, save_interval, save_on, δt = param
        # configuration
        f::Array{Float64,3} = initiate_uniform_pde(ϕa, ϕp, Nx, Nθ);
        t::Float64 = 0.;
        s::Float64 = T;
        loaded::Bool = false

        while s>0.
            try
                t, f = load_pde(param,s)
                loaded = true
                s = -1.
            catch
                loaded = false
                s += -save_interval
            end
        end

        if loaded
            println("load at t = $(t)")
            s = t + save_interval
            while t < T
                while t < s
                    t, f = time_step!(t, f; δt=δt, Nx=Nx, Nθ, Lx=Lx, DT=DT, v0=v0, DR=DR);
                end
                #save snapshot
                if save_on
                    filename    = pde_save_name(param,t)
                    data        = Dict("f" => f, "t" => t)
                    safesave(filename,data)
                end
                s += save_interval
            end
        else
            println("all loading failed; running new pde")
            t, f = run_new_pde(param)
        end
        return t, f
    end
#

println("v1.0")