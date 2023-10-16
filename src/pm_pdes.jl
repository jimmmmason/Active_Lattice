cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using TensorOperations, LinearAlgebra

## utility 
    function midpoint_bond_diff_1d(f::Vector{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0) 

        grad_f::Vector{Float64} = zeros(Nx)

        for x₁ in 1:Nx
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            grad_f[x₁] = Nx*( f[y₁] - f[x₁] )/Lx
        end
        return grad_f
    end

    function midpoint_bond_diff_2d(f::Matrix{Float64}; Nx::Int64 = 100,  Lx::Float64 = 1.0) 

        grad_f::Matrix{Float64} = zeros(Nx, 3)

        for x₁ in 1:Nx
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            grad_f[x₁,:] = Nx*( f[y₁,:] - f[x₁,:] )/Lx
        end
        return grad_f
    end

    function midpoint_bond_av_1d(f::Vector{Float64}; Nx::Int64 = 100) 
        av_f::Vector{Float64}= zeros(Nx)

        for x₁ in 1:Nx
            ## 1 direction
            y₁::Int64 = (x₁ +Nx)%Nx +1
            av_f[x₁] = ( f[y₁] + f[x₁] )/2
        end
        return av_f
    end

    function site_div_1d(f::Vector{Float64}; Nx::Int64 = 100, Lx::Float64 = 1.0) 

        div_f::Vector{Float64} = zeros(Nx)

        for x₁ in 1:Nx
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            div_f[y₁] += Nx*( f[y₁] - f[x₁] )/Lx
        end
        return div_f
    end

    function site_div_2d(f::Matrix{Float64}; Nx::Int64 = 100,  Lx::Float64 = 1.0) 

        div_f::Matrix{Float64} = zeros(Nx, 3)

        for x₁ in 1:Nx
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            div_f[y₁,:] += Nx*( f[y₁,:] - f[x₁,:] )/Lx
        end

        return div_f
    end

    function flip_term(f::Matrix{Float64}; Nx::Int64 = 100)
        flip_term::Matrix{Float64} = zeros(Nx,3)
        flip_term[:,1] = f[:,1] - f[:,2]
        flip_term[:,2] = f[:,2] - f[:,1]
        return flip_term
    end

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
    
    function mob(f::Matrix{Float64}, ρ::Vector{Float64})
        ds::Vector{Float64} = self_diff.(ρ)
        return f.*ds
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

    function mag(f::Matrix{Float64})
        return f[:,2] - f[:,1]
    end

    function coeff_s(rho::Float64,ds::Float64)
        return ((rho*ds)>0 ? (1-rho-ds)/(rho*ds) : 0)
    end

    function coeff_mag_s(f::Matrix{Float64},ρ::Vector{Float64})
        m    ::Vector{Float64} = mag(f);
        ds   ::Vector{Float64} = self_diff.(ρ);
        s    ::Vector{Float64} = coeff_s.(ρ,ds);
        mag_s::Vector{Float64} = s.*m
        return mag_s
    end
#

## running pde 
    function new_pde_param(DT::Float64, v0::Float64, DR::Float64, N::Int64, Lx::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
        param::Dict{String, Any} = Dict{String,Any}()
        Nx::Int64 = Int64(Lx*N ÷ 1)
        @pack! param = DT, v0, DR, N, Lx, ϕa, ϕp, δt, δ, T , name, Nx, save_interval, save_on
        return param
    end

    function initiate_uniform_pde(ϕa::Float64, ϕp::Float64, Nx::Int64 = 100)
        f::Matrix{Float64}     = zeros(Nx,3)
        f[:,1:2]= fill(ϕa/2,(Nx,2))
        f[:,3]  = fill(ϕp,(Nx))
        return f
    end

    function U_velocities(f::Matrix{Float64}, ρ::Vector{Float64}; Nx::Int64 =100, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10.)
        logtol::Float64 = log(1e-10);
    
        eθ:: Array{Float64,2} = reshape([-1 1 0],1,3) 
    
        logmf::Matrix{Float64} = map(x -> (x>0 ? log(x) : logtol), f);
        p_rho ::Vector{Float64} = p.(ρ) #functon p is labelled W in the pdf
    
        U::Matrix{Float64}  = -DT*midpoint_bond_diff_2d(logmf .+ p_rho; Nx=Nx, Lx=Lx) .+ v0*midpoint_bond_av_1d(coeff_mag_s(f,ρ); Nx =Nx ) .+ v0*eθ 

        return U
    end
    
    function F_fluxes(U::Matrix{Float64}, moba::Matrix{Float64}; Nx::Int64 =100)
        F ::Matrix{Float64} = zeros(Nx,3);

        for x₁ in 1:Nx
            local y₁
            ## 1 direction
            y₁ ::Int64 = (x₁ +Nx)%Nx +1
            F[x₁,:] = upwind.(U[x₁,:], moba[x₁,:], moba[y₁,:])
        end
        return F
    end

    function time_step!(t::Float64, f::Matrix{Float64}; δt::Float64 = δt, Nx::Int64 =100, Lx::Float64 = 1.0, DT::Float64 = 1.0, v0::Float64 = 10., DR::Float64 = 1.0)
        ρ::Vector{Float64} = sum(f; dims =2)[:,1];
        
        U::Matrix{Float64}     = U_velocities(f,ρ; Nx=Nx, Lx=Lx, DT=DT, v0=v0);
        mobf::Matrix{Float64}  = mob(f,ρ);
        F::Matrix{Float64}     = F_fluxes(U, mobf; Nx=Nx);
        
        a::Float64 = maximum(abs.(U));
        
        tempu::Float64 = 1/(6*a*Nx);
        dt::Float64= min(δt, tempu);
    
        f -= dt*( site_div_2d(F; Nx=Nx, Lx=Lx) + DR*flip_term(f; Nx=Nx) )
        t += dt
        return t
    end
#

## saving funcitons 
    function pde_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/pm_pdes_raw/$(name)/[DT,v0,DR,N,Lx,ϕa,ϕp]=$([DT,v0,DR,N,Lx,ϕa,ϕp])/t=$(s).jld2"
    end

    function pde_time_series_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/pm_pdes_pro/$(name)/[DT,v0,DR,N,Lx,ϕa,ϕp]=$([DT,v0,DR,N,Lx,ϕa,ϕp])/t=$(s).jld2"
    end

    function load_pde(param::Dict{String, Any},t::Float64)
        filename::String = pde_save_name(param::Dict{String, Any},t::Float64)
        data::Dict{String, Any} = load(filename)
        @unpack t, f = data 
        return t, f
    end

    function load_compress_pde(param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, N, save_interval, save_on = param
        t_saves::Vector{Float64}            = []
        f_saves::Vector{Matrix{Float64}}  = []
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
#

## pertubation
    function lin_pert_values(param)
        @unpack DT, v0, DR, Lx, ϕa, ϕp = param
        ω = 2*π/Lx;
        Pe = v0;
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

    function dist_from_unif(f, param)
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt, δ = param
        return sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2   ))
    end

    function perturb_pde!(f::Matrix{Float64}, param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt, δ = param
    
        ω, value, vector = lin_pert_values(param)

        wave::Vector{ComplexF64}   = exp.((1:Nx)*(im*2*π/Nx))
        pertf::Matrix{Float64}     = zeros(Nx,3)

        pertf[:,1] = real.( wave*(vector[1]- vector[2])/2 )
        pertf[:,2] = real.( wave*(vector[1]+ vector[2])/2 ) 
        pertf[:,3] = real.( wave*(vector[3]) )

        c = dist_from_unif(f,param)
        f += δ*pertf/c

        return f
    end
#

## complete functions 

    function run_new_pde(param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt = param
        # configuration
        f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
        perturb_pde!(f, param);
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
                t = time_step!(t, f; δt=δt, Nx=Nx, Lx=Lx, DT=DT, v0=v0, DR=DR);
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

    function run_current_pde(param::Dict{String, Any},dt::Float64, f::Matrix{Float64},t::Float64)
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt = param
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
                t = time_step!(t, f; δt=δt, Nx=Nx, Lx=Lx, DT=DT, v0=v0, DR=DR);
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
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt = param
        # configuration
        f::Matrix{Float64} = initiate_uniform_pde(ϕa, ϕp, Nx);
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
                    t = time_step!(t, f; δt=δt, Nx=Nx, Lx=Lx, DT=DT, v0=v0, DR=DR);
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

## proccessing funcitons 
    # function f_to_rgb(f; type = "rho" )
    #     Nx, Ny = size(ρa)
    #     rgb_image = ones(Ny, Nx,3)
    #     if type == "rho"
    #         rgb_image[:,:,3] = -ρa'[Ny:-1:1,:]  .+1
    #         rgb_image[:,:,1] = -ρp'[Ny:-1:1,:] .+1
    #         rgb_image[:,:,2] = -ρa'[Ny:-1:1,:]  -ρp'[Ny:-1:1,:]  .+1
    #     else
    #         rgb_image[:,:,1] = map((x ->  max(0,x)), m'[Ny:-1:1,:])
    #         rgb_image[:,:,2] = map((x -> -min(0,x)), m'[Ny:-1:1,:])
    #         rgb_image[:,:,3] = -ρa'[Ny:-1:1,:] -ρp'[Ny:-1:1,:] .+1
    #     end
    #     return rgb_image
    # end
#

println("v1.0")