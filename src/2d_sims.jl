cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using StatsBase, DrWatson, Random

## utility 
    function weight_index(i::Int64,j::Int64,k::Int64; Nx::Int64 = 1, Ny::Int64 = 1)
        return i + Nx*(j-1) + Nx*Ny*(k-1)
    end
#

#runnning simulation 
    using StatsBase, DataStructures, UnPack, LinearAlgebra, Random

    function uniform_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, Δt = 0.01, Dθ =10., T=1.0)
        param = Dict{String,Any}()   
        E = [[1,0],[0,1],[0,-1],[-1,0],]
        site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))
        POSITIONS = reshape(collect(1:(L^2*4)), (L,L,4))
        function angles(x,n)
            if n == 1
                return 2*π*rand()
            else
                return -1
            end
        end
        function rates(n,m,i)
            E = [[1,0],[0,1],[0,-1],[-1,0],]
            if m[1]>0.
                return 0.
            elseif n[1]==0.
                return 0.
            elseif n[2]==2.
                return L^2*D
            else
                return L^2*D + (L*λ*E[i]/2)⋅[cos(n[2]),sin(n[2])] 
            end
        end
        @pack! param = name, L, D, λ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ, T, POSITIONS
        return param
    end

    function initialize_model(param::Dict{String,Any})
        @unpack name, L, D, λ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ = param
        # create configuration, rates and jumps
        η = zeros(L,L,2)
        c = zeros(L,L,4)
        j = fill([],(L,L,4))
        #fill configuration
        for x₁ in 1:L, x₂ in 1:L
            local w, n
            w = Weights(site_distribution[x₁, x₂ ])
            #fill model
            n = sample([0,1,2],w)
            η[x₁, x₂, : ] = [n, angles([x₁, x₂],n)]
        end
        #fill rates and jumps
        for x₁ in 1:L, x₂ in 1:L 
            for i in 1:4
                local y
                #find adjacent site
                y₁, y₂  = ([x₁, x₂] + E[i] +[L-1,L-1]) .% L + [1,1]
                #fill rates 
                c[x₁, x₂ ,i] = rates(η[x₁, x₂,: ],η[y₁, y₂,: ],i)
                #fill jump vectors
                j[x₁, x₂ ,i] = [[x₁, x₂] ,[y₁, y₂]]
            end
        end
        #pack into model
        w     = weights( [(c...)...])
        α = sum(c)
        Δτ = Δt/α
        t = 0.
        model = Dict{String,Any}()
        @pack! model = η, w, j, t, α, Δτ
        @pack! param = name, L, D, λ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ
        return model
    end

    function get_rates(param, η)
        @unpack L, rates, Δt = param
        # create configuration, rates and jumps
        c = zeros(L,L,4)
        j = fill([],(L,L,4))
        #fill rates and jumps
        for x₁ in 1:L, x₂ in 1:L 
            for i in 1:4
                local y
                #find adjacent site
                y₁, y₂  = ([x₁, x₂] + E[i] +[L-1,L-1]) .% L + [1,1]
                #fill rates 
                c[x₁, x₂ ,i] = rates(η[x₁, x₂,: ],η[y₁, y₂,: ],i)
            end
        end
        #pack into model
        w     = weights( [(c...)...])
        α = sum(c)
        Δτ = Δt/α
        return w, j, α, Δτ
    end

    function model_step!(param::Dict{String,Any},model::Dict{String,Any})
        @unpack L, E, rates, Dθ, Δt, POSITIONS = param
        @unpack η, j, t, α, w = model
        if 2*Dθ*Δt> 0.05*2*π
            println("time step too large")
            Δt = 0
            t += 10e6
        end
        #run Gillespie algorithm until t+Δt
        δt::Float64 = 0.
        while δt < Δt
                #update total propensity
                α = sum(w)
                #select jump
                jump  = sample(j, w)
                #execute jump
                η[jump[2]...,:], η[jump[1]...,:] = η[jump[1]...,:], η[jump[2]...,:]
                #update time
                δt += randexp()/α
                #update propensity
                for (x₁, x₂) in jump
                    for i in 1:4
                            local y₁, y₂ 
                            #find adjacent site
                            y₁, y₂  = ([x₁, x₂] + E[i] +[L-1,L-1]) .% L + [1,1]
                            #correct new rates 
                            w[POSITIONS[x₁, x₂ ,i]]  = rates(η[x₁, x₂, : ],η[y₁, y₂, : ],i  )
                            w[POSITIONS[y₁, y₂ ,5-i]] = rates(η[y₁, y₂, : ],η[x₁, x₂, : ],5-i)
                    end
                end
        end
        #diffuse angles for interval t+δt
        for x₁ in 1:L, x₂ in 1:L
            if η[x₁,x₂,1] == 1.
                r = randn()
                η[x₁,x₂,2] = (η[x₁,x₂,2] + sqrt(2*Dθ*δt)*r + 2*π) % (2*π)
            end
        end
        t += δt
        @pack! model = η, w, t
        return t, model
    end
#

## saving funcitons 
    function sim_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/2d_sims_raw/$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])/t=$(s).jld2"
    end

    function sim_time_series_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        return "/store/DAMTP/jm2386/Active_Lattice/data/2d_sims_pro/$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])/T=$(s)_Δt=$(save_interval).jld2"
    end

    function load_sim(param::Dict{String, Any},t::Float64)
        filename::String = sim_save_name(param::Dict{String, Any},t::Float64)
        data::Dict{String, Any} = load(filename)
        @unpack η, t = data 
        return t, η
    end

    function load_compress_sim(param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, N₁, N₂, save_interval, save_on = param
        t_saves::Vector{Float64}            = []
        η_saves::Vector{Array{Float64, 3}}  = []
        data::Dict{String, Any} = Dict()

        try
            try
                filename::String = sim_time_series_save_name(param,T)
                data = load(filename)
            catch
                filename::String = sim_time_series_save_name(param,T-save_interval)
                data= load(filename)
            end
            @unpack t_saves, η_saves = data
            println("fast load")
        catch
            println("full load")
            s = 0.
            t = 0.
            while s<T
                try
                    t, η = load_sim(param,s)
                    push!(η_saves,η)
                    push!(t_saves,t)
                    s += save_interval
                catch
                    s += save_interval
                end
            end
            if t > 0.
                filename::String = sim_time_series_save_name(param,t)
                data = Dict("η_saves" => η_saves, "t_saves" => t_saves)
                safesave(filename,data)
                println("saved")
            end
        end

        return t_saves, η_saves
    end

    function sim_vid_save_name(param::Dict{String, Any},t::Float64)
        @unpack DT, v0, DR, N, Lx, Ly, ϕa, ϕp, name, save_interval = param
        s = round(t ; digits = Int64( -log10(save_interval) ÷ 1 ))
        filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/2d_sims_vids/$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])/t=$(s).mp4"
        pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/2d_sims_vids/$(name)/[DT,v0,DR,N,Lx,Ly,ϕa,ϕp]=$([DT,v0,DR,N,Lx,Ly,ϕa,ϕp])"
        return filename, pathname
    end
#

## simulation functions 

    function run_new_sim(param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, N₁, N₂, save_interval, save_on = param
        # configuration
        model = initialize_model(param);
        @unpack t, η = model
        s = save_interval

        #inital save
        if save_on
            filename::String        = sim_save_name(param,t)
            data::Dict{String, Any} = Dict("η" => η, "t" => t)
            safesave(filename,data)
        end

        while t < T
            while t < s
                t, model = model_step!(param, model);
            end
            #save snapshot
            if save_on
                @unpack t, η = model
                filename    = sim_save_name(param,t)
                data        = Dict("η" => η, "t" => t)
                safesave(filename,data)
            end
            s += save_interval
        end
        return t, η
    end

    function run_current_sim(param::Dict{String, Any},dt::Float64, model::Dict{String,Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp , name, N₁, N₂, save_interval, save_on = param
        # configuration
        @unpack t, η = model
        t_end   = t+dt
        s       = t +save_interval

        #inital save
        if save_on
            filename::String = sim_save_name(param,t)
            data::Dict{String, Any} = Dict("η" => η, "t" => t)
            safesave(filename,data)
        end

        while t < t_end
            while t < s
                t, model = model_step!(param, model);
            end
            #save snapshot
            if save_on
                @unpack t, η = model
                filename    = sim_save_name(param,t)
                data        = Dict("η" => η, "t" => t)
                safesave(filename,data)
            end
            s += save_interval
        end
        return t, η
    end

    function load_and_run_sim(param::Dict{String, Any})
        @unpack DT, v0, DR, N, Lx, ϕa, ϕp, T , name, N₁, N₂, save_interval, save_on = param
        # configuration
        model = initialize_model(param);
        @unpack t, η = model
        loaded::Bool = false
        s = 0.

        while s>0.
            try
                t, η = load_sim(param,s)
                w, j, α, Δτ = get_rates(param, η)
                @pack! model = η, w, j, t, α, Δτ
                loaded = true
                model = 
                s = -1.
            catch
                loaded = false
                #println("load failed at t = $(s)")
                s += -save_interval
            end
        end

        if loaded
            println("load at t = $(t)")
            s = t + save_interval
            while t < T
                while t < s
                    t, model = model_step!(param, model);
                end
                #save snapshot
                if save_on
                    @unpack t, η = model
                    filename    = sim_save_name(param,t)
                    data        = Dict("η" => η, "t" => t)
                    safesave(filename,data)
                end
                s += save_interval
            end
        else
            println("all loading failed; running new simulation")
            t, η = run_new_sim(param)
        end
        return t, η
    end
#

## proccessing funcitons 
    using PyPlot

    function site_ρ(ηx::Array{Float64,1})
        if ηx[1] == 1.
            return 1
        elseif ηx[1] == 2.
            return 1
        else
            return 0
        end
    end

    function local_density(η; N₁ = 100, N = 100, ϵ = 0.1) 
        r = Int64(N*ϵ ÷ 1)
        local_den = fill(0., ( N₁, N₁))
        E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
        S = length(E)
        for x₁ in 1: N₁, x₂ in 1: N₁
            for e ∈ E
                y₁, y₂  = ([x₁,x₂] + e +[ N₁-1, N₁-1]) .%  N₁ + [1,1]
                    local_den[x₁, x₂] += site_ρ(η[y₁, y₂,: ])/S
            end
        end
        return local_den
    end

    function local_passive_density(η; N₁ = 100, N = 100, ϵ = 0.1) 
        r = Int64(N*ϵ ÷ 1)
        local_den = fill(1., ( N₁, N₁))
        E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
        S = length(E)
        for x₁ in 1: N₁, x₂ in 1: N₁
            for e ∈ E
                y₁, y₂  = ([x₁,x₂] + e +[ N₁-1, N₁-1]) .%  N₁ + [1,1]
                if η[y₁, y₂,: ][1] == 2
                    local_den[x₁, x₂] += -1/S
                end
            end
        end
        return local_den
    end

    function local_active_density(η; N₁ = 100, N = 100, ϵ = 0.1) 
        r = Int64(N*ϵ ÷ 1)
        local_den = fill(1., ( N₁, N₁))
        E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
        S = length(E)
        for x₁ in 1: N₁, x₂ in 1: N₁
            for e ∈ E
                y₁, y₂  = ([x₁,x₂] + e +[ N₁-1, N₁-1]) .%  N₁ + [1,1]
                if η[y₁, y₂,: ][1] ==1
                    local_den[x₁, x₂] += -1/S
                end
            end
        end
        return local_den
    end

    function local_polarisation(η; N₁ = 100, N = 100, ϵ = 0.1) 
        r = Int64(N*ϵ ÷ 1)
        local_polarisatoin = fill(Complex(0.), (N₁,N₁))
        E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
        S = length(E)
        for x₁ in 1:N₁, x₂ in 1:N₁
            for e ∈ E
                y₁, y₂  = ([x₁,x₂] + e +[N₁-1,N₁-1]) .% N₁ + [1,1]
                if η[y₁, y₂,:][1] == 1
                    local_polarisatoin[x₁, x₂] += exp(η[y₁, y₂,:][2]*im)/S
                end
            end
        end
        return local_polarisatoin
    end

    function extract_points(η;  N₁ = 100, N = 100, ϵ = 0.1) 
        r = Int64(N*ϵ ÷ 1)
        passive = []
        active  = []
        directions  = []
        correction = [0.5,0.5]
        for x₁ in 1:N₁, x₂ in 1: N₁
            if η[x₁, x₂ ,1]== 2
                push!(passive,[x₁, x₂]-correction)
            elseif η[x₁, x₂ ,1]== 1
                push!(active,[x₁, x₂]-correction)
                push!(directions,η[x₁, x₂ ,2])
            end
        end
        
        passive = reshape([(passive...)...], (2,:))/N
        active  = reshape([(active...)...], (2,:))/N
        return passive, active, directions
    end
#

println("v1.1")