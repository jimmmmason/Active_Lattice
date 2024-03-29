cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")
##
#runnning simulation 
using StatsBase, DataStructures, UnPack, LinearAlgebra, Random

function uniform_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, d=2, Δt = 0.01, Dθ =10., T=1.0, γ = 0.)
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
            return L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
        end
    end
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ, T, POSITIONS
    return param
end

function extra_mixing_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=10, d=2, Δt = 0.001, γ = 0.01, Dθ =10., T=1.0)
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
        if n[1]== 0.
            return 0.
        elseif (n[1]==2.)&(m[1]==0.)
            return L^2*(D+γ) 
        elseif (n[1]==1.)&(m[1]==0.)
            return L^2*(D+γ) + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
        elseif (m[1]> 0.)
            return L^2*γ
        else
            return 69.
        end
    end
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ, T, POSITIONS
    return param
end

function initialize_model(param::Dict{String,Any})
    @unpack name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ = param
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
    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ
    return model
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
end

function run_model_until!(param::Dict{String,Any},model::Dict{String,Any},T; save_on =false)
    if save_on
        @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
        @unpack η, t = model
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
        data = Dict{String,Any}();
        @pack! data = η, t
        safesave(filename,data)
    end
    while model["t"] < T 
        model_step!(param,model)
        if save_on
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            @unpack η, t = model
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            data = Dict{String,Any}();
            @pack! data = η, t
            safesave(filename,data)
        end
    end
end

function run_and_dump_sim(param::Dict{String,Any},model::Dict{String,Any},T; dump_interval = 0.01, save_on =false)
    while model["t"] < T
        local t
        t = deepcopy(model["t"]+dump_interval)
        run_model_until!(param,model,t; save_on =false)
        if save_on
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            @unpack η, t = model
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            data = Dict{String,Any}();
            @pack! data = η, t
            safesave(filename,data)
        end
    end
end

function run_sim(param)
    @unpack T, Δt, χ, ρ, Pe = param
    model = initialize_model(param);
    println("starting run: χ = $(param["χ"]), ρ = $(param["ρ"]), Pe = $(param["Pe"])")
    s = 50*Δt
    while s<T 
        try
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            η, t = wload(filename, "η", "t")
            c = zeros(L,L,4)
            E = [[1,0],[0,1],[0,-1],[-1,0],]
            for x₁ in 1:L, x₂ in 1:L 
                for i in 1:4
                    #find adjacent site
                    y₁, y₂  = ([x₁, x₂] + E[i] +[L-1,L-1]) .% L + [1,1]
                    #fill rates 
                    c[x₁, x₂ ,i] = rates(η[x₁, x₂,: ],η[y₁, y₂,: ],i)
                end
            end
            #pack into model
            w    = weights( [(c...)...])
            @pack! model = η, w, t 
            println("loaded t = $(t)")
        catch
            println("load t = $(s) failed")
            s = T
        end
        s+= 50*Δt
    end
    run_model_until!(param, model, T; save_on = true);
    println("success: χ = $(param["χ"]), ρ = $(param["ρ"]), Pe = $(param["Pe"])")
end

function run_sim_dump(param)
    @unpack T, Δt = param
    dump_interval = 0.01
    model = initialize_model(param);
    println("starting run: γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])")
    while model["t"] < T
        local t
        t = deepcopy(model["t"]+dump_interval)
        run_model_until!(param,model,t; save_on =false)
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            @unpack η, t = model
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            data = Dict{String,Any}();
            @pack! data = η, t
            safesave(filename,data)
    end
    println("success: γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])")
end



function load_dump_sim_run(param)
    @unpack  T, Δt = param
    dump_interval = 0.01
    s = T
    i = 0
    η = []
    j = 0
    t = 0.
    while (0 ≤ s)&(j==0.) 
        try
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            η, t = wload(filename, "η", "t")
            println("load success at t = $(s)")
            i = 1
            j = 1
        catch
            println("load failed at t = $(s)")
        end
        s -= 0.1
    end
    if i==1
        println("loading success")
        model = Dict{String,Any}()
        @unpack name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ = param
        c = zeros(L,L,4)
        j = fill([],(L,L,4))
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
        @pack! model = η, w, j, t, α, Δτ
        while model["t"] < T
            local t
            t = deepcopy(model["t"]+dump_interval)
            run_model_until!(param,model,t; save_on =false)
                @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
                @unpack η, t = model
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
                data = Dict{String,Any}();
                @pack! data = η, t
                safesave(filename,data)
        end
    else
        println("loading failed")
    end
end

function load_dump_sim_run_home_save(param)
    @unpack  T, Δt = param
    dump_interval = 0.01
    s = T
    i = 0
    η = []
    j = 0
    t = 0.
    while (0 ≤ s)&(j==0.) 
        try
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            η, t = wload(filename, "η", "t")
            println("load success at t = $(s)")
            i = 1
            j = 1
        catch
            println("load failed at t = $(s)")
        end
        s -= 0.1
    end
    if i==1
        println("loading success")
        model = Dict{String,Any}()
        @unpack name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ = param
        c = zeros(L,L,4)
        j = fill([],(L,L,4))
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
        @pack! model = η, w, j, t, α, Δτ
        while model["t"] < T
            local t
            t = deepcopy(model["t"]+dump_interval)
            run_model_until!(param,model,t; save_on =false)
                @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
                @unpack η, t = model
                filename = "/home/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
                data = Dict{String,Any}();
                @pack! data = η, t
                safesave(filename,data)
        end
    else
        println("loading failed")
    end
end




function run_sim_old(param)
    @unpack T, Δt = param
    model = initialize_model(param);
    println("starting run: γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])")
    run_model_until!(param, model, T; save_on = true);
    println("success: γ = $(param["γ"]), ρ = $(param["ρa"]), λ = $(param["λ"])")
end

function load_etas(param::Dict{String,Any},T; dump_interval = 0.01, start_time= 0.)
    s = start_time
    t_saves = []
    η_saves = []
    while s < T
        try 
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
            η, t = wload(filename, "η", "t")
            push!(t_saves,t)
            push!(η_saves,η)
        catch
        end
        s += dump_interval
    end
    return t_saves, η_saves
end

function load_etas_1(param::Dict{String,Any},T; dump_interval = 0.01, start_time= 0.)
    s = start_time
    t_saves = []
    η_saves = []
    while s < T
        try 
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)_#1.jld2";
            η, t = wload(filename, "η", "t")
            push!(t_saves,t)
            push!(η_saves,η)
        catch
        end
        s += dump_interval
    end
    return t_saves, η_saves
end

function load_etas_av(param::Dict{String,Any},T; dump_interval = 0.01, start_time= 0., trials = 10)
    s = start_time
    t_saves = []
    η_saves = []
    no_abort = true
    while s < T
        try 
            t = s
            @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
            filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)_#1.jld2";
            η, t = wload(filename, "η", "t")
            push!(t_saves,t)
            push!(η_saves,η)
        catch
        end
        for i in 1:(trials-1)
            try 
                t = s
                @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
                filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 3))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)_#$(i).jld2";
                η, t = wload(filename, "η", "t")
                push!(t_saves,t)
                push!(η_saves,η)
            catch
                println("only $(i) trials at t = $(s), aborting...")
                s = T
                no_abort = false
            end
        end
        s += dump_interval
    end
    if no_abort
        l = Int64(length(t_saves)/(trials))
        t_saves = reshape(t_saves,(trials,l))
        η_saves = reshape(η_saves,(trials,l))
    end
    return t_saves, η_saves
end
##
#visualise data
using PyPlot, PyCall
@pyimport matplotlib.animation as anim

function animate_etas(param,t_saves,η_saves)
    frames = length(η_saves)-1
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    makeframe(i) = plot_eta(fig,ax,param, t_saves[i+1], η_saves[i+1])
    interval = Int64(round(10000/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    @unpack name, L, λ, γ, ρa, ρp = param
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/$(name)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)"
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/$(name)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)/start_time=$(round(t_saves[1]; digits=5))_end_time=$(round(t_saves[frames+1]; digits=5))_interval=$(interval)_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ).mp4"
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end

function make_sim_vid(param)
    @unpack T = param
    t_saves, η_saves = load_etas(param, T; dump_interval = 0.01)
    animate_etas(param,t_saves,η_saves)
end

function make_sim_vid_lite(param)
    @unpack T = param
    t_saves, η_saves = load_etas(param, T; dump_interval = 0.01)
    animate_etas(param,t_saves,η_saves)
end

function plot_eta(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Float64,3}; title = true, r = 3)
    @unpack name, L, λ, γ, ρa, ρp, E, Δt, site_distribution, angles, rates = param
    ax.clear()
    #collect data
    #y = [Int64(L/2),Int64(L/2)]-indexmin(local_polarisation(η, L; r = 3))
    #η = translate_η(η, L, y)
    passive, active, directions, polarisations = extract_points_v2(η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    Φ = round(translation_invariance(η;L= L); digits =4)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        if title
            ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t), Φ = $(Φ)")
        else
            #ax.set_title("Φ = $(Φ)")
        end
    # Plot points
    densities = reshape(local_density(η,L; r=r)',L^2)/(2*r+1)^2  #reshape(local_density(η,L; r=r),L^2) #local_polarisation(η,L; r=r)
    xs = []
    ys = []
    for x₁ in 1:L, x₂ in 1:L
        push!(xs, (x₁-0.5)/L)
        push!(ys, (x₂-0.5)/L)
    end
    ax.scatter(xs,ys, cmap = bone(),
            c = densities,
            s = 75, 
            marker = "s", 
            linewidths = 0.,
            alpha=0.5,
        )   
    ax.quiver(active[1,:],active[2,:], dx,dy, cmap = hsv(),
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

function plot_eta_old(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Float64,3}; title = true)
    @unpack name, L, λ, γ, ρa, ρp,  E, Δt, site_distribution, angles, rates, Pe, ρ, χ = param
    ax.clear()
    #collect data
    passive, active, directions, polarisations = extract_points(η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        if title
            if χ == 1.0
                ax.set_title("ϕ = $(ρ), Pe = $(Pe), t = $(round(t; digits =3)), Φ = $(Φ)")
            else 
                ax.set_title("ϕ = $(ρ),  χ = $(χ), Pe = $(Pe), t = $(t), Φ = $(Φ)")
            end
        else
            ax.set_title("Φ = $(Φ)")
        end
    # Plot points
    ax.scatter(active[1,:],active[2,:], cmap = gray(),
            c = polarisations/maximum(polarisations),
            s = 75, 
            marker = "s", 
            linewidths = 0.,
            alpha=0.5,
        )   
    ax.quiver(active[1,:],active[2,:], dx,dy, cmap = hsv(),
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

function plot_eta_old_old(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Float64,3}; title = true)
    @unpack name, L, λ, γ, ρa, ρp,  E, Δt, site_distribution, angles, rates, Pe, ρ, χ = param
    ax.clear()
    #collect data
    passive, active, directions, polarisations = extract_points(η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =3)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        if title
            if χ == 1.0
                ax.set_title("ϕ = $(ρ), Pe = $(Pe), t = $(round(t; digits =3)), Φ = $(Φ)")
            else 
                ax.set_title("ϕ = $(ρ),  χ = $(χ), Pe = $(Pe), t = $(t), Φ = $(Φ)")
            end
        else
            ax.set_title("Φ = $(Φ)")
        end
    # Plot points
    ax.quiver(active[1,:],active[2,:], dx,dy, cmap = hsv(),
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

function plot_eta_polarisation(fig::Figure, ax::PyObject, param::Dict{String,Any}, t::Float64, η::Array{Float64,3}; title = true)
    @unpack name, L, λ, γ, ρa, ρp,  E, Δt, site_distribution, angles, rates, Pe, ρ, χ = param
    ax.clear()
    #collect data
    passive, active, directions, polarisations = extract_points(Ω,η,L)
    dx = cos.(directions)
    dy = sin.(directions)
    t = round(t; digits=5)
    Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =3)
    #figure configuration
    ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.axis([0., 1., 0., 1.])
        ax.set_aspect("equal")
        if title
            if χ == 1.0
                ax.set_title("ϕ = $(ρ), Pe = $(Pe), t = $(round(t; digits =3)), Φ = $(Φ)")
            else 
                ax.set_title("ϕ = $(ρ),  χ = $(χ), Pe = $(Pe), t = $(t), Φ = $(Φ)")
            end
        else
            ax.set_title("Φ = $(Φ)")
        end
    # Plot points
    ax.quiver(active[1,:],active[2,:], dx,dy, cmap = gray(),
        polarisations,
        scale_units = "x",
        pivot = "mid",
        minlength = 0.1,
        minshaft = 1,
        width =1/L,
        headlength = 5,
        headaxislength=5,
        scale = L,
        alpha = 0.5,
    )
    ax.errorbar(passive[1,:],passive[2,:], 
        markersize = 400/L, 
        fmt= "o", 
        color = "black",
        alpha=0.8,
    )
    return fig
end

function extract_points(η::Array{Float64,3},L; r=3)
    passive = []
    active  = []
    directions  = []
    polarisations = []
    correction = [0.5,0.5]
    polarisation = local_polarisation(η,L; r=r)
    for x₁ in 1:L, x₂ in 1:L
        if η[x₁, x₂ ,1]== 2
            push!(passive,[x₁, x₂]-correction)
        elseif η[x₁, x₂ ,1]== 1
            push!(active,[x₁, x₂]-correction)
            push!(directions,η[x₁, x₂ ,2])
            push!(polarisations, polarisation[x₁, x₂ ])
        end
    end
    polarisations = (- polarisations .+ maximum(polarisations) )/maximum(polarisations)
    passive = reshape([(passive...)...], (2,:))/L
    active  = reshape([(active...)...], (2,:))/L
    return passive, active, directions, polarisations
end

function extract_points_v2(η::Array{Float64,3},L; r=3)
    passive = []
    active  = []
    directions  = []
    polarisations = []
    correction = [0.5,0.5]
    polarisation = local_polarisation(η,L; r=r)
    for x₁ in 1:L, x₂ in 1:L
        if η[x₁, x₂ ,1]== 2
            push!(passive,[x₁, x₂]-correction)
        elseif η[x₁, x₂ ,1]== 1
            push!(active,[x₁, x₂]-correction)
            push!(directions,η[x₁, x₂ ,2])
        end
        push!(polarisations, polarisation[x₁, x₂ ])
    end
    polarisations = (- polarisations .+ maximum(polarisations) )/maximum(polarisations)
    passive = reshape([(passive...)...], (2,:))/L
    active  = reshape([(active...)...], (2,:))/L
    return passive, active, directions, polarisations
end
##
#phase seperation metircs

function site_ρ(ηx::Array{Float64,1})
    if ηx[1] == 1.
        return 1
    elseif ηx[1] == 2.
        return 1
    else
        return 0
    end
end

function site_θ(ηx::Array{Float64,1})
    if ηx[1] == 1.
        return ηx[2]
    elseif ηx[1] == 2.
        return false
    else
        return false
    end
end

function fourier_config(η::Array{Float64,3};L::Int64 = 1, k::Vector{Float64}=[0,0])
    ϕLL = 0.
    for x₁ in 1:L, x₂ in 1:L
        ϕLL += ( 1-site_ρ(η[x₁, x₂,: ]))*exp(- im* k⋅ [x₁,x₂] )
    end
    return ϕLL/L^2
end

function translation_invariance(η::Array{Float64,3};L::Int64 = 1)
    return norm(fourier_config(η; L= L, k =[2*π/L,0.]))+norm(fourier_config(η; L= L, k =[0.,2*π/L]))
end

function density_hist(param::Dict{String,Any}, η::Array{Float64,3}; r = 2) 
    @unpack name, L, λ, γ, ρa, ρp,  E, Δt, site_distribution, angles, rates = param
    local_density = []
    B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]#if (i^2+j^2) ≤ r^2]
    n = length(B)
    for x₁ in 1:L, x₂ in 1:L
        local ρr 
        ρr = 0.
        for e ∈ B
            y₁, y₂  = ([x₁, x₂] + e +[L-1,L-1]) .% L + [1,1]
            ρr += site_ρ(η[y₁, y₂,: ])/n
        end
        push!(local_density,ρr) 
    end
    return local_density
end

function pol_hist(param::Dict{String,Any}, η::Array{Float64,3}; r = 2) 
    @unpack name, L, λ, γ, ρa, ρp,  E, Δt, site_distribution, angles, rates = param
    local_density = []
    B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
    for x₁ in 1:L, x₂ in 1:L
        local ρr
        ρr = [0., 0.]
        for e ∈ B
            y₁, y₂  = ([x₁, x₂] + e +[L-1,L-1]) .% L + [1,1]
            θ= site_θ(η[y₁, y₂,: ])
            if θ == false
            else
                ρr += [cos(θ), sin(θ)]/(2*r +1)^2
            end
        end
        push!(local_density,ρr) 
    end
    return local_density
end

function time_density_hist(fig::Figure, ax::PyObject, param::Dict{String,Any}, t_saves, η_saves; r = 3, bins = 3, name= "label", add_label = false)
    h = [];
    for η ∈ η_saves
        append!(h, density_hist(param, η; r = r));
    end
    edges = collect((-1/(2*numbins)):(1/(numbins)):(1+(1/(2*numbins))))
    if add_label
        ax.hist(h; bins = edges, histtype = "step", density = true, label = name)
    else
        ax.hist(h; bins = edges, histtype = "step", density = true)
    end
    ax.xaxis.set_ticks(0:0.5:1)
end

function local_density(η, L; r = 3) 
    local_den = fill(1., (L,L))
    E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
    S = length(E)
    for x₁ in 1:L, x₂ in 1:L
        for e ∈ E
            y₁, y₂  = ([x₁,x₂] + e +[L-1,L-1]) .% L + [1,1]
                local_den[x₁, x₂] += -site_ρ(η[y₁, y₂,: ])/S
        end
    end
    return local_den
end


function local_polarisation(η, L; r = 3) 
    local_polarisatoin = fill(Complex(0.), (L,L))
    E = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
    for x₁ in 1:L, x₂ in 1:L
        for e ∈ E
            y₁, y₂  = ([x₁,x₂] + e +[L-1,L-1]) .% L + [1,1]
            if site_θ(η[y₁, y₂,:]) == false
            else
                local_polarisatoin[x₁, x₂] += exp(site_θ(η[y₁, y₂,: ])*im)
            end
        end
    end
    return abs.(local_polarisatoin)
end

function translate_η(η, L, y)
    η2 = zeros(L,L,2)
    for x₁ in 1:L, x₂ in 1:L
        x = [x₁,x₂]
        x2 = (x + y +[L-1,L-1]) .% L + [1,1]
        η2[x2...,:] =  η[x₁, x₂,: ]
    end
    return η2
end 

function indexmin(A)
    n = length(size(A))
    x = argmin(A)
    X = []
    for i in 1:n
        push!(X,x[i])
    end 
    return X
end
#

println("booted")

#Example 
#=
#Parameters
param = uniform_initial_param(L=32, λ = 16, ρa = 0.9, ρp = 0.0, Δt = 0.01)
param = extra_mixing_initial_param(name = "test10", L=32, λ = 16, ρa = 0.7, ρp = 0.0, Δt = 0.01, γ = 0.00, T = 2.2)

run_sim(param)
model = initialize_model(param)
run_model_until!(param,model,0.02; save_on = false)
@profview model_step!(param, model)
#expand variables
@unpack name, L, λ, γ, ρa, ρp,  E, site_distribution, angles, rates = param
@unpack η, w, j, t = model
#Run and save
using BenchmarkTools
T = 0.05
@time model_step!(param, model);
run_model_until!(param,model, T; save_on =true)
#Loading
T = 1.0
t_saves, η_saves = load_etas(param, T; dump_interval = 0.01);
#plotting
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = length(t_saves)
plot_eta(fig,ax,param, t_saves[n], η_saves[n])
display(fig)
#video
animate_etas(param,t_saves,η_saves)
#symmetry
y = translation_invariance.(η_saves;Ω = Ω,L= L)
n = length(t_saves)
translation_invariance(η_saves[n];Ω = Ω,L= L)
clf()
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(y)
display(fig)
clf()
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = 6000
plot_eta(fig,ax,param, t_saves[n], η_saves[n])
display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_eta(fig,ax,param, t, η)
display(fig)
#
"/store/DAMTP/jm2386/Active_Lattice/data/sims_raw/$(name)/size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ)/time=$(round(t; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Dθ=$(Dθ).jld2";
#

fig, ax = PyPlot.subplots(figsize =(10, 10))
h =  density_hist(param, η; r= 5)
ax.hist(h; bins = 1000)
fig
=#

