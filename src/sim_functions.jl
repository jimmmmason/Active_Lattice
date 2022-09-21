cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using StatsBase, DataStructures, UnPack, LinearAlgebra


function uniform_initial_param(; name = "test", D =1. , λ =1. ,ρa = 0.1, ρp = 0.1, L=4, d=2)
    param = Dict{String,Any}()  
    #this is the only dimension dependent part:
    Ω = [[i,j] for i in 1:L for j in 1:L] 
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))

    function angles(x,n) 
        if n == 1
            return 2*π*rand()
        else
            return -1
        end
    end
    function rates(n,m,i)
        if m[1]>0
            return 0
        elseif n[1]==0
            return 0.
        elseif n[2]==2
            return L^2*D
        else
            return L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
        end
    end

    @pack! param = name, L, D, λ, ρa, ρp, Ω, E, site_distribution, angles, rates
    return param
end

function initialize(param::Dict{String,Any})
    @unpack name, L, λ, ρa, ρp, Ω, E, site_distribution, angles, rates = param
    # create configuration, rates and jumps
    η = fill([],(L,L))
    c = fill(0.,(L,L,4))
    j = fill([],(L,L,4))
    #fill configuration
    for x ∈ Ω
        local w, n
        w = Weights(site_distribution[x...])
        #fill model
        n = sample([0,1,2],w)
        η[x...] = [n, angles(x,n)]
    end
    #fill rates and jumps
    for x ∈ Ω 
        for i in 1:4
            local y
            #find adjacent site
            y  = (x + E[i] +[L-1,L-1]) .% L + [1,1]
            #fill rates 
            c[x...,i] = rates(η[x...],η[y...],i)
            #fill jump vectors
            j[x...,i] = [x,y]
        end
    end
    #pack into model
    t = 0.
    model = Dict{String,Any}() 
    @pack! model = η, c, j, t
    return model
end

function model_step!(param::Dict{String,Any},model::Dict{String,Any})
    @unpack name, L, λ, ρa, ρp, Ω, E, site_distribution, angles, rates = param
    @unpack η, c, j, t = model
    #select jump
    w     = Weights( [(c...)...])
    jump  = sample(j, w)
    # increase time
    t += randexp()/sum(c)
    #execute jumps
    η[jump[2]...] = deepcopy(η[jump[1]...])
    η[jump[1]...] = [0,-1]
    #correct propensity
    for x in jump
        for i in 1:4
            local y
            #find adjacent site
            y  = (x + E[i] +[L-1,L-1]) .% L + [1,1]
            #correct new rates 
            c[x...,i]   = rates(η[x...],η[y...],i  )
            c[y...,5-i] = rates(η[y...],η[x...],5-i)
        end
    end
    @pack! model = η, c, j, t
    return t
end

function run_model_until!(param::Dict{String,Any},model::Dict{String,Any},T; return_all = false)
    if return_all
        η_saves = []
        t_saves = []
        while model["t"] < T 
            model_step!(param,model)
            push!(η_saves,model["η"])
            push!(t_saves,t)
        end
        return t_saves, η_saves
    else
        while model["t"] < T 
            model_step!(param,model)
        end
    end
end

function run_model_intervals!(param::Dict{String,Any},model::Dict{String,Any},T; interval = 0.001)
    η_saves = []
    t_saves = []
    while model["t"] < T
        run_model_until!(param,model,model["t"]+interval; return_all = false)
        push!(η_saves, model["η"])
        push!(t_saves, model["t"])
    end 
end