cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using StatsBase, DataStructures, UnPack

function initialize(site_distribution::Array{Float64,2}, angles::Array{Float64,2}, rates::Array{} ; L=10,)
    η = fill(0,(L,2))
    c = zeros(((L,2),4))
    j = fill([],(L,4))
    Ω = [[i,j] for i in 1:L for j in 1:L] 
    E = [[1,0],[0,1],[-1,0],[0,-1],]

    for x ∈ Ω
        local w, n
        w = Weights(site_distribution)
        #fill model
        n = sample([0,1,2],w)
        η[x...] = [n, angles[x...]]
    end
    for x ∈ Ω
        for i in 1:4
        #find adjacent sites
        y  = (x + E[i] -[1,1]) .% L + [1,1]
        #fill rates
        c[x...,i]   = rates[η[x...]+1,η[y...]]
        #fill jump vectors
        j[x...,i] = [x,y]
    end
    return η, c, j, Ω
end

function model_step!(config::Array{Int64,1},propensities::Array{Float64,1},jumps::Array{},rates::Array{};k =2, d=1, L=10)
    #select jump
    w     = Weights(propensities)
    jump  = sample(jumps, w)
    #execute jumps
    #println(η)
    #println(jump)
    config[jump[1]] += -1
    config[jump[2]] += 1
    #correct propensity
    for i in jump
        #find adjacent sites
        iplus  = (i % L)+1
        iminus = L - ((L+1-i)% L)
        # correct jumps out
        propensities[i]   = rates[config[i]+1,config[iplus]+1]
        propensities[i+L] = rates[config[i]+1,config[iminus]+1]
        # correct jumps in
        propensities[iplus+L]   = rates[config[iplus]+1,config[i]+1]
        propensities[iminus] = rates[config[iminus]+1,config[i]+1]
    end
    return config, propensities, jumps
end