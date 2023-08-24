cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#

### run all 1d pde simulations
## params for figure 3+4
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
        k = 20
        save_interval = 0.01
        χ = 1.0
        Dx = 1. 
        Dθ = 4.0
        Nx = 64
        Nθ = 32
name = "article_pde_run_1d_δ=$(δ)_l=$(1/sqrt(Dθ))"
#
Pes = 1:1:100
ρs = collect(0.4:0.02:0.99)
ρ = 0.7
for ρ in ρs, χ in [χ], Dθ in [Dθ]
        try
                critical_Pe = crit_Pe(ρ,χ; Dx = Dx ,Pe_max = 100., Dθ = Dθ, k=k)
        catch
                critical_Pe = -100.
        end
        for Pe in [Pe for Pe in Pes if abs(Pe - critical_Pe) ≤ 2]
                local param
                #
                param = pde_param_fraction(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                                save_interval = save_interval, max_steps = 1e7,
                                pert = pert, k =k, δ = δ,
                        )
                #
                push!(params,param)
        end
end
# 
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
#
pmap(perturb_pde_run_1d, params; distributed = true, batch_size=1, on_error=nothing,)
###
# #make videos
# pmap(pmap_make_phase_video_1d , params; distributed = true, batch_size=1, on_error=nothing,)
#
# compute stability 
param = params[1]
stab_type = "full"
stabdata = find_stab_data(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
stabdata = fillout_stab_data_horizontal(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
stabdata = fillout_stab_data(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
#


### run all 2d pde simulations
## params for figure 3 & 7
params = []
        pert = "rand"
        T  = 4.0
        δ  = 1e-2
        k = 20
        save_interval = 0.01
        χ = 1.0
        Dx = 1. 
        Dθ = 4.0
        Nx = 64
        Nθ = 32
name = "article_pde_run_2d_δ=$(δ)_l=$(1/sqrt(Dθ))"
#
T = 4.0
Pe = 30.
ρ = 0.5
param = pde_param_fraction(; name = name, 
        ρ = ρ, Pe = Pe, χ = χ, T = T, 
        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
        save_interval = save_interval, max_steps = 1e7,
        pert = pert, k =k, δ = δ,
)
push!(params,param)
#
#
T = 2.0
Pe = 12.
ρ = 0.7
param = pde_param_fraction(; name = name, 
        ρ = ρ, Pe = Pe, χ = χ, T = T, 
        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
        save_interval = save_interval, max_steps = 1e7,
        pert = pert, k =k, δ = δ,
)
push!(params,param)
# 
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
#
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
###
# #make videos
# pmap(pmap_make_phase_video_1d , params; distributed = true, batch_size=1, on_error=nothing,)
#
# compute stability 
param = params[1]
stab_type = "full"
stabdata = find_stab_data(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
stabdata = fillout_stab_data_horizontal(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
stabdata = fillout_stab_data(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
#

