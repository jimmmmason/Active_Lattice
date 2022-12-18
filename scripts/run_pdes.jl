cd("/home/jm2386/Active_Lattice/")
using DrWatson, Roots
@quickactivate "Active_Lattice"
###
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
###
#varying parameters
params = []
name = "high_density_stability_v101"
        Pes = collect(4.:0.1:6.)
        T  = 1.0
        max_runs = 5
for ρa in collect(0.45:0.05:0.95), Dθ in [100., 10000.]
        local param
        λ_step = (0.1)*sqrt(Dθ)
        λs = Pes*sqrt(Dθ)
        λmax = (11.)*sqrt(Dθ)
        #
        param = pde_param(; name = name, ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = λ_step, λmax = λmax, λs=λs)
        #
        push!(params,param)
end
params = []
pert = "n=1"
T  = 0.003
f(x) = ap_lin_stab_line(x,0.5; Dx =1. ,Pe = 350., Dθ =100.,k=50 )
root = find_zero(f, (0.22,  0.23))
for ρa in [root], ρp in [0.5], Pe in [350.], Dθ in [100.], δ in [1e-5],k in [10,20,30,40]
        name = "stability_2d_actpass_k=$(k)"
        local param
                #
                param = pde_param_k(; name = name, 
                        ρa = ρa, Pe = Pe, ρp = ρp, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = 100, Nθ = 50, 
                        save_interval = 1e-5, max_steps = 1e8,
                        pert = pert, δ = δ, k=k
                )
                #
                push!(params,param)
end
#run pdes
@everywhere include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
@everywhere include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
###
length(params)/nworkers()
pmap(perturb_pde_run, params; distributed = true, batch_size=1, on_error=nothing,)
#