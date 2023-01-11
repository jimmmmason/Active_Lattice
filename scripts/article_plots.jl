cd("/home/jm2386/Active_Lattice/")
using DrWatson
using Roots
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#
@everywhere include("/home/jm2386/Active_Lattice/src/article_src.jl")
###
PyPlot.close("all")
###


###
#generic simulation plots
#
###


###
#generic 2d plots
#
###


###
#stability plots
#parameters
χ = 0.5
params = []
    pert = "n=1"
    T  = 1.0
    δ = 0.01
    Nx= 50
    Nθ= 20
    Dθ= 1.
    name = "article_stability_1d_δ=$(δ)"
    ρs = 0.01:0.01:0.99
    Pes = 5.:5.:100.
    stab_type = "full"
for ρ in [0.5], χ in [0.25,0.5,0.75], Pe in [100.]
        local param
        #
        name = "article_stability_1d_δ=$(δ)"
        param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
                )
        #
        push!(params,param)
    end
    param = params[3]
#compute stability
for param in params
    stabdata = find_stab_data(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 1.0, stab_type = "full")
end
#stab plot
χ = 0.5
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    stabdata = wload(filename)
    plot_stab_frac(fig, ax, stabdata; xs=xs, ρs = ρs, param = param, χ = χ)
display(fig)
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/article_stabiliy/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/article_stabiliy/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#imag plot
χ = 0.75
Dθ= 4.0
    ρ = 0.5
    Pe = 100.
    pert = "n=1"
    T  = 1.0
    δ = 0.01
    Nx= 50
    Nθ= 20
    name = "article_stability_1d_δ=$(δ)"
    ρs = 0.01:0.01:0.99
    Pes = 5.:5.:100.
    stab_type = "full"
    param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
    )
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    plot_imaginary_frac(fig, ax; param = param, χ = χ)
display(fig)
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/article_imaginary/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/article_imaginary/$(name)/χ=$(χ)_Dθ=$(Dθ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###


###
#dist from uniform plots
#
χ = 0.75
ρ = 0.55
Pe = 10.
params = []
    pert = "n=1"
    T  = 1.0
    δ = 0.01
    Nx= 50
    Nθ= 20
    Dθ= 100.
    ρs = 0.05:0.05:0.95
    Pes = 5.:5.:100.
    stab_type = "full"
    
    param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 0.01, max_steps = 1e7,
                        pert = pert, k =40, δ = δ,
                )
#
t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.01)
    dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    ax.plot(t_saves,dist_saves)
    ax.set_xlabel("t")
    ax.set_ylabel("‖ρ-ϕ‖₂")
    ax.set_title("ℓ = $(sqrt(Dθ)) ρ = $(ρ) Pe = $(Pe) χ=$(χ)")
display(fig)
###



###
#travelling wave plots 1d
#
T  = 0.5
pert = "n=1"
    χ = 0.25
    using Roots
    f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    root = find_zero(f, (0.6,  0.8))
    ρ = -0.04+root
    Pe = 20. 
    Dθ = 100.
    δ  = 1e-6
    Nx = 128
    Nθ = 64
    name = "article_wave_1d_δ=$(δ)"
param = pde_param_fraction(; name = name, 
                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                save_interval = 1e-4, max_steps = 1e8,
                pert = pert, δ = δ, k = 40
)
#
frames = 100.
save_interval = 0.1*T/frames
t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
fig, axs = plt.subplots(3, 2, figsize=(12,12))
i =27
test_vid_phase_pde_plot_1d(fig, axs, param, t_saves, fa_saves, fp_saves, i)
display(fig)


minimum(ρa_saves[270])-ρa
minimum((minimum.(ρa_saves).-ρa))

minimum((minimum(ρa_saves[1]).-ρa))

spatial_fourier_mode2(ρ; Nx = Nx)


spatial_fourier_mode2s(ρa_saves,ρa_saves; Nx = Nx)



mean(ρp_saves[1])-ρp

maximum(maximum(ρa_saves))-ρa
minimum(minimum(ρa_saves))-ρa


animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves; frames = 99)
###


###
#sim pde comparison plot
#sim hist
T  = 2.0
        χ  = 0.75
        ρ  = 0.8
        Pe = 40.
        Dθ = 100.
        name = "article_sim"
param = sim_param_fraction(; name = name, 
        ρ = ρ, Pe = Pe, χ = χ, T = T, 
        Dθ = Dθ, L = 128, Δt = 0.001
)
fig, ax = PyPlot.subplots(figsize =(10, 10))
t_saves , η_saves = load_etas(param, t_end; dump_interval = 0.1, start_time = t_start);
r = 5
time_density_hist(fig, ax, param, t_saves, η_saves; r = r, bins = (2*r+1)^2 )
#pde hist 
pert = "rand"
        T  = 2.0
        χ  = 0.75
        ρ  = 0.8
        Pe = 40. 
        Dθ = 100. 
        δ  = 1e-2
        Nx = 50
        Nθ = 20
        name = "article_rand_2d_δ=$(δ)"
param = pde_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                        save_interval = 1e-4, max_steps = 1e8,
                        pert = pert, δ = δ, k = 40
)
t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = 0.01, start_time = T-0.02)
pde_density_hist(fig, ax, param, fa_saves[1], fp_saves[1]; bins = 1000)

