cd("/home/jm2386/Active_Lattice/")
using DrWatson
using Roots
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
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
χ = 0.75
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
    param = params[1]
#compute stability
for param in params
    stabdata = find_stab_data(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 1.0, stab_type = "full")
end
#stab plot
χ = 0.75
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    stabdata = wload(filename)
    plot_stab_frac(fig, ax, stabdata; xs=xs, param = param, χ = χ)
display(fig)
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/article_stabiliy/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/article_stabiliy/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#imag plot
χ = 0.75
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    plot_imaginary_frac(fig, ax; param = param, χ = χ)
display(fig)
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/article_imaginary/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/article_imaginary/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
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
#travelling wave plots
#
###


###
#sim pde comparison plots
#
###