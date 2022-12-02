cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
name = "linear_stability"
Pes = collect(50.:10.:100.)
ρs  = collect(0.8:0.05:1.)
Dθ = 100.
δ = 1e-4
δt = 1e-7
T  = 0.001
Nx = 100
Nθ = 25

pert_dists = []
for ρa in ρs, Pe in Pes
    local param, density, y1, y0
    param = pde_param(; 
        name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = 0.0001, δ= δ
    )
    density = initialize_density(param)
    #perturb_pde!(param,density; pert = "rand",δ = δ);
    perturb_pde!(param,density; pert = "n=1",δ = δ);

    @unpack fa, fp, t = density
    y0 = dist_from_unif(param, fa, fp)
    pde_step!(param,density)
    @unpack fa, fp, t = density
    y1 = dist_from_unif(param, fa, fp)

    Δy = y1 - y0

    push!(pert_dists,Δy)
end
pert_dists = reshape(pert_dists, (6,5))
fig, ax = PyPlot.subplots(figsize =(10, 10))
    colmap = PyPlot.plt.cm.viridis_r
    norm1 = matplotlib.colors.Normalize(vmin=0.01*minimum(pert_dists), vmax= 0.01*maximum(pert_dists) );
    ax.contourf(ρs, Pes, pert_dists; levels = 25, norm = norm1, cmap = colmap )#
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
display(fig)

pert_dists_save = deepcopy(pert_dists)


pert_dists = []
for ρa in ρs, Pe in Pes
    local param, density, y1, y0
    param = pde_param(; 
        name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = 50, Nθ = 20, save_interval = 0.0001, δ= δ
    )
    density = initialize_density(param)
    #perturb_pde!(param,density; pert = "rand",δ = δ);
    perturb_pde!(param,density; pert = "n=1",δ = δ);

    @unpack fa, fp, t = density
    y0 = dist_from_unif(param, fa, fp)
    pde_step_sym!(param,density)
    @unpack fa, fp, t = density
    y1 = dist_from_unif(param, fa, fp)

    Δy = y1 - y0

    push!(pert_dists,Δy)
end
pert_dists = reshape(pert_dists, (20,20))





