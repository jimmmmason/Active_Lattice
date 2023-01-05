cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
###
#=

Stability parameters

=#
###
name = "stability_1d_actpass_2"
T  = 1.0
Dθ = 1.
ρp = 0.3
ρs  = collect(0.4:0.05:0.95)
Pes = collect(5.:5.:100.)
#Pes = collect(0.:0.5:10.)
    Nx = 50
    Nθ = 20
    λs = Pes*sqrt(Dθ)
###
###
name = "stability_1d_actpass_2"
T  = 1.0
Dθ = 100.
ρp = 0.4
ρs  = collect(0.4:0.05:1.0)
Pes = collect(0.:5.0:100.)
    Nx = 50
    Nθ = 20
    λs = Pes*sqrt(Dθ)
###
###
name = "stability_2d_actpass_2"
T  = 2.0
Dθ = 100.
ρp = 0.4
ρs  = collect(ρp:0.1:1.0)
Pes = collect(0.:10.0:100.)
    Nx = 50
    Nθ = 20
    λs = Pes*sqrt(Dθ)
###
#=

Get stability data from saves

=#
###
stabdata = refresh_stab_data_1d(; ρs = ρs, Dθ = Dθ, Nx = Nx, Nθ = Nθ, λs = λs, name = name, ρp = ρp)
###
#=

Load previously saved stability

=#
###
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ).jld2"
stabdata = wload(filename)
###
###
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_ρp=$(ρp).jld2"
stabdata = wload(filename)

for ρ in ρs
    try
        stable = stabdata["ρ = $(ρ)"]["stable"]
        push!(stable,0.)
        stabdata["ρ = $(ρp)"]["stable"] = stable
    catch
        if ρp>0.4
            stable = []
            push!(stable,0.)
            stabdata["ρ = $(ρp)"]["stable"] = stable
        end
    end
end
stable = []
for λ in λs
    push!(stable,λ)
end
if ρp>0.4
    stabdata["ρ = $(ρp)"]["stable"] = stable
end
stabdata["ρ = $(1.)"]["stable"] = stable

#stabdata = Dict{String,Any}()

###
#=

Plot stability

=#
###
fig, ax = PyPlot.subplots(figsize =(10, 10))
xs = append!(collect(0.01:0.01:0.9),collect(0.901:0.001:0.999))
ys = collect(0:1:maximum(Pes))
plot_stab(fig, ax, stabdata; ρs = ρs ,xs = xs, ys =ys, xtic = (minimum(ρs)-ρp):0.1:(maximum(ρs)-ρp), ytic = 0:10:maximum(Pes), axlim = [minimum(ρs)-ρp, maximum(ρs)-ρp, minimum(Pes), maximum(Pes)], Dθ = Dθ,ρp=ρp)
display(fig)



#Dθ = (40.)^2
Dx = 1.
fig, ax = PyPlot.subplots(figsize =(10, 10))
xs = append!(collect(0.01:0.001:0.9),collect(0.901:0.001:0.999))
ys = collect(0:0.1:maximum(Pes))
xtic = (minimum(ρs)-ρp):0.1:(maximum(ρs)-ρp)
ytic = 0:10:maximum(Pes)
axlim = [minimum(ρs)-ρp, maximum(ρs)-ρp, minimum(Pes), maximum(Pes)]
ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    ax.set_xlabel("ρa")
    ax.set_ylabel("Pe")
    ax.set_title("Dθ= $(Dθ), ρᵖ = $(ρp)")
    x = xs
    y = ys
    n = length(x)
    m = length(y)
    z = zeros(m,n)
        for i in 1:m, j in 1:n
            z[i,j] = ap_lin_stab_line(x[j]-ρp,ρp; Dx =Dx ,Pe = y[i], Dθ = Dθ,k =20)
        end
ax.contour(x.+(-ρp),y,z; levels = [0])

n = length(x)
    m = length(y)
    z = zeros(m,n);
    Dx =1 
        for i in 1:m, j in 1:n
            z[i,j] = abs(ap_lin_stab_imaginary(x[j]-ρp,ρp; Dx =Dx ,Pe = y[i], Dθ = Dθ, k =20))
        end

    colmap = PyPlot.plt.cm.viridis_r
        norm1 = matplotlib.colors.Normalize(vmin=0., vmax= maximum(z) );
        ax.contourf(x.-ρp, y, z; levels = 200, norm = norm1, cmap = colmap )
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
#fig
display(fig)
###
#save
pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/stabiliy_imaginary/$(name)";
mkpath(pathname)
filename = "/store/DAMTP/jm2386/Active_Lattice/plots/stabiliy_imaginary/$(name)/Nx=$(Nx)_Nθ=$(Nθ).pdf";
PyPlot.savefig(filename,dpi = 100, format = "pdf")
###

ap_lin_stab_line(0.225,0.5; Dx =1. ,Pe = 377.037855, Dθ =100.,k=50 )

ap_lin_stab_imaginary(0.225,0.5; Dx =1. ,Pe = 377.1, Dθ =100.,k=50  )

1/40

x = xs
    y = ys
    n = length(x)
    m = length(y)
    z = zeros(m,n);
    Dx =1 
        for i in 1:m, j in 1:n
            z[i,j] = abs(ap_lin_stab_imaginary(x[j]-ρp,ρp; Dx =Dx ,Pe = y[i], Dθ = Dθ))
        end

    colmap = PyPlot.plt.cm.viridis_r
        norm1 = matplotlib.colors.Normalize(vmin=0., vmax= maximum(z) );
        ax.contourf(x.-ρp, y, z; levels = 200, norm = norm1, cmap = colmap )
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
#fig

    display(fig)

    ax.contour(x.+(-ρp),y,z; levels = [0.0001])



ap_lin_stab_imaginary(0.22,0.5; Dx =1. ,Pe = 500., Dθ = 10000.,k=50 )

using Roots

