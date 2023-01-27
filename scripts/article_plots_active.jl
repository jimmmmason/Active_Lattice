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



### fig 1
#convergence of λₙ
χ = 1.0
    Dθ = 4.0
    Dx = 1.
    xtic = 0.4:0.2:1
    ytic = 0:10:100 
    axlim = [0.4, 1., 0., 40.]
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    stabdata = Dict()
    for k in [2,4,6,40]
        X=[]
        Y=[]
        for x in xs
            try
                local f
                f(y) = lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k)
                Pe = find_zero(f, (0.,  100.))
                push!(Y,Pe)
                push!(X,x)
            catch
            end
        end
        if k ==40
            ax.plot(X,Y,color = "black", label = "n = $(k)")
        else
            ax.plot(X,Y,color = "black",linestyle="dashed",label = "n = $(k)")
        end
    end
    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    ax.set_title(L"\lambda_n = 0",fontsize=15)
    ax.set_xlabel(L"\rho",fontsize=15)
    ax.set_ylabel(L"\mathrm{Pe}", fontsize=15)
    ax.legend(loc = "upper right", fontsize=15)
    #ax.set_title("ℓ = $(1/sqrt(Dθ))")
display(fig)
ax.scatter(0.55,14)
name = "figure_1_λ_convergence"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#
###



### fig 2
#rand pde
#
###



### fig 3
#pde dist from uniform
#
###



### fig 4
#stability with pde
#
###



### fig 5
#sim examples
#
###



### fig 6
#sim vs pde hist
params = []
        pert = "rand"
        T  = 4.0
        χ = 1.0
        L = 128
        Δt = 0.01
        #using Roots
        #f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        #root = find_zero(f, (0.6,  0.8))
for ρ in [0.7], Pe in [10., 12.], Dθ in [4.]
        name = "article_sim_safe"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = L, Δt = Δt
                )
                #
                push!(params,param)
end
#ß
fig, axs = plt.subplots(1, 2, figsize=(12,5))
r = 5
t_end = 4.0
start_time = 3.0
for i in 1:3
    param = params[i]
    t_saves , η_saves = load_etas(param, t_end; dump_interval = 0.01, start_time = start_time);
    time_density_hist(fig, axs[i], param, t_saves, η_saves; r = r, bins = (2*r+1)^2 )
end
display(fig)
#
###