cd("/home/jm2386/Active_Lattice/")
using DrWatson
using Roots
using LaTeXStrings
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
    xtic = 0.5:0.1:1
    ytic = 0:10:100 
    axlim = [0.45, 1., 0., 40.]
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
rc("text", usetex=true)
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))
    stabdata = Dict()
    for k in [2,4,6,40]
        X=[]
        Y=[]
        for x in xs
            try
                local f
                f(y) = a_lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k)
                Pe = find_zero(f, (0.,  100.))
                push!(Y,Pe)
                push!(X,x)
            catch
            end
        end
        if k ==40
            ax.plot(X,Y, label = latexstring("\$ n= $(k) \$")) #,color = "black"
        else
            ax.plot(X,Y,linestyle="dashed",label = latexstring("\$ n = $(k) \$"))
        end
    end
    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    ax.set_title(L"\lambda_n = 0",fontsize=20)
    ax.set_xlabel(L"\rho",fontsize=20)
    ax.set_ylabel(L"\mathrm{Pe}", fontsize=20)
    ax.legend(loc = "upper right", fontsize=20)
    #ax.set_title("ℓ = $(1/sqrt(Dθ))")
display(fig)
name = "figure_1_λ_convergence"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/lambdanplot.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#
###



### fig 2
#rand pde
ρ  = 0.7
Pe = 10. #12.
Dθ = 4.
        pert = "rand"
        T  = 4.0
        χ = 1.0
        L = 128
        Δt = 0.001 # need to run sims for this
        pert = "n=1"
        T  = 1.0
        δ  = 1e-2
name = "article_rand_2d_δ=$(δ)" #need to change name
param = pde_param_fraction(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                                save_interval = 0.01, max_steps = 1e7,
                                pert = pert, k =40, δ = δ,
                        )
                #
#
t_saves, fa_saves, fp_saves = load_pdes(param,2.0; save_interval = 1.0, start_time = 0.0)
#
rc("text", usetex=true)
fig, axs = plt.subplots(2, 3, figsize=(15,10), constrained_layout = true)
for i in 1:3
    fa = fa_saves[i]
    fp = fp_saves[i]
    t = t_saves[i]
    density = Dict{String,Any}()
    @pack! density = fa, fp, t
    plot_pde_mass(fig,axs[2*i-1],param,density; scale = "relative", cbar = false)
    plot_pde_mag( fig,axs[2*i],param,density; scale = [0., 0.6], cbar = false)
    if i ==1
        axs[2*i].yaxis.set_ticks(0.2:0.2:1.0)
        axs[2*i-1].yaxis.set_ticks(0.2:0.2:1.0)
    else
        axs[2*i-1].yaxis.set_ticks([])
        axs[2*i].yaxis.set_ticks([])
    end
    axs[2*i-1].set_title( latexstring("\$ t = $(round(t; digits = 2)) \$"),fontsize=20)
    axs[2*i-1].xaxis.set_ticks([])
    axs[2*i].xaxis.set_ticks(0.:0.2:1.0)
    axs[2*i-1].yaxis.set_tick_params(labelsize=15)
    axs[2*i].yaxis.set_tick_params(labelsize=15)
    axs[2*i-1].xaxis.set_tick_params(labelsize=15)
    axs[2*i].xaxis.set_tick_params(labelsize=15)
end
colmap = PyPlot.plt.cm.viridis_r
norm1 = matplotlib.colors.Normalize(vmin=0.0, vmax= 1.0)
norm2 = matplotlib.colors.Normalize(vmin=0.0, vmax= 0.8)
fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = axs[1:2:5], shrink=0.9)
fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm2, cmap = colmap), ax = axs[2:2:6], shrink=0.9)
title = latexstring("\$ \\phi = $(ρ), \\mathrm{Pe} = $(Pe),  \\ell = $(round(1/sqrt(Dθ); digits = 2)) \$")
fig.suptitle(title,fontsize=20)
display(fig)
name = "figure_2_PDE_timeseries"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/pdetimeseries_Pe=$(Pe)_.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###



### fig 3
#pde dist from uniform
#
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
        name = "article_stability_10_1d_δ=$(δ)"
for ρ in [0.7], χ in [1.0], Pe in [6., 8., 10., 12.], Dθ in [4.]
    local param
    #
    param = pde_param_fraction(; name = name, 
                    ρ = ρ, Pe = Pe, χ = χ, T = T, 
                    Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                    save_interval = 0.01, max_steps = 1e7,
                    pert = pert, k =40, δ = δ,
            )
    #
    push!(params,param)
end
#
#
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
        name = "article_stability_10_1d_δ=$(δ)"
for ρ in [0.47], χ in [1.0], Pe in [ 28., 30., 32.,34.], Dθ in [4.]
    local param
    #
    param = pde_param_fraction(; name = name, 
                    ρ = ρ, Pe = Pe, χ = χ, T = T, 
                    Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                    save_interval = 0.01, max_steps = 1e7,
                    pert = pert, k =40, δ = δ,
            )
    #
    push!(params,param)
end
#
params = []
        pert = "n=1"
        T  = 4.0
        δ  = 1e-4
        name = "article_stability_10_1d_δ=$(δ)"
for ρ in [0.99], χ in [1.0], Pe in [8., 10., 12., 14.], Dθ in [4.]
    local param
    #
    param = pde_param_fraction(; name = name, 
                    ρ = ρ, Pe = Pe, χ = χ, T = T, 
                    Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                    save_interval = 0.01, max_steps = 1e7,
                    pert = pert, k =40, δ = δ,
            )
    #
    push!(params,param)
end
#
rc("text", usetex=true)
fig, ax = plt.subplots(1, 1, figsize=(10,10))
#ax.plot(fa_saves[1])
#display(fig)
#i=1
for i in 1:4
    param = params[i]
    @unpack ρ, Pe, χ, Dθ = param
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = 0.01)
    dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)
    label = latexstring("\$ \\mathrm{Pe} = $(Pe) \$")
    ax.plot(t_saves,dist_saves, label = label)
end
ax.axis([0., 4.0 ,0., .0004 ],fontsize=15)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.set_xlabel(L"t", fontsize=20)
    ax.set_ylabel(L"\Vert \rho - \phi \Vert_{H_2}", fontsize=20)
    title = latexstring(
        "\$ \\ell  =  $(1/sqrt(Dθ))\$, \$ \\phi =  $(ρ) \$"
    )
    ax.set_title(title,fontsize=20)
    ax.legend(loc= "upper left", fontsize=20)
    ax.plot(0:0.01:1,5*δ*ones(100),linestyle="dashed")
display(fig)
name = "figure_3_PDE_stab_timeseries"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/pdetimeseries.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###



### fig 4
#stability with pde
#
params = []
        pert = "n=1"
        T  = 1.0
        δ  = 1e-4
        name = "article_stability_10_1d_δ=$(δ)"
        Nx = 128
        Nθ = 64
for ρ in [0.7], χ in [1.0], Pe in [6.], Dθ in [4.]
    local param
    #
    param = pde_param_fraction(; name = name, 
                    ρ = ρ, Pe = Pe, χ = χ, T = T, 
                    Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                    save_interval = 0.01, max_steps = 1e7,
                    pert = pert, k =40, δ = δ,
            )
    #
    push!(params,param)
end
param = params[1]
#
ρs = collect(0.45:0.01:1.0)
Pes = 2.:2.:40.
stab_type = "full"
for param in params
    stabdata = find_stab_data(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = "full")
end
#
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
stabdata = wload(filename)
stabdata = find_stab_data(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = "full")
#
χ = 1.0
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
fig, ax = PyPlot.subplots(figsize =(10, 10))
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    stabdata = wload(filename)
    plot_stab_frac(fig, ax, stabdata; xs=xs, ρs = ρs, xtic = 0.4:0.1:1.0,  axlim = [minimum(ρs), maximum(ρs), minimum(Pes), maximum(Pes)], param = param, χ = χ)
    ax.set_xlabel(L"\rho")
    ax.set_ylabel("Pe")
    #ax.set_title("ℓ = $(1/sqrt(Dθ)), χ = $(χ)")
display(fig)
name = "figure_4_sim_stab_pde"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###


### fig 5
#sim examples
params = []
ρ  = 0.7
Dθ = 4.
        pert = "rand"
        T  = 4.0
        χ = 1.0
        L = 32
        Δt = 0.001 # need to run sims for this
        pert = "n=1"
        T  = 1.0
        δ  = 1e-2
        name = "article_sim_safe"
for Pe in [8.,10.,12.]
        param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = L, Δt = Δt
                )
                push!(params,param)
end
#
rc("text", usetex=true)
fig, axs = plt.subplots(3, 3, figsize=(15,15))
for j in 1:3
    param = params[j]
    for k in 1:3
        i = 3*k + j -3
        #
        t_end = 4.1
        start_time = 0.0
        t_saves , η_saves = load_etas(param, t_end; dump_interval = 0.5, start_time = start_time);
        #   
        η = η_saves[k]
        t = t_saves[k]
        density = Dict{String,Any}()
        plot_eta(fig, axs[i], param, t, η; title = false, r = 3)
        if k ==1
            axs[i].yaxis.set_ticks(0.0:0.2:1.0)
        else
            axs[i].yaxis.set_ticks([])
        end
        if j == 3
            axs[i].xaxis.set_ticks(0.0:0.2:1.0)
        else
            axs[i].xaxis.set_ticks([])
        end
        if j ==1
            axs[i].set_title( latexstring("\$ t = $(round(t; digits = 2)) \$"),fontsize=20)
        end
        if k ==1 
            @unpack Pe =param
            axs[i].set_ylabel(latexstring("\$ \\mathrm{Pe} = $(Pe) \$"),fontsize=20)
        end
        axs[i].yaxis.set_tick_params(labelsize=15)
        axs[i].xaxis.set_tick_params(labelsize=15)
    end
end
@unpack ρ, Dθ =params[1]
title = latexstring("\$ \\phi = $(ρ), \\ell = $(round(1/sqrt(Dθ); digits = 2)) \$")
fig.suptitle(title,fontsize=20)
display(fig)
name = "figure_5_sim_timeseries"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/simtimeseries.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###



### fig 5.1
#sim examples density plot
params = []
        pert = "rand"
        T  = 4.0
        χ = 1.0
        #L = 128
        Δt = 0.001 # need to run sims for this
        pert = "n=1"
        T  = 1.0
        δ  = 1e-2
        Nx = 50
        Nθ = 20
for ρ in [0.7], Pe in [12.], Dθ in [4.], L in [128]
        name = "article_sim_safe"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = L, Δt = Δt
                )
                #
                push!(params,param)
                #
end
param = params[1]
#
t_end = 1.1
start_time = 0.0
    @unpack ρ, Pe, Dθ, L = param
    r = 4
    numbins = (2*r+1)^2
    t_saves , η_saves = load_etas_1(param, t_end; dump_interval = 0.5, start_time = start_time);
    for i in 1:3
        t = t_saves[i]
        η = η_saves[i]
        time_density_hist_time_save(param, t, η; r = r, bins = numbins)
        time_density_polarisation_save(param, t, η; r = r, bins = numbins)
    end
#
rc("text", usetex=true)
fig, axs = plt.subplots(2, 3, figsize=(15,10), constrained_layout = true)
for i in 1:3
    t = t_saves[i]
    plot_sim_mass(fig,axs[2*i-1],param,t; scale = "relative", cbar = false)
    plot_sim_mag( fig,axs[2*i],param,t; scale = [0., 0.6], cbar = false)
    if i ==1
        axs[2*i].yaxis.set_ticks(0.2:0.2:1.0)
        axs[2*i-1].yaxis.set_ticks(0.2:0.2:1.0)
    else
        axs[2*i-1].yaxis.set_ticks([])
        axs[2*i].yaxis.set_ticks([])
    end
    axs[2*i-1].set_title( latexstring("\$ t = $(round(t; digits = 2)) \$"),fontsize=20)
    axs[2*i-1].xaxis.set_ticks([])
    axs[2*i].xaxis.set_ticks(0.:0.2:1.0)
    axs[2*i-1].yaxis.set_tick_params(labelsize=15)
    axs[2*i].yaxis.set_tick_params(labelsize=15)
    axs[2*i-1].xaxis.set_tick_params(labelsize=15)
    axs[2*i].xaxis.set_tick_params(labelsize=15)
end
@unpack Pe, ρ, Dθ= param
    colmap = PyPlot.plt.cm.viridis_r
    norm1 = matplotlib.colors.Normalize(vmin=0.0, vmax= 1.0)
    norm2 = matplotlib.colors.Normalize(vmin=0.0, vmax= 0.8)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = axs[1:2:5], shrink=0.9)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm2, cmap = colmap), ax = axs[2:2:6], shrink=0.9)
    title = latexstring("\$ \\phi = $(0.7), \\mathrm{Pe} = $(Pe),  \\ell = $(round(1/sqrt(Dθ); digits = 2)) \$")
    fig.suptitle(title,fontsize=20)
display(fig)
name = "figure_5.5_sim_density_timeseries"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/simdentimeseries_Pe=$(Pe)_r=$(r)_.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###
###


### fig 6
#sim timeseries
params = []
        pert = "rand"
        T  = 2.0
        χ = 1.0
        L = 128
        Δt = 0.001
        #using Roots
        #f(x) = lin_stab_line_fraction(x,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
        #root = find_zero(f, (0.6,  0.8))
for ρ in [0.7], Pe in [6., 8., 10., 12.], Dθ in [4.]
        name = "article_sim_safe"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = L, Δt = Δt
                )
                #
                push!(params,param)
                local param
    #
end
#
rc("text", usetex=true)
fig, ax = plt.subplots(1, 1, figsize=(10,10))
t_end = 4.0
start_time = 0.0
for i in 1:4
    param = params[i]
    @unpack ρ, Pe, χ, Dθ, L = param
    t_saves , η_saves = load_etas_1(param, t_end; dump_interval = 0.01, start_time = start_time);
    invar_saves = translation_invariance.(η_saves;L = L)
    label = latexstring("\$ \\mathrm{Pe} = $(Pe) \$")
    ax.plot(t_saves,invar_saves, label = label)
end
ax.axis([0., 2. ,0., .3 ],fontsize=15)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.set_xlabel(L"t", fontsize=20)
    ax.set_ylabel(L"\Phi", fontsize=20)
    @unpack Dθ,ρ = params[1]
    title = latexstring(
        "\$ \\ell  =  $(1/sqrt(Dθ))\$, \$ \\phi =  $(ρ) \$"
    )
    ax.set_title(title,fontsize=20)
    ax.legend(loc= "center right", fontsize=20, bbox_to_anchor = [1., 0.35])
display(fig)
name = "figure_6_sim_stab_timeseries"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/simtimeseries.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#
###



### fig 7
#sim vs pde hist
params = []
params_pde = []
        pert = "rand"
        T  = 4.0
        χ = 1.0
        #L = 128
        Δt = 0.001 # need to run sims for this
        pert = "n=1"
        T  = 1.0
        δ  = 1e-2
        Nx = 128
        Nθ = 64
for ρ in [0.7], Pe in [8., 12.], Dθ in [4.], L in [32,64,128]
        name = "article_sim_safe"
        local param
                #
                param = sim_param_fraction(; name = name, 
                        ρ = ρ, Pe = Pe, χ = χ, T = T, 
                        Dθ = Dθ, L = L, Δt = Δt
                )
                #
                push!(params,param)
                #
                name = "article_rand_actually_2d_δ=$(δ)" #need to run pde in 2d 
                param = pde_param_fraction(; name = name, 
                                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                                save_interval = 0.01, max_steps = 1e8,
                                pert = pert, k =40, δ = δ,
                        )
                #
                push!(params_pde,param)
end
# Collect hist vectors
for i in 1:2
    param_pde = params_pde[3*i]
    t_saves, fa_saves, fp_saves = load_pdes_rand(param_pde,6.0; save_interval = 0.02, start_time = 2.0)
    pde_density_hist_av_save(param_pde, fa_saves, fp_saves; bins = 17^2, smoothing = false)
end
t_end = 2.0
start_time = 1.0
for j in 1:6
    param = params[j]
    @unpack ρ, Pe, Dθ, L = param
    r = Int64(L/16)
    numbins = (2*r+1)^2
    t_saves , η_saves = load_etas_1(param, t_end; dump_interval = 0.1, start_time = start_time);
    time_density_hist_save(param, t_saves, η_saves; r = r, bins = numbins)
end
#
rc("text", usetex=true)
fig, axs = plt.subplots(1, 2, figsize=(15,10))
t_end = 4.0
start_time = 3.5
for i in 1:2
    param_pde = params_pde[3*i]
    ax = axs[i]
    if i==1
        add_label = true
    else
        add_label = false
    end
    for j in (3*i-2):1:(3*i)
        param = params[j]
        @unpack ρ, Pe, Dθ, L = param
        name = latexstring(
            "\$ N = $(L) \$"
        )
        r = Int64(L/16)
        numbins = (2*r+1)^2
        time_density_hist_load(fig, axs[i], param; r = r, bins = numbins, label_name=name, add_label = add_label )
    end
    #
    name = latexstring(
        "\$ N = \\infty \$"
    )
    #pde_density_hist_av_load(fig, ax, param_pde; bins = 17^2, smoothing = false, label_name=name, add_label = add_label ) # bins = 17^2
    #
    ax.axis(fontsize=15)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.set_xlabel(L"\phi", fontsize=20)
    #ax.set_ylabel(L"\Phi", fontsize=20)
    @unpack Pe = param_pde
    title = latexstring(
        "\$ \\mathrm{Pe} = $(Pe) \$"
    )
    ax.set_title(title,fontsize=20)
    if i == 1
        ax.legend(loc = "upper left", fontsize=20)
        ax.axis([0., 1.0 ,0., 50. ],fontsize=15)
    else
        ax.axis([0., 1.0 ,0., 8. ],fontsize=15)
    end
end
    @unpack Dθ,ρ = params[1]
    title = latexstring(
    "\$ \\ell  =  $(1/sqrt(Dθ)),  \\phi =  $(ρ) \$"
    )
    fig.suptitle(title,fontsize=20)
display(fig)
name = "figure_7_histograms"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_paper/$(name)/hist.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#
###




#=
17/25
18/25

ρ_sim = []
for x₁ in 1:L, x₂ in 1:L
    local ρr 
    #ρr = 0.
    #=
    for e ∈ B
        y₁, y₂  = ([x₁, x₂] + e +[L-1,L-1]) .% L + [1,1]
        ρr += site_ρ(η[y₁, y₂,: ])/(2*r +1)^2
    end
    =#
    ρr = site_ρ(η[x₁, x₂,: ])
    push!(ρ_sim,ρr) 
end
sum(ρ_sim)/L^2

r = 2
ρ_sim_h = []
B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
for x₁ in 1:L, x₂ in 1:L
    local ρr 
    ρr = 0.
    for e ∈ B
        y₁, y₂  = ([x₁, x₂] + e +[L-1,L-1]) .% L + [1,1]
        ρr += site_ρ(η[y₁, y₂,: ])/(2*r +1)^2
    end
    push!(ρ_sim_h,ρr) 
end
sum(ρ_sim_h)/L^2
mode(ρ_sim_h)

fig, ax = plt.subplots(1, 1, figsize=(15,5))
r = 2
numbins = (2*r+1)^2
edges = collect((-1/(2*numbins)):(1/(numbins)):(1+(1/(2*numbins))))
ax.hist(ρ_sim_h; bins = edges, histtype = "step", density = true)
r = 10
numbins = (2*r+1)^2
edges = collect(0.01:(1/(numbins)):1.0)
ax.hist(ρ_sim_h; bins = edges, histtype = "step", density = true)
display(fig)

ρ = sum(fa; dims =3)
2*π*minimum(ρ)/64
2*π*maximum(ρ)/64
=#