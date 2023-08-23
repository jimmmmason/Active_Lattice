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


### λ_imag
# 
χ = 1.0
    Dθ = 100.0
    Dx = 1.
    ρ_start = 0.4
    ρ_end = 1.0
    Pe_end = 100
    Pe_start = 0
    norm2 = matplotlib.colors.Normalize(vmin=minimum(0), vmax= maximum(40) )
    xtic = 0.4:0.1:1
    ytic = Pe_start:((Pe_end-Pe_start)/10):Pe_end
xs = collect(0.4:0.001:0.999) # append!(collect(0.01:0.01:0.99),collect(0.99:0.0001:0.9999))
rc("text", usetex=true)
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))

    k= 10
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
    ax.plot(X,Y,color = "black") #label = latexstring("\$ n= $(k) \$")) #,color = "black"

    λ_imag = zeros(801,201)
    Δρ = (ρ_end-ρ_start)/200
    ΔPe = (Pe_end-Pe_start)/800

    for i in 0:1:200, j in 0:1:800
        ρ = Δρ*i+ρ_start
        Pe = ΔPe*j +Pe_start
        λ_imag[801-j,i+1] = abs.(lin_imaginary_fraction(ρ,χ; Dx =Dx ,Pe = Pe, Dθ =Dθ,k=k ))
    end
    colmap = PyPlot.plt.cm.viridis
    norm1 = matplotlib.colors.Normalize(vmin=minimum(λ_imag), vmax= maximum(λ_imag) )
    ax.matshow(λ_imag; norm = norm2,  cmap = colmap, extent = [ρ_start, ρ_end, Pe_start, Pe_end])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm2, cmap = colmap), ax = ax)# fraction = 0.0455)
ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    axlim = [ρ_start, ρ_end, Pe_start, Pe_end]
    ax.xaxis.set_tick_params(labelsize=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    ax.set_xlabel(L"\phi",fontsize=20)
    ax.set_ylabel(L"\mathrm{Pe}", fontsize=20)
    #ax.legend(loc = "upper left", fontsize=20)
    ax.set_aspect(0.25*Δρ/ΔPe)
    title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ) \$")
    ax.set_title(title,fontsize=20)
display(fig)
#save fig
name = "λ_imag_plot_χ=$(χ)_Dθ=$(Dθ)_Dx=$(Dx)"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)/lambdanplot.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#extract Pe
Pemin = minimum(Y)
i = argmin(Y)
ϕmin = X[i]
###



### pde video
params = []
        pert = "n=1"
        T  = 24.0
        δ  = 1e-4
        k = 20
        save_interval = 0.01
        χ = 0.1
        Dθ = 4.0
        ϕmin = 0.929
        Pemin = 41.1445070772628
        Pe = Pemin
name = "ap_wave_1d_δ=$(δ)_χ=$(χ)_l=$(1/sqrt(Dθ))"
param = pde_param_fraction(; name = name, 
                                ρ = ϕmin, Pe = Pe, χ = χ, T = T, 
                                Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                                save_interval = save_interval, max_steps = 1e7,
                                pert = pert, k =k, δ = δ,
)
for ρ in [ϕmin], χ in [χ], Pe in [Pemin, (Pemin +1.0)], Dθ in [Dθ]
                            local param
                            #
                            param = pde_param_fraction(; name = name, 
                                            ρ = ρ, Pe = Pe, χ = χ, T = T, 
                                            Dθ = Dθ, δt = 1e-5, Nx = 128, Nθ = 64, 
                                            save_interval = save_interval, max_steps = 1e7,
                                            pert = pert, k =k, δ = δ,
                                    )
                            #
                            push!(params,param)
end
# video frame
i=1
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    fig, axs = plt.subplots(2, 1, figsize=(8,8))
    ldg = vid_pde_plot_1d(fig, axs, param, t_saves, fa_saves, fp_saves, i+1)
display(fig)
axs[1].legend(loc = "upper right", fontsize = 20)
#fig.tight_layout()
name = "wave_frame_χ=$(χ)_Dθ=$(Dθ)_Dx=$(Dx)_ϕ=$(ϕmin)_Pe=$(Pe)"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)/lambdanplot.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf",  bbox_extra_artists=( ldg,))
# create video
make_phase_video_1d_plus(param; frames = 240)
###



### λ_check
# 
χ = 0.5
    Dθ = 100.0
    Dx = 1.
    ρ_start = 0.0
    ρ_end = 1.0
    Pe_end = 100
    Pe_start = 0
    norm2 = matplotlib.colors.Normalize(vmin=minimum(0), vmax= maximum(40) )
    xtic = 0.4:0.1:1
    ytic = Pe_start:((Pe_end-Pe_start)/10):Pe_end
xs = collect(0.0:0.001:0.999) # append!(collect(0.01:0.01:0.99),collect(0.99:0.0001:0.9999))
rc("text", usetex=true)
k= 10
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))

    X=[]
    Y=[]
    # using Roots
    # for x in xs
    #         try
    #             local f
    #             f(y) = lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k)
    #             Pe = find_zero(f, (0.,  100.))
    #             push!(Y,Pe)
    #             push!(X,x)
    #         catch
    #         end
    # end
    # ax.plot(X,Y,color = "black") #label = latexstring("\$ n= $(k) \$")) #,color = "black"

    using Roots
    for ω in [ [i*2*π,j*2*π] for j in 0:1:3 for i in 1:2]
        X=[]
        Y=[]
        for x in xs
                try
                    local f
                    f(y) = lin_stab_line_fraction_checker(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k, ω = ω )
                    Pe = find_zero(f, (0.,  500.))
                    push!(Y,Pe)
                    push!(X,x)
                catch
                end
        end
        lab = latexstring("\$ \\omega = $(ω) \$")
        ax.plot(X,Y, label = lab) #color = "red"
    end
    # for ω in [ [2*π, j*2*π] for j in 0:1:2 ]
    #     X=[]
    #     Y=[]
    #     for x in xs
    #             try
    #                 local f
    #                 f(y) = lin_stab_line_fraction_checker(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k, ω = ω )
    #                 Pe = find_zero(f, (0.,  500.))
    #                 push!(Y,Pe)
    #                 push!(X,x)
    #             catch
    #             end
    #     end
    #     lab = latexstring("\$ \\omega = $(ω) \$")
    #     ax.plot(X,Y, label = lab) #color = "red"
    # end
    #=
    λ_imag = zeros(201,201)
    Δρ = (ρ_end-ρ_start)/200
    ΔPe = (Pe_end-Pe_start)/200
    for i in 0:1:200, j in 0:1:200
        ρ = Δρ*i+ρ_start
        Pe = ΔPe*j +Pe_start
        λ_imag[201-j,i+1] = abs.(lin_imaginary_fraction(ρ,χ; Dx =Dx ,Pe = Pe, Dθ =Dθ,k=k ))
    end
    colmap = PyPlot.plt.cm.viridis
    norm1 = matplotlib.colors.Normalize(vmin=minimum(λ_imag), vmax= maximum(λ_imag) )
    ax.matshow(λ_imag; norm = norm2,  cmap = colmap, extent = [ρ_start, ρ_end, Pe_start, Pe_end])
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm2, cmap = colmap), ax = ax)# fraction = 0.0455)
    =#
ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    axlim = [ρ_start, ρ_end, Pe_start, Pe_end]
    ax.xaxis.set_tick_params(labelsize=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    ax.set_xlabel(L"\phi",fontsize=20)
    ax.set_ylabel(L"\mathrm{Pe}", fontsize=20)
    ax.legend(loc = "upper left", fontsize=20)
    #ax.set_aspect(Δρ/ΔPe)
    title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ) \$")
    ax.set_title(title,fontsize=20)
display(fig)
#save fig
name = "λ_imag_plot_χ=$(χ)_Dθ=$(Dθ)_Dx=$(Dx)"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)/lambdanplot.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
#extract Pe
Pemin = minimum(Y)
i = argmin(Y)
ϕmin = X[i]
###



###
params = []
        pert = "n=1"
        T  = 12.0
        δ  = 1e-2
        k = 20
        save_interval = 0.01
        χ = 0.5
        Dθ = 100.0
        Nx = 128
        Nθ = 8
name = "ap_wave_1d_δ=$(δ)_χ=$(χ)_l=$(1/sqrt(Dθ))"
#
ρ = 0.7
Pe = 5
param = pde_param_fraction(; name = name, 
                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                save_interval = save_interval, max_steps = 1e7,
                pert = pert, k =k, δ = δ,
        )
#
# plot phase diagram
ρs = collect(0.5:0.05:0.99)
Pes = 5:5:40
stab_type = "full"
stabdata = find_stab_data_ap(;stabdata = Dict{String,Any}(), ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 20.0, t_start = 10.0, stab_type = stab_type, save_interval = 0.1)
#
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
stabdata = wload(filename)
#
# stabdata = fillout_stab_data(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
# stabdata = fillout_stab_data_horizontal(;stabdata = stabdata, ρs = ρs, Pes = Pes,  param = param, save_on = true, t_end = 4.0, stab_type = stab_type, save_interval = 0.25)
# 
#
rc("text", usetex=true)
using Roots
xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
fig, ax = PyPlot.subplots(figsize =(10, 10))
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    stabdata = wload(filename)
    plot_stab_frac_ap(fig, ax, stabdata; xs=xs, ρs = ρs, xtic = 0.4:0.1:1.0,  axlim = [minimum(ρs), maximum(ρs), minimum(Pes), maximum(Pes)], param = param, χ = χ)
    xtic = 0.0:0.1:1
    ytic = 0:10:100
    axlim = [minimum(ρs), maximum(ρs), minimum(Pes), maximum(Pes)]
    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    #ax.set_title(L"\lambda_n = 0",fontsize=20)
    #ax.set_xlabel(L"\phi",fontsize=20)
    #ax.set_ylabel(L"\mathrm{Pe}", fontsize=20)
    ax.legend(loc = "upper right", fontsize=20)
    #ax.set_title("ℓ = $(1/sqrt(Dθ)), χ = $(χ)")
    #title = latexstring("\\ell = $(round(1/sqrt(Dθ); digits = 2)) \$")
    #fig.suptitle("",fontsize=20)
display(fig)
name = "ap_stab+pde_mix"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive/$(name)/χ=$(χ)_Nx=$(Nx)_Nθ=$(Nθ).pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf")
###
###
###
params = []
        pert = "n=1"
        T  = 12.0
        δ  = 1e-2
        k = 20
        save_interval = 0.01
        χ = 0.5
        Dθ = 100.0
        Nx = 128
        Nθ = 8
name = "ap_repeating_pert_1d_δ=$(δ)_χ=$(χ)_l=$(1/sqrt(Dθ))"
name = "ap_wave_repeatingpert_1d_δ=$(δ)_χ=$(χ)_l=$(1/sqrt(Dθ))"
#
#
ρ = 0.95
Pe = 12.5
T = 36
rc("text", usetex=false)
PyPlot.close("all")
param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = Pe, χ = χ, T = T, 
                            Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = 1e8,
                            pert = pert, k =k, δ = δ, frames = 200., pert_interval = 0.5
)

@unpack T, save_interval,frames = param
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)

    fig, ax = plt.subplots(1, 1, figsize=(10,10))

    skip = 10
    frames = length(t_saves)

    skip_dist_saves = time_dist_from_past_1d(param, fa_saves, fp_saves, skip)
    unif_dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)

    ax.plot(t_saves,unif_dist_saves, color = "red", label = L"dist from unif")
    ax.plot(t_saves[(skip+1):1:frames],skip_dist_saves, color = "black", label = L"time skip")
    ax.legend(loc = "center left", fontsize=20, bbox_to_anchor = (1,0.5) )

display(fig)
###
make_phase_video_1d_plus(param)
###

function stationary_type(params; 
    skips = [5,10,15], 
    frames = 200,
    reset = true,
    )
    param = params[1]
    @unpack T, save_interval,name, δ = param

    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability.jld2"

    homo_stable = []
    inhomo_stable = []
    double_unstable = []
    unclassified = []

    stabdata = Dict{String,Any}()
    if reset 
    else
        try 
            stabdata = wload(filename)
            @unpack homo_stable, inhomo_stable, double_unstable, unclassified = stabdata
        catch
        end
    end

    for param in params
        save_interval = T/frames
        t_saves, fa_saves, fp_saves = load_pdes_1d(param,T; save_interval = save_interval, start_time = T/4)

        skip_dist_max = 0.
        for skip in skips
            skip_dist_saves = time_dist_from_past_1d(param, fa_saves, fp_saves, skip)
            skip_dist_max = max(skip_dist_max,maximum(skip_dist_saves))
        end

        unif_dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)

        
        unif_dist_max = maximum(unif_dist_saves)

        if (unif_dist_max< 2*δ)&(param ∉ homo_stable)
            push!(homo_stable,param)
        elseif (skip_dist_max< 0.1)&(unif_dist_max > 4*δ)&(param ∉ inhomo_stable)
            push!(inhomo_stable,param)
        elseif (skip_dist_max > 0.14)&(unif_dist_max > 4*δ)&(param ∉ double_unstable)
            push!(double_unstable,param)
        elseif (param ∉ homo_stable)&(param ∉ inhomo_stable)&(param ∉ double_unstable)&(param ∉ unclassified)
            push!(unclassified,param)
        end
    end
    
    @pack! stabdata = homo_stable, inhomo_stable, double_unstable, unclassified
    wsave(filename,stabdata)

    return stabdata
end

function plot_stationary_data(param; 
    xs = collect(0.4:0.001:1.0), 
    axlim = [0.4,1.0, 0., 20.],
    xtic = 0.0:0.1:1.0,
    ytic = 0.:10.:100.,
    )
    @unpack T, save_interval,frames,name, δ, Dθ, χ = param

    chi = χ
    Dtheta = Dθ
    delta = δ
    
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability.jld2"
    stabdata = wload(filename)
   
    index = Dict{String,Any}("homo_stable"=>1, "inhomo_stable"=>2, "double_unstable"=> 3, "unclassified" => 4 )
    fmts = ["s","o","*","v" ]
    colors = ["blue","green","red","purple"]
    labels = [L"\mathrm{stable-unstable}",L"\mathrm{unstable-stable}",L"\mathrm{unstable-unstable}",L"\mathrm{unclassified}"]
    fig, ax = plt.subplots(1, 1, figsize=(10,10))

    for (element, params) in stabdata
        
        x_points = []
        y_points = []

        for param_k in params
                @unpack Dθ, χ, δ, Pe, ρ  = param_k
                if (Dθ, χ, δ) == (Dtheta, chi, delta)
                    push!(x_points,ρ)
                    push!(y_points, Pe)
                end
        end

        INDEX = index[element]

        ax.errorbar(x_points,y_points, 
                fmt= fmts[INDEX], 
                color = colors[INDEX],
                alpha=0.8,
                label = labels[INDEX],
        )

    end

    X=[]
    Y=[]
    for x in xs
        try
            local f
            f(y) = lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ)
            Pe = find_zero(f, (0.,  100.))
            push!(Y,Pe)
            push!(X,x)
        catch
        end
    end
    ax.plot(X,Y,color = "black", label = L"\mathrm{Linear}")

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    ax.set_xlabel(L"\phi",fontsize=20)
    ax.set_ylabel(L"\mathrm{Pe}", fontsize=20)
    ax.legend(loc = "upper left", fontsize=20)
    #ax.set_aspect(1.)
    title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ) \$")
    ax.set_title(title,fontsize=20)

    return  fig, ax
end
     
stationary_type(params;     
skips = [5,10,15], 
frames = 200,
reset = true,
)
param = params[1]
fig, ax = plot_stationary_data(param)
display(fig)


# @unpack T, save_interval, max_steps, pert, δ, pert_interval = param
    # density = initialize_density_1d(param)
    # nudge_pde_1d!(param, density; δ = 0.01)

    # @unpack fa, fp = density
    # ρ = fp + sum(fa; dims = 2)[:,1]*2*π/Nθ
    # fig, ax = plt.subplots(1, 1, figsize=(10,10)) 
    # ax.plot(ρ)
    # fig
    # repeating_perturb_pde_run_1d(param)
# e
