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
    Dθ = 4.0
    Dx = 1.
    ρ_start = 0.4
    ρ_end = 1.0
    Pe_end = 50
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
    ax.set_aspect(Δρ/ΔPe)
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
make_phase_video_1d(param; frames = 2400)
###



### λ_check
# 
χ = 0.8
    Dθ = 4.0
    Dx = 1.
    ρ_start = 0.4
    ρ_end = 1.0
    Pe_end = 50
    Pe_start = 0
    norm2 = matplotlib.colors.Normalize(vmin=minimum(0), vmax= maximum(40) )
    xtic = 0.4:0.1:1
    ytic = Pe_start:((Pe_end-Pe_start)/10):Pe_end
xs = collect(0.4:0.001:0.999) # append!(collect(0.01:0.01:0.99),collect(0.99:0.0001:0.9999))
rc("text", usetex=true)
fig, ax = fig, ax = PyPlot.subplots(figsize =(10, 10))

    k= 5
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


    for ω in [ [j*2*π, 2*π] for j in 0:1:2 ]
        X=[]
        Y=[]
        for x in xs
                try
                    local f
                    f(y) = lin_stab_line_fraction_checker(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k, ω = ω )
                    Pe = find_zero(f, (0.,  100.))
                    push!(Y,Pe)
                    push!(X,x)
                catch
                end
        end
        lab = latexstring("\$ \\omega = $(ω) \$")
        ax.plot(X,Y, label = lab) #color = "red"
    end
    for ω in [ [2*π, j*2*π] for j in 0:1:2 ]
        X=[]
        Y=[]
        for x in xs
                try
                    local f
                    f(y) = lin_stab_line_fraction_checker(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k, ω = ω )
                    Pe = find_zero(f, (0.,  100.))
                    push!(Y,Pe)
                    push!(X,x)
                catch
                end
        end
        lab = latexstring("\$ \\omega = $(ω) \$")
        ax.plot(X,Y, label = lab) #color = "red"
    end
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
    ax.set_aspect(Δρ/ΔPe)
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