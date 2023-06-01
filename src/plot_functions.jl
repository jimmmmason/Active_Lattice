cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")


#visualise data
using TensorOperations, PyPlot, PyCall
@pyimport matplotlib.animation as anim

function animate_pdes(param,t_saves,fa_saves,fp_saves; video_length = 10000.)
    @unpack name, λ, ρa, ρp,Nx, Nθ, δt,χ, cbar_min, cbar_max = param
    frames = length(t_saves)-1
    fig, axs = plt.subplots(1, 2, figsize=(12,5))
    function makeframe(i)
        clf()
        fa = fa_saves[i+1]
        fp = fp_saves[i+1]
        t  = t_saves[i+1]
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        @unpack λ, Nx, Nθ, S, ρa, ρp, Dθ= param
        density = Dict{String,Any}()
        @pack! density = fa, fp, t
        plot_pde_mass(fig,ax1,param,density; cbar_max = cbar_max, cbar_min = cbar_min)
        plot_pde_mag( fig,ax2,param,density)
        l = 1/sqrt(Dθ)
        fig.suptitle("ϕ = $(round(ρa+ρp; digits = 3)),  χ = $(round(χ; digits = 3)), Pe = $(round(λ/sqrt(Dθ); digits = 3)), l = $(round(l; digits = 3)), t = $(round(t; digits = 3))",size  =15. )
        return fig
    end
    interval = Int64(round(video_length/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    T = t_saves[frames+1]
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(T; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function plot_pde_mag(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any}; scale = "local", cbar = true)
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #fig, ax = PyPlot.subplots(figsize =(10, 10))
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := 2π *fa[a,b,θ]*eθ[θ,d]/Nθ
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
        #ax.yaxis.set_ticks([])
        ax.axis([Δx, 1., Δx, 1.])
        ax.set_aspect("equal")
        #ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    x = Δx:Δx:1
    y = Δx:Δx:1
    xx = [x̃ for x̃ ∈ x, ỹ ∈ y]'
    yy = [ỹ for x̃ ∈ x, ỹ ∈ y]'
    ax.streamplot(xx, yy, m[:,:,1]', m[:,:,2]', color = absmag', cmap = colmap, density = 1)#2.5
    if scale == "local"
        norm1 = matplotlib.colors.Normalize(vmin=minimum(absmag), vmax= maximum(absmag));
    else
        norm1 = matplotlib.colors.Normalize(vmin=scale[1], vmax= scale[2]);
    end
    if cbar
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    end
    #display(fig)
    return fig
end

function plot_pde_mass(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any}; scale = "not_relative", cbar = true, cbar_max = 1.0, cbar_min = 0.0)
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #=
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    =#
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,b,d] := 2π *fa[a,b,θ]*eθ[θ,d]/Nθ
    end
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    #absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
    #ax.yaxis.set_ticks([])
    ax.axis([Δx, 1., Δx, 1.])
    ax.set_aspect("equal")
    #ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    if scale == "relative"
        norm1 = matplotlib.colors.Normalize(vmin= cbar_min, vmax= cbar_max);
    else
        norm1 = matplotlib.colors.Normalize(vmin=minimum(ρ), vmax= maximum(ρ) );
    end
    ax.contourf(Δx:Δx:1, Δx:Δx:1,ρ'; levels = 25, norm = norm1, cmap = colmap )
    if cbar
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    end
    #fig
    return fig
end 

function plot_error(fig::Figure, ax::PyObject,param, t_saves, fa_saves, fp_saves; t_max = 0.3, y_max = 0.2)
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
    ax.plot(t_saves,dist_saves)
    ax.set_title("ρₐ = $(ρa), Pe = $(λ)")
    ax.axis([0, t_max, 0., y_max])
end

#=
fig, axs = plt.subplots(1, 2, figsize=(12,5))
plot_pde(fig,axs,param,density)
display(fig)
=#
function plot_pde(fig::Figure, axs ::Array{PyObject,1} , param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp, Dθ= param
    @unpack fa, fp, t = density
    plot_pde_mass(fig,axs[1],param,density)
    plot_pde_mag( fig,axs[2],param,density)
    l = 1/sqrt(Dθ)
    ρ = ρa+ρp
    χ = ρa/ρ
    fig.suptitle("ϕ = $(round(ρ; digits = 2)),  χ = $(round(χ; digits = 2)), Pe = $(round(λ/sqrt(Dθ); digits = 3)), l = $(round(l; digits = 3)), t = $(round(t; digits = 3))",size  =15. )
    return fig
end 

function plot_pde_lite!(fig,param::Dict{String,Any}, t, fa, fp )
    #PyPlot.close("all")
    clf()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    @unpack λ, Nx, Nθ, S, ρa, ρp, Dθ= param
    density = Dict{String,Any}()
    @pack! density = fa, fp, t
    plot_pde_mass(fig,ax1,param,density)
    plot_pde_mag( fig,ax2,param,density)
    l = 1/sqrt(Dθ)
    fig.suptitle("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(round(λ/sqrt(Dθ); digits = 3)), l = $(round(l; digits = 3)), t = $(round(t; digits = 3))",size  =15. )
    
    return fig
end 

function Stability_condition(ϕ,Dx,Pe,Dθ)
    v0 = Pe*sqrt(Dθ);
    ds = self_diff(ϕ);
    dp = self_diff_prime(ϕ);
    c1 = Dx *ds;
    c2 = Dx *(1 - ds)/(2*π);
    c3 = -   v0 * (1 - ϕ - ds);
    c4 = -   v0 * (dp*ϕ)/2;
    c5 = -   v0*(ds)/2;
    c6 = Dθ; 
    ω2 = (2*π)^2; 
    μ = c1*ω2 + c6;
    μ0 = c1*ω2;
    var1 = 1 + (μ^2)/( 2*ω2*c5^2);
    Λ = -var1 + sqrt(var1^2 - 1); 
    return μ - (- ω2*(c4 + c5)*(c3 + c5)/(μ0 + c2*ω2) - ω2*(1 + Λ )/μ)
end

function make_video(param)
    @unpack T, save_interval, video_length = param
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    animate_pdes(param,t_saves,fa_saves,fp_saves; video_length=video_length)
end


###

function make_video_1d(param; frames = 100.)
    @unpack T, save_interval = param
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    animate_pdes_1d(param,t_saves,fa_saves,fp_saves)
end

function animate_pdes_1d(param,t_saves,fa_saves,fp_saves)
    @unpack name, λ, ρa, ρp,Nx, Nθ, δt = param
    frames = length(t_saves)-1
    fig, axs = plt.subplots(1, 2, figsize=(12,5))
    function makeframe(i)
        clf()
        fa = fa_saves[i+1]
        fp = fp_saves[i+1]
        t  = t_saves[i+1]
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        @unpack λ, Nx, Nθ, S, ρa, ρp, Dθ= param
        density = Dict{String,Any}()
        @pack! density = fa, fp, t
        plot_pde_mass_1d(fig,ax1,param,density)
        plot_pde_mag_1d( fig,ax2,param,density)
        l = 1/sqrt(Dθ)
        fig.suptitle("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(round(λ/sqrt(Dθ); digits = 3)), l = $(round(l; digits = 3)), t = $(round(t; digits = 3))",size  =15. )
        return fig
    end
    interval = Int64(round(10000/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    T = t_saves[frames+1]
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(T; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function plot_pde_1d(fig::Figure, axs ::Array{PyObject,1} , param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp, Dθ= param
    @unpack fa, fp, t = density
    plot_pde_mas_1d(fig,axs[1],param,density)
    plot_pde_mag_2d( fig,axs[2],param,density)
    l = 1/sqrt(Dθ)
    fig.suptitle("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(round(λ/sqrt(Dθ); digits = 3)), l = $(round(l; digits = 3)), t = $(round(t; digits = 3))",size  =15. )
    return fig
end

function plot_pde_mag_1d(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #fig, ax = PyPlot.subplots(figsize =(10, 10))
    #collect data
    eθ = [cos.(S*2π/Nθ) sin.(S*2π/Nθ)]
    @tensor begin
        m[a,d] := 2π *fa[a,θ]*eθ[θ,d]/Nθ
    end
    #ρ = fp + sum(fa; dims =2)[:,1].*(2*π/Nθ)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
        #ax.yaxis.set_ticks([])
        ax.axis([Δx, 1., Δx, 1.])
        ax.set_aspect("equal")
        #ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    x = Δx:Δx:1
    y = Δx:Δx:1
    xx = [x̃ for x̃ ∈ x, ỹ ∈ y]'
    yy = [ỹ for x̃ ∈ x, ỹ ∈ y]'
    mm = zeros(Nx,Nx,2)
    for x in 1:Nx
        mm[:,x,:] = m
    end
    absmag  = sqrt.(mm[:,:,1].^2+mm[:,:,2].^2)
    ax.streamplot(xx, yy, mm[:,:,1]', mm[:,:,2]', color = absmag', cmap = colmap, density = 1)#2.5
    norm1 = matplotlib.colors.Normalize(vmin=minimum(absmag), vmax= maximum(absmag));
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    #display(fig)
    return fig
end

function plot_pde_mass_1d(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #=
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    =#
    #collect data
    ρ = fp + sum(fa; dims =2)[:,1].*(2*π/Nθ)
    rho = zeros(Nx,Nx)
    for x in 1:Nx
        rho[:,x] = ρ
    end
    #absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    #Φ = round(translation_invariance(η;Ω = Ω,L= L); digits =4)
    Δx = 1/Nx
    #figure configuration
    #ax.xaxis.set_ticks([])
    #ax.yaxis.set_ticks([])
    ax.axis([Δx, 1., Δx, 1.])
    ax.set_aspect("equal")
    #ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    norm1 = matplotlib.colors.Normalize(vmin=minimum(ρ), vmax= maximum(ρ) );
    ax.contourf(Δx:Δx:1, Δx:Δx:1,rho'; levels = 25, norm = norm1, cmap = colmap )
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    #fig
    return fig
end 

# 1d wave plots 

function make_phase_video_1d(param; frames = 100)
    @unpack T, save_interval = param
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    frames = Int64(round(length(t_saves)))-1
    animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves; frames = frames-1)
end

function pmap_make_phase_video_1d(param)
    @unpack T, save_interval = param
    frames = T*100
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    frames = Int64(round(length(t_saves)))-1
    animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves; frames = frames-1)
end

function animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves; frames = 99)
    @unpack name, λ, ρa, ρp, Nx, Nθ, δt, Dθ, χ, γ = param
    fig, axs = plt.subplots(2, 1, figsize=(10,10))
    function makeframe(i)
        clf()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        axs = ax1, ax2
        vid_pde_plot_1d(fig, axs, param, t_saves, fa_saves, fp_saves, i+1)
        return fig
    end
    interval = 5*Int64(round(20000/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    T = t_saves[Int64(round((frames+1)))]
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(T; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function test_vid_phase_pde_plot_1d(fig::Figure, axs, param::Dict{String,Any}, t_saves, fa_saves, fp_saves, i; speed_factor = 0.1)
    @unpack Nx, Nθ, ρa, ρp, χ, Dθ, Dx, k, γ = param
    ρa_saves, ρp_saves = spatial_densities(fa_saves, fp_saves; Nx= Nx, Nθ = Nθ)
    phasea_saves, phasep_saves = spatial_fourier_modes(ρa_saves, ρp_saves; Nx = Nx) 
    phase2a_saves, phase2p_saves = spatial_fourier_mode2s(ρa_saves, ρp_saves; Nx = Nx) 

    j = Int64(round(i/speed_factor))

    axs[1].plot((1:Nx)/Nx,ρa_saves[j].-ρa, color = "red")
    axs[1].plot((1:Nx)/Nx,ρp_saves[j].-ρp, color = "black")
    #axs[1].xaxis.set_ticks(xtic)
    #axs[1].yaxis.set_ticks(ytic)
    #axs[1].axis(axlim)
    #axs[1].axis([0., 1., minimum(minimum.(ρa_saves))-ρa,maximum(maximum.(ρa_saves))-ρa])
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("ρₐ-ϕₐ,ρₚ-ϕₚ")
    axs[1].set_title("t = $(round(t_saves[j]; digits = 3))")


    #=
    K = collect(0:1:(k-1))
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k=k, γ= γ)
    ω = 2*π
    a, A = ap_MathieuEigen(matrix)
    lin_fa = zeros(Nx, Nθ)
    lin_fp = zeros(Nx)
    a = a[k+1]
    for x in 1:Nx, θ in 1:Nθ
        lin_fa[x...,θ...] = real.( dot(A[2:1:(k+1),k+1],cos.(θ*K*(2*π/Nθ)))*exp(a*t_saves[j]+im*x*ω/Nx) )
        lin_fp[x...] = real.(A[1,k+1]*exp(a*t_saves[j]+im*x*ω/Nx));
    end
    c = dist_from_unif_1d(param, lin_fa.+ρa/(2*π), lin_fp.+ρp)
    lin_ρa = sum(lin_fa; dims =2)[:,1].*(2*π/Nθ)
    axs[1].plot((1:Nx)/Nx,exp(real(a)*t_saves[j])*δ*lin_ρa/c, color = "green")
    axs[1].plot((1:Nx)/Nx,exp(real(a))*t_saves[j]*δ*lin_fp/c, color = "blue")
    
    
        K = collect(0:1:(k-1))
        matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k=k, γ= γ)
        ω = 2*π
        a, A = ap_MathieuEigen(matrix)
        Pa = (x,θ) -> real.( dot(A[2:1:(k+1),k+1],cos.(θ*K*(2*π/Nθ)))*exp(a*t_saves[j]-im*x*ω/Nx) )
        Pp = (x) -> real.(A[1,k+1]*exp(a*t_saves[j]-im*x*ω/Nx));
        a = a[k]

        perta = zeros(Nx,Nθ)
        pertp = zeros(Nx)
        for x₁ in 1:Nx, θ in 1:Nθ
            perta[x₁, θ] = Pa(x₁, θ);
        end
        for x₁ in 1:Nx
            pertp[x₁] = Pp(x₁);
        end

        c = dist_from_unif_1d(param, perta.+ρa/(2*π), pertp.+ρp)
        lin_fa = δ*perta/c
        lin_fp = δ*pertp/c
        lin_ρa = sum(lin_fa; dims =2)[:,1].*(2*π/Nθ)


        axs[1].plot((1:Nx)/Nx,exp(real(a)*t_saves[j])*lin_ρa, color = "pink")
        axs[1].plot((1:Nx)/Nx,exp(real(a)*t_saves[j])*lin_fp, color = "blue")
    =#

    axs[2].plot((1:Nx)/Nx,mag_1d(fa_saves[j]; Nθ = Nθ), color = "red")
    #axs[1].xaxis.set_ticks(xtic)
    #axs[1].yaxis.set_ticks(ytic)
    #axs[1].axis(axlim)
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("m")
    #axs[2].axis([0., 1., minimum(minimum.(mag_1d.(fa_saves; Nθ = Nθ))),maximum(maximum.(mag_1d.(fa_saves; Nθ = Nθ)))])

    axs[3].plot(t_saves, abs.(phasea_saves), color = "red", label = "Active")
    axs[3].plot(t_saves, abs.(phasep_saves), color = "black", label = "Passive")
    #axs[3].xaxis.set_ticks(xtic)
    #axs[3].yaxis.set_ticks(ytic)
    #axs[3].axis(axlim)
    axs[3].set_xlabel("t")
    #ax.set_ylabel("Pe")
    axs[3].set_ylabel("1ˢᵗ mode amplitude")
    #axs[3].legend(loc = "upper right")
    axs[3].scatter(t_saves[j],abs.(phasea_saves)[j],color = "red")
    axs[3].scatter(t_saves[j],abs.(phasep_saves)[j],color = "black")
    #axs[3].set_aspect(maximum(t_saves)/maximum(abs.(phasea_saves)))

    axs[4].plot(t_saves, phase_args(angle.(phasea_saves)), color = "red")
    axs[4].plot(t_saves, phase_args(angle.(phasep_saves)), color = "black")
    #axs[2].xaxis.set_ticks(xtic)
    #axs[2].yaxis.set_ticks(ytic)
    #axs[2].axis(axlim)
    axs[4].set_xlabel("t")
    axs[4].set_ylabel("1ˢᵗ mode phase")
    axs[4].scatter(t_saves[j],phase_args(angle.(phasea_saves))[j],color = "red")
    axs[4].scatter(t_saves[j],phase_args(angle.(phasep_saves))[j],color = "black")
    #axs[4].set_aspect(maximum(t_saves)/(maximum(phase_args(angle.(phasea_saves)))-minimum(phase_args(angle.(phasea_saves)))))

    axs[5].plot(t_saves, abs.(phase2a_saves), color = "red", label = "Active")
    axs[5].plot(t_saves, abs.(phase2p_saves), color = "black", label = "Passive")
    #axs[3].xaxis.set_ticks(xtic)
    #axs[3].yaxis.set_ticks(ytic)
    #axs[3].axis(axlim)
    axs[5].set_xlabel("t")
    #ax.set_ylabel("Pe")
    axs[5].set_ylabel("2ⁿᵈ mode amlitude")
    axs[5].legend(loc = "upper right")
    axs[5].scatter(t_saves[j],abs.(phase2a_saves)[j],color = "red")
    axs[5].scatter(t_saves[j],abs.(phase2p_saves)[j],color = "black")
    #axs[3].set_aspect(maximum(t_saves)/maximum(abs.(phasea_saves)))

    axs[6].plot(t_saves, phase_args(angle.(phase2a_saves)), color = "red")
    axs[6].plot(t_saves, phase_args(angle.(phase2p_saves)), color = "black")
    #axs[2].xaxis.set_ticks(xtic)
    #axs[2].yaxis.set_ticks(ytic)
    #axs[2].axis(axlim)
    axs[6].set_xlabel("t")
    axs[6].set_ylabel("2ⁿᵈ mode phase")
    axs[6].scatter(t_saves[j],phase_args(angle.(phase2a_saves))[j],color = "red")
    axs[6].scatter(t_saves[j],phase_args(angle.(phase2p_saves))[j],color = "black")

    fig.suptitle("ℓ= $(1/sqrt(Dθ)), χ = $(χ), ρ = $(round(ρa+ρp;digits = 3))", fontsize = 20, y = 1.05)
    fig.tight_layout()
end

function vid_pde_plot_1d(fig::Figure, axs, param::Dict{String,Any}, t_saves, fa_saves, fp_saves, i)
    @unpack Nx, Nθ, ρa, ρp, χ, Dθ, Dx, k, γ,Pe = param
    ρa_saves, ρp_saves = deepcopy(spatial_densities(fa_saves, fp_saves; Nx= Nx, Nθ = Nθ))

    push!(ρa_saves[i], ρa_saves[i][1])
    push!(ρp_saves[i], ρp_saves[i][1])

    ρsum = ρp_saves[i]+ρa_saves[i]

    #sum(length.(fp_saves).-128)
    #=
    axs[1].plot((1:Nx)/Nx,ρa_saves[i], color = "red")
    axs[1].plot((1:Nx)/Nx,ρp_saves[i]+ρa_saves[i], color = "black")
    axs[1].plot((1:Nx)/Nx,ρp_saves[i], color = "blue")
    =#
    axs[1].plot((0:1:Nx)/Nx,ρa_saves[i].-ρa, color = "red", label = L"\rho^a - \phi^a")
    axs[1].plot((0:1:Nx)/Nx,ρsum.-ρa.-ρp, color = "black", label = L"\rho - \phi")
    axs[1].plot((0:1:Nx)/Nx,ρp_saves[i].-ρp, color = "blue", label = L"\rho^p - \phi^p")

    axs[1].xaxis.set_ticks(0.:0.2:1.0)
    axs[1].xaxis.set_tick_params(labelsize=15)
    axs[1].yaxis.set_tick_params(labelsize=15)
    rhomax = maximum(maximum(ρa_saves).-ρa)+maximum(maximum(ρp_saves).-ρp)
    axs[1].axis([0., 1., -rhomax , rhomax])
    #axs[1].axis([0., 1., min(minimum(minimum(ρa_saves)),minimum(minimum(ρp_saves))),maximum(maximum( ρa_saves+ρp_saves ))])
    axs[1].set_xlabel(L"x",fontsize=20)
    #axs[1].set_ylabel(L"\rho,",fontsize=20)
    title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ), \\phi = $(ρa+ρp), \\mathrm{Pe} = $(round(Pe; digits = 3)), t = $(round(t_saves[i]; digits = 3))\$")
    axs[1].set_title(title,fontsize=20)

    mat1 = zeros(1, Nx+1)
    mat2= zeros(1, Nx+1)
    mags = mag_1d(fa_saves[i]; Nθ = Nθ)
    push!(mags,mags[1])
    mat1[1, :] = mags
    mat2[1, :] = mags.*(-ρsum.+1)

    #colmap = PyPlot.plt.cm.seismic
    colmap = PyPlot.plt.cm.PRGn
    norm1 = matplotlib.colors.Normalize(vmin= -maximum(abs.(mags)) , vmax= maximum(abs.(mags)) )
    #norm2 = matplotlib.colors.Normalize(vmin= minimum(mags/10) , vmax= maximum(mags)/10 )

    axs[2].matshow(mat1; norm = norm1,  cmap = colmap, extent = [0., 1., 0., 0.1])
    #axs[3].matshow(mat2; norm = norm2,  cmap = colmap, extent = [0., 1., 0., 0.1])

    axs[2].set_aspect(1.)
    #axs[3].set_aspect(1.)

    axs[2].xaxis.set_ticks(0.:0.2:1.0)
    axs[2].yaxis.set_ticks([])
    axs[2].xaxis.set_tick_params(labelsize=15)
    axs[2].xaxis.tick_bottom()
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    #ax.set_xlabel(L"x",fontsize=20)

    axs[2].set_ylabel(L"\mathbf{p}", fontsize=20, rotation=0)
    axs[2].yaxis.set_label_coords(-.05, .5)

    lines, labels = axs[1].get_legend_handles_labels()
    fig.tight_layout()
    ldg = fig.legend(lines, labels, loc = "center right", fontsize=20, bbox_to_anchor = (0.25, 0.25, 1, 1),
    bbox_transform = plt.gcf().transFigure)

    #fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = axs[2:1:3], fraction = 0.0455)
    return ldg
end

function phase_pde_plot_1d_test(fig::Figure, axs ::Array{PyObject,2}, param::Dict{String,Any}, t_saves, fa_saves, fp_saves)
    @unpack Nx, Nθ = param
    ρa_saves, ρp_saves = spatial_densities(fa_saves, fp_saves; Nx= 50, Nθ = 20)
    phasea_saves, phasep_saves = spatial_fourier_modes(ρa_saves, ρp_saves; Nx = 50) 

    axs[1].plot(t_saves,ρa_saves,color = "red")
    axs[1].plot(t_saves,ρa_saves,color = "black")

    axs[2].plot(t_saves,ρa_saves,color = "red")
    axs[2].plot(t_saves,ρa_saves,color = "black")
end

function spatial_densities(fa_saves, fp_saves; Nx= 50, Nθ = 20)
    ρa_saves, ρp_saves = [], fp_saves
    for fa in fa_saves
        ρa = sum(fa; dims =2)[:,1].*(2*π/Nθ)
        push!(ρa_saves,ρa)
    end
    return ρa_saves, ρp_saves
end

function phase_args(phase_saves)
    cnts_phase = [phase_saves[1]]
    L = length(phase_saves)
    k = 0
    for i in 2:L
        x = phase_saves[i]-phase_saves[i-1]
        if abs.(x)>3
            push!(cnts_phase,phase_saves[i]-2*π*sign(x)-2*π*k)
            k+= sign(x)
        else 
            push!(cnts_phase,phase_saves[i]-2*π*k)
        end
    end
    return cnts_phase
end

function spatial_fourier_modes(ρa_saves, ρp_saves; Nx = 50)
    phasea_saves, phasep_saves = [], []
    for ρa in ρa_saves
        phase = spatial_fourier_mode(ρa; Nx = 50) 
        push!(phasea_saves,phase)
    end
    for ρp in ρp_saves
        phase = spatial_fourier_mode(ρp; Nx = 50) 
        push!(phasep_saves,phase)
    end
    return phasea_saves, phasep_saves
end

function spatial_fourier_mode(ρ; Nx = 50)
    phase = 0.
    for x in 1:Nx
        phase += ρ[x]*exp(-2*π*im*x/Nx)/Nx
    end
    return phase
end

function spatial_fourier_mode2s(ρa_saves, ρp_saves; Nx = 50)
    phasea_saves, phasep_saves = [], []
    for ρa in ρa_saves
        phase = spatial_fourier_mode2(ρa; Nx = Nx) 
        push!(phasea_saves,phase)
    end
    for ρp in ρp_saves
        phase = spatial_fourier_mode2(ρp; Nx = Nx) 
        push!(phasep_saves,phase)
    end
    return phasea_saves, phasep_saves
end

function spatial_fourier_mode2(ρ; Nx = 50)
    phase = 0.
    for x in 1:Nx
        phase += ρ[x]*exp(-2*π*im*2*x/Nx)/Nx
    end
    return phase
end

function mag_1d(fa; Nθ = 20)
    eθ = cos.(1:Nθ*2π/Nθ)
    @tensor begin
        m[a] := 2π *fa[a,θ]*eθ[θ]/Nθ
    end
    return m
end
##


function plot_sim_mag(fig::Figure, ax::PyObject, param::Dict{String,Any}, t; scale = "local", cbar = true)
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ= param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/pol_save/r=$(r)_size=$(L)/t=$(round(t;digits=2))_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    data = wload(filename)
    @unpack h = data
    h = reshape(h,(L,L))
    m = zeros(L,L,2)
    for x in 1:L, y in 1:L
        m[x,y,:] = h[x,y]
    end
    absmag  = sqrt.(m[:,:,1].^2+m[:,:,2].^2)
    t = round(t; digits=5)
    Δx = 1/L
    ax.axis([Δx, 1., Δx, 1.])
    ax.set_aspect("equal")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    x = Δx:Δx:1
    y = Δx:Δx:1
    xx = [x̃   for x̃ ∈ x, ỹ ∈ y]'
    yy = [1+Δx- ỹ for x̃ ∈ x, ỹ ∈ y]'
    #ax.matshow(ρ; norm = norm1,  cmap = colmap, extent = [Δx, 1., Δx, 1.])
    ax.streamplot(xx, yy, m[:,:,1], -m[:,:,2], color = absmag, cmap = colmap, density = 1)#2.5
    if scale == "local"
        norm1 = matplotlib.colors.Normalize(vmin=minimum(absmag), vmax= maximum(absmag));
    else
        norm1 = matplotlib.colors.Normalize(vmin=scale[1], vmax= scale[2]);
    end
    if cbar
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    end
    #display(fig)
    return fig
end

function plot_sim_mass(fig::Figure, ax::PyObject, param::Dict{String,Any},t; scale = "relative", cbar = true)
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/hist_save/r=$(r)_size=$(L)/t=$(round(t;digits=2))_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    data = wload(filename)
    @unpack h = data
    ρ = reshape(h,(L,L))
    #
    t = round(t; digits=5)
    Δx = 1/L
    ax.axis([Δx, 1., Δx, 1.])
    ax.set_aspect("equal")
    # Plot points
    colmap = PyPlot.plt.cm.viridis
    if scale == "relative"
        norm1 = matplotlib.colors.Normalize(vmin=0., vmax= 1.);
    else
        norm1 = matplotlib.colors.Normalize(vmin=minimum(ρ), vmax= maximum(ρ) );
    end
    ax.matshow(ρ; norm = norm1,  cmap = colmap, extent = [Δx, 1., Δx, 1.]) #Δx:Δx:1, Δx:Δx:1,
    if cbar
        fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    end
    #fig
    return fig
end 
println("booted")