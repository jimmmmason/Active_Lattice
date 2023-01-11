cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")


#visualise data
using TensorOperations, PyPlot, PyCall
@pyimport matplotlib.animation as anim

function animate_pdes(param,t_saves,fa_saves,fp_saves)
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
        plot_pde_mass(fig,ax1,param,density)
        plot_pde_mag( fig,ax2,param,density)
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

function plot_pde_mag(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
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
    colmap = PyPlot.plt.cm.viridis_r
    x = Δx:Δx:1
    y = Δx:Δx:1
    xx = [x̃ for x̃ ∈ x, ỹ ∈ y]'
    yy = [ỹ for x̃ ∈ x, ỹ ∈ y]'
    ax.streamplot(xx, yy, m[:,:,1]', m[:,:,2]', color = absmag', cmap = colmap, density = 1)#2.5
    norm1 = matplotlib.colors.Normalize(vmin=minimum(absmag), vmax= maximum(absmag));
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    #display(fig)
    return fig
end

function plot_pde_mass(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
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
    colmap = PyPlot.plt.cm.viridis_r
    norm1 = matplotlib.colors.Normalize(vmin=minimum(ρ), vmax= maximum(ρ) );
    ax.contourf(Δx:Δx:1, Δx:Δx:1,ρ'; levels = 25, norm = norm1, cmap = colmap )
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
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
    @unpack T, save_interval = param
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    animate_pdes(param,t_saves,fa_saves,fp_saves)
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
    colmap = PyPlot.plt.cm.viridis_r
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
    colmap = PyPlot.plt.cm.viridis_r
    norm1 = matplotlib.colors.Normalize(vmin=minimum(ρ), vmax= maximum(ρ) );
    ax.contourf(Δx:Δx:1, Δx:Δx:1,rho'; levels = 25, norm = norm1, cmap = colmap )
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
    #fig
    return fig
end 

# 1d wave plots 

function make_phase_video_1d(param; frames = 100.)
    @unpack T, save_interval = param
    save_interval = T/frames
    t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = save_interval)
    animate_pdes_1d(param,t_saves,fa_saves,fp_saves)
end


function animate_phase_pdes_1d(param,t_saves,fa_saves,fp_saves)
    @unpack name, λ, ρa, ρp, Nx, Nθ, δt = param
    frames = length(t_saves)-1
    fig, axs = plt.subplots(1, 2, figsize=(12,5))
    function makeframe(i)
        clf()
        vid_phase_pde_plot_1d(fig, axs, param, t_saves, fa_saves, fp_saves, i)
        return fig
    end
    interval = Int64(round(10000/frames))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    T = t_saves[frames+1]
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/vids/pde_phase_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)/time=$(round(T; digits = 5))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dpi= 100, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
end 

function vid_phase_pde_plot_1d(fig::Figure, axs ::Array{PyObject,2}, param::Dict{String,Any}, t_saves, fa_saves, fp_saves, i)
    @unpack Nx, Nθ, ρa, ρp = param
    ρa_saves, ρp_saves = spatial_densities(fa_saves, fp_saves; Nx= Nx, Nθ = Nθ)
    phasea_saves, phasep_saves = spatial_fourier_modes(ρa_saves, ρp_saves; Nx = 50) 

    axs[1].plot((1:Nx)/Nx,ρa_saves[i].-ρa, color = "red")
    axs[1].plot((1:Nx)/Nx,ρp_saves[i].-ρp, color = "black")
    #axs[1].xaxis.set_ticks(xtic)
    #axs[1].yaxis.set_ticks(ytic)
    #axs[1].axis(axlim)
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("ρₐ-ϕₐ,ρₚ-ϕₚ")
    axs[1].set_title("ℓ= $(1/sqrt(Dθ)), χ = $(χ)")

    axs[2].plot((1:Nx)/Nx,mag_1d(fa_saves[i]; Nθ = Nθ), color = "red")
    #axs[1].xaxis.set_ticks(xtic)
    #axs[1].yaxis.set_ticks(ytic)
    #axs[1].axis(axlim)
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("m")

    axs[3].plot(t_saves, abs.(phasea_saves), color = "red", label = "Active")
    axs[3].plot(t_saves, abs.(phasep_saves), color = "black", label = "Passive")
    #axs[3].xaxis.set_ticks(xtic)
    #axs[3].yaxis.set_ticks(ytic)
    #axs[3].axis(axlim)
    axs[3].set_xlabel("t")
    #ax.set_ylabel("Pe")
    axs[3].set_ylabel("Amlitude")
    axs[3].legend(loc = "upper right")

    axs[4].plot(t_saves, phase_args(angle.(phasea_saves)), color = "red")
    axs[4].plot(t_saves, phase_args(angle.(phasep_saves)), color = "black")
    #axs[2].xaxis.set_ticks(xtic)
    #axs[2].yaxis.set_ticks(ytic)
    #axs[2].axis(axlim)
    axs[4].set_xlabel("t")
    axs[4].set_ylabel("Phase")
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

function mag_1d(fa; Nθ = 20)
    eθ = cos.(1:Nθ*2π/Nθ)
    @tensor begin
        m[a] := 2π *fa[a,θ]*eθ[θ]/Nθ
    end
    return m
end

println("booted")