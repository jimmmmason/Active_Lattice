cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")


#visualise data
using TensorOperations, PyPlot, PyCall
@pyimport matplotlib.animation as anim

function animate_pdes(param,t_saves,pde_saves)
    @unpack name, L, λ, ρa, ρp, Δt = param
    frames = length(η_saves)-1
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    makeframe(i) = plot_eta(fig,ax,param, t_saves[i+1], η_saves[i+1])
    interval = Int64(round(5/Δt))
    myanim = anim.FuncAnimation(fig, makeframe, frames=frames, interval=interval)
    # Convert it to an MP4 movie file and saved on disk in this format.
    pathname = "/home/jm2386/Active_Lattice/plots/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)";
    mkpath(pathname)
    filename = "/home/jm2386/Active_Lattice/plots/pde_vids/$(name)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_Δt=$(Δt)/time=$(round(T; digits = 5))_size=$(L)_active=$(ρa)_passive=$(ρp)_lamb=$(λ).mp4";
    myanim[:save](filename, bitrate=-1, dπ= 100, extra_args=["-vcodec", "libx264", "-πx_fmt", "yuv420p"])
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
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    cp = ax.contourf(Δx:Δx:1, Δx:Δx:1,absmag; levels = 100, set_cmap = winter() )
    fig.colorbar(cp)
    #fig
    return fig
end 

function plot_pde_mass(fig::Figure, ax::PyObject, param::Dict{String,Any}, density::Dict{String,Any})
    @unpack λ, Nx, Nθ, S, ρa, ρp= param
    @unpack fa, fp, t = density
    #fig, ax = PyPlot.subplots(figsize =(10, 10))
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
        ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), Pe = $(λ), t = $(t)")
    # Plot points
    cp = ax.contourf(Δx:Δx:1, Δx:Δx:1,ρ; levels = 100, set_cmap = winter() )
    fig.colorbar(cp)
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

function plot_stab(fig, ax, stabdata; ρs = 0.05:0.05:0.95 ,xs = collect(0.01:0.01:0.99), xtic = 0:0.2:1, ytic = 0:10:100, axlim = [0., 1., 0., 100.] )
    stable_points = []
    unstable_points = []
    unsure_points = []
    for ρ in ρs
            @unpack stable, unstable, unsure = stabdata["ρ = $(ρ)"]
            for λ ∈ stable
                    append!(stable_points,  [ρ; λ])
            end
            for λ ∈ unstable
                    append!(unstable_points, [ρ; λ])
            end
            for λ ∈ unsure
                    append!(unsure_points, [ρ; λ])
            end
    end             
    stable_points   = reshape(stable_points, (2,Int64(length(stable_points)/2)) )
    unstable_points = reshape(unstable_points, (2,Int64(length(unstable_points)/2)) )
    unsure_points = reshape(unsure_points, (2,Int64(length(unsure_points)/2)) )

    lower_bound = λsym.(xs; Dθ = 10.)

    ax.plot(xs,lower_bound)

    ax.errorbar(stable_points[1,:],stable_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "green",
    alpha=0.8,
    )

    ax.errorbar(unstable_points[1,:],unstable_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "red",
    alpha=0.8,
    )

    ax.errorbar(unsure_points[1,:],unsure_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "blue",
    alpha=0.8,
    )

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    #ax.set_title("ρₐ = $(ρa),  ρₚ = $(ρp), λ = $(λ), t = $(t), Φ = $(Φ)")

end

println("booted")