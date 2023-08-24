cd("/home/jm2386/Active_Lattice/")
using DrWatson
using Distributed
@quickactivate "Active_Lattice"
println("Loading ...")
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/sim_functions.jl")
include("/home/jm2386/Active_Lattice/src/lin_stab_solver.jl")
##
function pde_param_fraction(; name = "test", D =1., Dx = 1., Pe =1., Dθ = 10, ρ= 0.5, χ = 1.0, Nx = 100, Nθ = 20, δt = 1e-5, T= 0.001, save_interval = 0.01, max_steps = 1e8, max_runs = 6, λ_step = 10., λmax = 100., λs = 20.:20.:100., pert = "n=1", δ = 0.01, k=20,γ = 0.0, video_length = 10000., cbar_max = 1.0, cbar_min = 0.0, frames = 1000, pert_interval = 5.)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    ρp = 0.0
    λ = Pe*sqrt(Dθ)
    ρp = (1-χ)*ρ
    ρa = χ*ρ
    param = Dict{String,Any}()
    @pack! param = k, name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E, Dθ, T, save_interval, max_steps, max_runs, λ_step, λmax, λs, pert, δ, Pe, Dx, χ, ρ, γ, video_length, cbar_max, cbar_min, frames, pert_interval
    return param
end
function sim_param_fraction(;  name = "test", D =1. , Pe =1. ,ρ = 0.5, χ = 0.5, L=10, d=2, Δt = 0.01, Dθ =10., T=1.0, γ = 0.)
    param = Dict{String,Any}()
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    λ = Pe*sqrt(Dθ)/2 # see erignoux paper this lambda is half the pde lambda but gives same Pe in final model 
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    site_distribution = fill([1-ρa-ρp, ρa, ρp],(L,L))
    POSITIONS = reshape(collect(1:(L^2*4)), (L,L,4))
    function angles(x,n) 
        if n == 1
            return 2*π*rand()
        else
            return -1
        end
    end
    rates(n,m,i) = 0.
    if γ>0.
        function rates(n,m,i)
            E = [[1,0],[0,1],[0,-1],[-1,0],]
            rate = 0.
            if n[1]==0.
                rate = 0.
            elseif n[1]==2.
                rate = L^2*D
            elseif n[1]==1.
                rate = L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
            end
            if m[1]>0.
                return rate*γ
            else
                return rate
            end
        end
    else
        function rates(n,m,i)
            E = [[1,0],[0,1],[0,-1],[-1,0],]
            if m[1]>0.
                return 0.
            else
                if n[1]==0.
                    return 0.
                elseif n[1]==2.
                    return L^2*D
                elseif n[1]==1.
                    return L^2*D + L*λ*E[i]⋅[cos(n[2]),sin(n[2])] 
                end
            end
        end
    end

    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ, T, POSITIONS, χ, ρ, Pe
    return param
end
#
function lin_stab_line_fraction(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return real(amin)
end
function a_lin_stab_line_fraction(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = a_MathieuEigen_lite(matrix; k = k)
    return real(amin)
end
using Roots
function crit_Pe(ρ,χ; Dx =1. ,Pe_max = 40, Dθ =100.,k=40)
    f(x) = a_lin_stab_line_fraction(ρ,χ; Dx =Dx ,Pe = x, Dθ =Dθ, k=k)
    return find_zero(f,(0,Pe_max))
end
function lin_imaginary_fraction(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end
function a_lin_imaginary_fraction(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = a_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end
#
function find_stab_data(;stabdata = Dict{String,Any}(), ρs = 0.4:0.05:1.0, Pes = 5.:5.:100.,  param = param, save_on = true, t_end = 1.0, stab_type = "full", save_interval = 0.1)
    @unpack Dθ, Nx, Nθ, χ, name,T,pert,k,δ,max_steps,δt = param
    λs = sqrt(Dθ)*Pes
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    
    if save_on
        try 
            stabdata = wload(filename)
        catch
        end
    end

    for ρ in ρs
        stable = []
        unstable = []
        unclassified = []
        data = Dict{String,Any}()
        #load ans
        for λ ∈ λs
            if (stab_type == "lin")&(λ ∉ stable)&(λ ∉ unstable)&(λ ∉ unclassified)
                try
                    local t_saves, fa_saves, fp_saves, stab_dsit0, stab_dsit1
                        param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = λ/sqrt(Dθ), χ = χ, T = T, 
                            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = max_steps,
                            pert = pert, k =k, δ = δ
                        )
                        t_saves, fa_saves, fp_saves = load_pdes_1d(param,0.01; save_interval = 0.001, start_time = 0.0)

                        stab_dsit0 = dist_from_unif_1d(param, fa_saves[1], fp_saves[1])
                        stab_dsit1 = dist_from_unif_1d(param, fa_saves[2], fp_saves[2])

                        if (stab_dsit0>stab_dsit1)&(λ ∉ stable)
                            push!(stable, λ)
                        elseif (λ ∉ unstable)
                            push!(unstable, λ)
                        end
                catch
                end
            elseif (stab_type == "full")&(λ ∉ stable)&(λ ∉ unstable)&(λ ∉ unclassified)
                try     
                        local t_saves, fa_saves, fp_saves, stab_dsit0, stab_dsit1, dist_saves, n, end_slope
                        param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = λ/sqrt(Dθ), χ = χ, T = T, 
                            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = max_steps,
                            pert = pert, k =k, δ = δ
                        )
                        t_saves, fa_saves, fp_saves = load_pdes_1d(param,t_end; save_interval = save_interval, start_time = 0.0)
                        dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)
                        n = length(t_saves)

                        stab_dsit0 = dist_saves[1]
                        stab_dsit1 = maximum(dist_saves)

                        end_slope = dist_saves[n] - dist_saves[n-1]

                        if (2*stab_dsit0>stab_dsit1)&(λ ∉ stable)&(end_slope ≤ 1e-8)
                            push!(stable, λ)
                        elseif (λ ∉ unstable)&(2*stab_dsit0<stab_dsit1)
                            push!(unstable, λ)
                        else (λ ∉ unclassified)&(λ ∉ stable)&(λ ∉ unstable)
                            push!(unclassified, λ)
                        end
                catch
                end
            elseif (stab_type == "asymptotic")&(λ ∉ stable)&(λ ∉ unstable)
                try     
                        param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = λ/sqrt(Dθ), χ = χ, T = T, 
                            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = max_steps,
                            pert = pert, k =k, δ = δ
                        )

                        t_saves, fa_saves, fp_saves = load_pdes_1d(param,t_end; save_interval = 0.01, start_time = (t_end-0.1))
                        dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)

                        stab_dsit0 = dist_from_unif_1d(param, fa_saves[1], fp_saves[1])
                        stab_dsit1 = maximum(dist_saves)

                        if (5*stab_dsit0>stab_dsit1)&(λ ∉ stable)
                            push!(stable, λ)
                        elseif (λ ∉ unstable)
                            push!(unstable, λ)
                        end
                catch
                end
            end
        end
        #fillout ans
        push!(stable,0.)
        @pack! data = stable, unstable, unclassified
        stabdata["ρ = $(ρ)"] = data
    end

    #fillout ans
    stable = λs
    unstable = []
    data = Dict{String,Any}()
    @pack! data = stable, unstable
    stabdata["ρ = $(0.)"] = data
    stabdata["ρ = $(1.)"] = data

    if save_on
        wsave(filename,stabdata)
    end
    return stabdata
end
function fillout_stab_data(;stabdata = Dict{String,Any}(), ρs = 0.4:0.05:1.0, Pes = 5.:5.:100.,  param = param, save_on = true, t_end = 1.0, stab_type = "full", save_interval = 0.1)
    @unpack Dθ, Nx, Nθ, χ, name,T,pert,k,δ,max_steps,save_interval,δt = param
    λs = sqrt(Dθ)*Pes
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    
    if save_on
        try 
            stabdata = wload(filename)
        catch
        end
    end

    for ρ in ρs
        try
            local maxstab, minunstab
            data =stabdata["ρ = $(ρ)"]
            @unpack stable, unstable, unclassified = data
            minunstab = minimum(unstable)
            maxstab = maximum(stable)
            for λ ∈ λs
                if (λ < maxstab)&(λ ∉ stable)
                    push!(stable, λ)
                elseif (λ > minunstab)&(λ ∉ unstable)
                    push!(unstable, λ)
                end
            end
            #fillout ans
            data = Dict{String,Any}()
            @pack! data = stable, unstable, unclassified
            stabdata["ρ = $(ρ)"] = data
        catch
        end
    end

    #fillout ans
    stable = λs
    unstable = []
    data = Dict{String,Any}()
    @pack! data = stable, unstable
    stabdata["ρ = $(0.)"] = data
    stabdata["ρ = $(1.)"] = data

    if save_on
        wsave(filename,stabdata)
    end
    return stabdata
end
function fillout_stab_data_horizontal(;stabdata = Dict{String,Any}(), ρs = 0.0:0.05:1.0, Pes = 5.:5.:100.,  param = param, save_on = true, t_end = 1.0, stab_type = "full", save_interval = 0.1)
    @unpack Dθ, Nx, Nθ, χ, name,T,pert,k,δ,max_steps,save_interval,δt = param
    λs = sqrt(Dθ)*Pes
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    
    if save_on
        try 
            stabdata = wload(filename)
        catch
        end
    end

    for ρ in ρs
        if ρ ≤ 0.45
        try
            local maxstab
            data =stabdata["ρ = 0.45"]
            @unpack stable, unstable, unclassified = data
            maxstab = maximum(stable)
            for λ ∈ λs
                if (λ ≤ maxstab)&(λ ∉ stable)
                    push!(stable, λ)
                end
            end
            #fillout ans
            data = Dict{String,Any}()
            @pack! data = stable, unstable, unclassified
            stabdata["ρ = $(ρ)"] = data
        catch
            local maxstab
            data =stabdata["ρ = 0.45"]
            @unpack stable, unstable, unclassified = data
            maxstab = maximum(stable)
            for λ ∈ λs
                if (λ ≤ maxstab)&(λ ∉ stable)
                    push!(stable, λ)
                end
            end
            #fillout ans
            data = Dict{String,Any}()
            @pack! data = stable, unstable, unclassified
            merge!(stabdata,Dict{String,Any}("ρ = $(ρ)"=>data))
        end
        end
    end

    if save_on
        wsave(filename,stabdata)
    end
    return stabdata
end
function find_extreme_params(param,stabdata; ρs = 0.4:0.05:1.0)
    @unpack k, name, D, δt, Nx, Nθ, S, E, Dθ, T, save_interval, max_steps, max_runs, pert, δ, Pe, Dx, γ, χ = param
    params = []
    for ρ in ρs
        try
            local maxstab, minunstab
            data =stabdata["ρ = $(ρ)"]
            @unpack stable, unstable, unclassified = data
            minunstab = minimum(unstable)/sqrt(Dθ)
            maxstab = maximum(stable)/sqrt(Dθ)

            param_min = pde_param_fraction(; name = name, 
                ρ = ρ, Pe = minunstab, χ = χ, T = T, 
                Dθ = Dθ, δt = δt, Nx = Nx,Nθ = Nθ,
                save_interval = save_interval, max_steps = max_steps,
                pert = pert, k =k, δ = δ,
            )
            param_max = pde_param_fraction(; name = name, 
            ρ = ρ, Pe = maxstab, χ = χ, T = T, 
            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ,
            save_interval = save_interval, max_steps = max_steps,
            pert = pert, k =k, δ = δ,
            )
            push!(params,param_min)
            push!(params,param_max)
            
            for λ in unclassified
                Pe = λ/sqrt(Dθ)
                param_unlcass = pde_param_fraction(; name = name, 
                ρ = ρ, Pe = Pe, χ = χ, T = T, 
                Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                save_interval = save_interval, max_steps = max_steps,
                pert = pert, k =k, δ = δ,
                )
                push!(params,param_unlcass)
            end
        catch
        end
    end
    return params
end
function plot_stab_frac(fig, ax, stabdata; ρs = 0.4:0.05:1.0 ,xs = collect(0.4:0.001:1.0), xtic = 0.4:0.2:1, ytic = 0:10:100, axlim = [0.4, 1., 0., 40.], param = param, χ = 0.5 )
    @unpack Dθ, Dx,ρp = param
    stable_points = []
    unstable_points = []
    unclass_points = []

    for ρ in ρs
        try
            @unpack stable, unstable, unclassified = stabdata["ρ = $(ρ)"]
            for λ ∈ stable
                    append!(stable_points,  [ρ; λ/sqrt(Dθ)])
                    
            end
            for λ ∈ unstable
                    append!(unstable_points, [ρ; λ/sqrt(Dθ)])
            end
            for λ ∈ unclassified
                #append!(unstable_points, [ρ; λ/sqrt(Dθ)])
                append!(unclass_points, [ρ; λ/sqrt(Dθ)])
        end
        catch
        end
    end             
    stable_points   = reshape(stable_points, (2,Int64(length(stable_points)/2)) )
    unstable_points = reshape(unstable_points, (2,Int64(length(unstable_points)/2)) )
    unclass_points = reshape(unclass_points, (2,Int64(length(unclass_points)/2)) )
    #=
    n = length(xs)
    m = length(ys)
    zs = zeros(m,n)
    for i in 1:m, j in 1:n
            zs[i,j] = lin_stab_line_fraction(xs[j],χ; Dx =Dx ,Pe = ys[i], Dθ = Dθ)
    end
    ax.contour(xs,ys,zs; levels = [0])
    =#

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

    ax.errorbar(unclass_points[1,:],unclass_points[2,:], 
    #markersize = 400/L, 
    fmt= "*", 
    color = "purple",
    alpha=0.8,
    label = L"\mathrm{Unclassified}",
    )


    ax.errorbar(stable_points[1,:],stable_points[2,:], 
    #markersize = 400/L, 
    fmt= "s", 
    color = "blue",
    alpha=0.8,
    label = L"\mathrm{Stable}",
    )

    ax.errorbar(unstable_points[1,:],unstable_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "red",
    alpha=0.8,
    label = L"\mathrm{Unstable}",
    )

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
end
function plot_imaginary_frac(fig, ax; ρs = 0.4:0.05:1.0 ,xs = collect(0.4:0.001:1.0), ys =collect(0.0:0.5:100.0), xtic = 0.4:0.2:1, ytic = 0:10:100, axlim = [0.4, 1., 0., 100.], param = param, χ = 0.5 )
    @unpack Dθ, Dx, ρp = param

    n = length(xs)
    m = length(ys)
    zs = zeros(m,n)

    X=[]
    Y=[]
    Xs = append!(append!(collect(0.401:0.001:0.99),collect(0.99:0.0001:0.9999)),collect(0.4:0.00001:0.401))
    for x in Xs
        try
            local f
            f(y) = lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ)
            Pe = find_zero(f, (0.,  500.))
            push!(Y,Pe)
            push!(X,x)
        catch
        end
    end
    ax.plot(X,Y,color = "black")

    for i in 1:m, j in 1:n
        zs[i,j] = abs(lin_imaginary_fraction(xs[j],χ; Dx =Dx ,Pe = ys[i], Dθ = Dθ, k =20))
    end

    colmap = PyPlot.plt.cm.viridis_r
    norm1 = matplotlib.colors.Normalize(vmin=0., vmax= maximum(zs) );
    ax.contourf(xs, ys, zs; levels = 200, norm = norm1, cmap = colmap )
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)

    #ax.xaxis.set_ticks(xtic)
    #ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    ax.set_xlabel("ρ")
    ax.set_ylabel("Pe")
    ax.set_title("ℓ = $(1/sqrt(Dθ)), χ = $(χ)")
end
#
function pde_density_hist(fig::Figure, ax::PyObject, param::Dict{String,Any}, fa, fp; r =5, bins = 3, smoothing = false)
    @unpack Nx, Nθ = param
    edges = collect(0.:(1/(bins)):1.0)
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    if smoothing 
        h = zeros(Nx,Nx)
        B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
        for x₁ in 1:Nx, x₂ in 1:Nx
            for e in B
                y₁, y₂  = ([x₁, x₂] + e +[Nx-1,Nx-1]) .% Nx + [1,1]
                h[x₁,x₂] += ρ[y₁,y₂]/(2*r+1)^2
            end
        end
        h = reshape(h, Nx*Nx)
    else
        h = reshape(ρ, Nx*Nx)
    end
    ax.hist(h; bins = edges, histtype = "step", density = true)
    ax.xaxis.set_ticks(0:0.25:1)
end
function pde_density_hist_av(fig::Figure, ax::PyObject, param::Dict{String,Any}, fa_saves, fp_saves; r =1, bins = 200, smoothing = false, add_label = false, name = "label")
    @unpack Nx, Nθ, Pe = param
    edges = collect(0.:(1/(bins)):1.0)
    N = length(fa_saves)
    H = []
    for i in 1:N
        fa = fa_saves[i]
        fp = fp_saves[i]
        ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
        if smoothing
            h = zeros(Nx,Nx)
            B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]
            for x₁ in 1:Nx, x₂ in 1:Nx
                for e in B
                    y₁, y₂  = ([x₁, x₂] + e +[Nx-1,Nx-1]) .% Nx + [1,1]
                    h[x₁,x₂] += ρ[y₁,y₂]/(2*r+1)^2
                end
            end
            h = reshape(h, Nx*Nx)
        else
            h = reshape(ρ, Nx*Nx)
        end
        append!(H,h)
    end
    if add_label
        ax.hist(H; bins = edges, histtype = "step", density = true, label = name)
    else
        ax.hist(H; bins = edges, histtype = "step", density = true)
    end
    ax.xaxis.set_ticks(0:0.2:1)
end
function pde_density_hist_1d(fig::Figure, ax::PyObject, param::Dict{String,Any}, fa, fp; r = 5, bins = 3)
    @unpack Nx, Nθ = param
    edges = collect(0.:(1/(bins)):1.0)
    ρ = fp + sum(fa; dims =2)[:,1].*(2*π/Nθ)
    h = zeros(Nx)
    B = [i for i in (-r):1:r ]
    for x₁ in 1:Nx
        for e in B
            y₁  = (x₁ + e +Nx-1) % Nx + 1
            h[x₁] += ρ[y₁]/(2*r+1)
        end
    end
    ax.hist(h; bins = edges, histtype = "step", density = true)
    ax.xaxis.set_ticks(0:0.25:1)
end
function pde_density_hist_av_save(param::Dict{String,Any}, fa_saves, fp_saves; r =1, smoothing = false)
    @unpack Nx, Nθ, Pe, name, ρa, ρp, λ, δt, Dθ = param
    N = length(fa_saves)
    H = []
    for i in 1:N
        fa = fa_saves[i]
        fp = fp_saves[i]
        ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
        if smoothing
            h = zeros(Nx,Nx)
            B = [[i,j] for i in (-r):1:r for j in (-r):1:r ]#if (i^2+j^2) ≤ r^2]
            n = length(B)
            for x₁ in 1:Nx, x₂ in 1:Nx
                for e in B
                    y₁, y₂  = ([x₁, x₂] + e +[Nx-1,Nx-1]) .% Nx + [1,1]
                    h[x₁,x₂] += ρ[y₁,y₂]/n
                end
            end
            h = reshape(h, Nx*Nx)
        else
            h = reshape(ρ, Nx*Nx)
        end
        append!(H,h)
    end
    data = Dict{String,Any}()
    @pack! data = H
    @unpack Nx, Nθ, Pe, name, ρa, ρp, λ, δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/av_hist/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2"
    wsave(filename,data)
end
function pde_density_hist_av_load(fig::Figure, ax::PyObject, param::Dict{String,Any}; r =1, bins = 200, smoothing = false, add_label = false, label_name = "label", linewidth = 1.2)
    @unpack Nx, Nθ, Pe = param
    @unpack Nx, Nθ, Pe, name, ρa, ρp, λ, δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/av_hist/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2"
    data = wload(filename)
    @unpack H = data
    edges = collect((-1/(2*bins)):(1/(bins)):(1.0+(1/(2*bins))))
    if add_label
        ax.hist(H; bins = edges, histtype = "step", density = true, label = label_name, linewidth=linewidth)
    else
        ax.hist(H; bins = edges, histtype = "step", density = true)
    end
    ax.xaxis.set_ticks(0:0.2:1)
end

function time_density_hist_save(param::Dict{String,Any}, t_saves, η_saves; r = 3, bins = 3, name= "label", add_label = false)
    h = [];
    for η ∈ η_saves
        append!(h, density_hist(param, η; r = r));
    end
    data = Dict{String,Any}()
    @pack! data = h
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/hist_save/r=$(r)_size=$(L)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    wsave(filename,data)
end

function time_density_hist_time_save(param::Dict{String,Any}, t, η; r = 3, bins = 3, name= "label", add_label = false)
    h = [];
    append!(h, density_hist(param, η; r = r));
    data = Dict{String,Any}()
    @pack! data = h
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/hist_save/r=$(r)_size=$(L)/t=$(round(t;digits=2))_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    wsave(filename,data)
end

function time_density_polarisation_save(param::Dict{String,Any}, t, η; r = 3, bins = 3, name= "label", add_label = false)
    h = [];
    append!(h, pol_hist(param, η; r = r));
    data = Dict{String,Any}()
    @pack! data = h
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/pol_save/r=$(r)_size=$(L)/t=$(round(t;digits=2))_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    wsave(filename,data)
end

function time_density_hist_load(fig::Figure, ax::PyObject,param::Dict{String,Any}; r = 3, bins = 3, label_name= "label", add_label = false, ignore_zero = false)
    @unpack name, L, λ, γ, ρa, ρp, Δt, Dθ = param
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/sims_pro/hist_save/r=$(r)_size=$(L)/active=$(ρa)_passive=$(ρp)_lamb=$(λ)_gamma=$(γ)_Δt=$(Δt)_Dθ=$(Dθ).jld2"
    data = wload(filename)
    @unpack h = data
    edges = []
    if ignore_zero
        edges = collect((1/(2*bins)):(1/(bins)):(1+(1/(2*bins))))
    else
        edges = collect((-1/(2*bins)):(1/(bins)):(1+(1/(2*bins))))
    end
    if add_label
        ax.hist(h; bins = edges, histtype = "step", density = true, label = label_name)
    else
        ax.hist(h; bins = edges, histtype = "step", density = true)
    end
end
###

function find_stab_data_ap(;stabdata = Dict{String,Any}(), ρs = 0.4:0.05:1.0, Pes = 5.:5.:100.,  param = param, save_on = true, t_end = 1.0, t_start = 0.0, stab_type = "full", save_interval = 0.1)
    @unpack Dθ, Nx, Nθ, χ, name,T,pert,k,δ,max_steps,δt = param
    λs = sqrt(Dθ)*Pes
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_type=$(stab_type)_Nx=$(Nx)_Nθ=$(Nθ)_Dθ=$(Dθ)_χ=$(χ).jld2"
    
    if save_on
        try 
            stabdata = wload(filename)
        catch
        end
    end

    for ρ in ρs
        stable = []
        unstable = []
        unclassified = []
        wave =[]
        data = Dict{String,Any}()
        #load ans
        for λ ∈ λs
            if (stab_type == "full")&(λ ∉ stable)&(λ ∉ unstable)&(λ ∉ unclassified)&(λ ∉ wave)
                local t_saves, fa_saves, fp_saves, stab_dsit0, stab_dsit1, dist_saves, n, end_slope
                try     
                        Pe = λ/sqrt(Dθ)
                        param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = Pe, χ = χ, T = T, 
                            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = max_steps,
                            pert = pert, k =k, δ = δ
                        )
                        t_saves, fa_saves, fp_saves = load_pdes_1d(param,t_end; save_interval = save_interval, start_time = t_start)

            
                        dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)
                        n = length(t_saves)

                        stab_dsit0 = δ
                        stab_dsit1 = maximum(dist_saves)

                        end_slope = dist_saves[n] - dist_saves[n-1]

                        # condition for travelling wave 

                        Δfa_Δt = fa_saves[n]-fa_saves[n-1]
                        Δfp_Δt = fp_saves[n]-fp_saves[n-1]

                        L2_dist = sqrt( 2*π*sum( (Δfa_Δt ).^2)/(Nx*Nθ) + sum( (Δfp_Δt).^2)/(Nx) )

                        speed = L2_dist


                        if (2*stab_dsit0>stab_dsit1)&(λ ∉ stable)&(end_slope ≤ 1e-8)
                            push!(stable, λ)
                        elseif (λ ∉ unstable)&(2*stab_dsit0<stab_dsit1)&(speed<0.01)
                            push!(unstable, λ)
                        elseif (λ ∉ unstable)&(2*stab_dsit0<stab_dsit1)&(speed≥0.01)
                            push!(wave, λ)
                        elseif (λ ∉ unclassified)&(λ ∉ stable)&(λ ∉ unstable)&(λ ∉ wave)
                            push!(unclassified, λ)
                        end

                catch
                    println("load error")
                    println("ρ = $(ρ) χ = $(χ) Pe = $(Pe)")
                    println(t_saves)
                end
            end
        end
        #fillout ans
        # push!(stable,0.)
        @pack! data = stable, unstable, unclassified, wave
        stabdata["ρ = $(ρ)"] = data
    end

    # #fillout ans
    # stable = λs
    # unstable = []
    # data = Dict{String,Any}()
    # @pack! data = stable, unstable
    # stabdata["ρ = $(0.)"] = data
    # stabdata["ρ = $(1.)"] = data

    if save_on
        wsave(filename,stabdata)
    end
    return stabdata
end

function plot_stab_frac_ap(fig, ax, stabdata; ρs = 0.4:0.05:1.0 ,xs = collect(0.4:0.001:1.0), xtic = 0.4:0.2:1, ytic = 0:10:100, axlim = [0.4, 1., 0., 40.], param = param, χ = 0.5 )
    @unpack Dθ, Dx,ρp = param
    stable_points = []
    unstable_points = []
    wave_points = []
    unclass_points = []

    for ρ in ρs
        try
            @unpack stable, unstable, unclassified, wave = stabdata["ρ = $(ρ)"]
            for λ ∈ stable
                    append!(stable_points,  [ρ; λ/sqrt(Dθ)])
                    
            end
            for λ ∈ unstable
                    append!(unstable_points, [ρ; λ/sqrt(Dθ)])
            end
            for λ ∈ wave
                append!(wave_points, [ρ; λ/sqrt(Dθ)])
            end
            for λ ∈ unclassified
                #append!(unstable_points, [ρ; λ/sqrt(Dθ)])
                append!(unclass_points, [ρ; λ/sqrt(Dθ)])
            end
        catch
        end
    end             
    stable_points   = reshape(stable_points,    (2,Int64(length(stable_points)/2))      )
    unstable_points = reshape(unstable_points,  (2,Int64(length(unstable_points)/2))    )
    wave_points     = reshape(wave_points,      (2,Int64(length(wave_points)/2))        )
    unclass_points  = reshape(unclass_points,   (2,Int64(length(unclass_points)/2))     )
    #=
    n = length(xs)
    m = length(ys)
    zs = zeros(m,n)
    for i in 1:m, j in 1:n
            zs[i,j] = lin_stab_line_fraction(xs[j],χ; Dx =Dx ,Pe = ys[i], Dθ = Dθ)
    end
    ax.contour(xs,ys,zs; levels = [0])
    =#

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

    ax.errorbar(unclass_points[1,:],unclass_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "purple",
    alpha=0.8,
    label = L"\mathrm{Unclassified}",
    )


    ax.errorbar(stable_points[1,:],stable_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "blue",
    alpha=0.8,
    label = L"\mathrm{Stable}",
    )

    ax.errorbar(unstable_points[1,:],unstable_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "red",
    alpha=0.8,
    label = L"\mathrm{Unstable}",
    )

    ax.errorbar(wave_points[1,:],wave_points[2,:], 
    #markersize = 400/L, 
    fmt= "o", 
    color = "green",
    alpha=0.8,
    label = L"\mathrm{Wave}",
    )

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
end