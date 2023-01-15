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
function pde_param_fraction(; name = "test", D =1., Dx = 1., Pe =1., Dθ = 10, ρ= 0.5, χ = 1.0, Nx = 100, Nθ = 20, δt = 1e-5, T= 0.001, save_interval = 0.01, max_steps = 1e8, max_runs = 6, λ_step = 10., λmax = 100., λs = 20.:20.:100., pert = "n=1", δ = 0.01, k=20)
    Ω  = [[i,j] for i in 1:Nx for j in 1:Nx ] 
    S  = [ θ for θ in 1:Nθ]
    E = [[1,0],[0,1],[0,-1],[-1,0],]
    ρp = 0.0
    λ = Pe*sqrt(Dθ)
    ρp = (1-χ)*ρ
    ρa = χ*ρ
    param = Dict{String,Any}()
    @pack! param = k, name, D, λ, ρa, ρp, δt, Nx, Nθ, S,  E, Dθ, T, save_interval, max_steps, max_runs, λ_step, λmax, λs, pert, δ, Pe, Dx, χ, ρ
    return param
end
function sim_param_fraction(;  name = "test", D =1. , Pe =1. ,ρ = 0.5, χ = 0.5, L=10, d=2, Δt = 0.01, Dθ =10., T=1.0, γ = 0.)
    param = Dict{String,Any}()
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    λ = Pe*sqrt(Dθ)
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

    @pack! param = name, L, D, λ, γ, ρa, ρp, Δt, E, site_distribution, angles, rates, Dθ, T, POSITIONS
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
function lin_imaginary_fraction(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40 )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end
#
function find_stab_data(;stabdata = Dict{String,Any}(), ρs = 0.4:0.05:1.0, Pes = 5.:5.:100.,  param = param, save_on = true, t_end = 1.0, stab_type = "full")
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
        stable = []
        unstable = []
        data = Dict{String,Any}()
        #load ans
        for λ ∈ λs
            if (stab_type == "lin")&(λ ∉ stable)&(λ ∉ unstable)
                try
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
            elseif (stab_type == "full")&(λ ∉ stable)&(λ ∉ unstable)
                try
                        param = pde_param_fraction(; name = name, 
                            ρ = ρ, Pe = λ/sqrt(Dθ), χ = χ, T = T, 
                            Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, 
                            save_interval = save_interval, max_steps = max_steps,
                            pert = pert, k =k, δ = δ
                        )
                        t_saves, fa_saves, fp_saves = load_pdes_1d(param,t_end; save_interval = 0.01, start_time = 0.0)
                        dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)

                        stab_dsit0 = dist_from_unif_1d(param, fa_saves[1], fp_saves[1])
                        stab_dsit1 = maximum(dist_saves)

                        if (2*stab_dsit0>stab_dsit1)&(λ ∉ stable)
                            push!(stable, λ)
                        elseif (λ ∉ unstable)
                            push!(unstable, λ)
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
        @pack! data = stable, unstable
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
end
function plot_stab_frac(fig, ax, stabdata; ρs = 0.4:0.05:1.0 ,xs = collect(0.4:0.001:1.0), xtic = 0.4:0.2:1, ytic = 0:10:100, axlim = [0.4, 1., 0., 40.], param = param, χ = 0.5 )
    @unpack Dθ, Dx,ρp = param
    stable_points = []
    unstable_points = []

    for ρ in ρs
        try
            @unpack stable, unstable = stabdata["ρ = $(ρ)"]
            for λ ∈ stable
                    append!(stable_points,  [ρ; λ/sqrt(Dθ)])
                    
            end
            for λ ∈ unstable
                    append!(unstable_points, [ρ; λ/sqrt(Dθ)])
            end
        catch
        end
    end             
    stable_points   = reshape(stable_points, (2,Int64(length(stable_points)/2)) )
    unstable_points = reshape(unstable_points, (2,Int64(length(unstable_points)/2)) )
    
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
    ax.plot(X,Y,color = "black")

    
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

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    ax.set_xlabel("ρ")
    ax.set_ylabel("Pe")
    ax.set_title("ℓ = $(1/sqrt(Dθ)), χ = $(χ)")
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
            Pe = find_zero(f, (0.,  100.))
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

    ax.xaxis.set_ticks(xtic)
    ax.yaxis.set_ticks(ytic)
    ax.axis(axlim)
    ax.set_xlabel("ρ")
    ax.set_ylabel("Pe")
    ax.set_title("ℓ = $(1/Dθ), χ = $(χ)")
end
#
function pde_density_hist(fig::Figure, ax::PyObject, param::Dict{String,Any}, fa, fp; bins = 3)
    @unpack Nx, Nθ = param
    edges = collect((-1/(2*bins)):(1/(bins)):(1+1/(2*bins)))
    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
    h = reshape(ρ, Nx*Nx)
    ax.hist(h; bins = edges, histtype = "step", density = true)
    ax.xaxis.set_ticks(0:0.25:1)
end
