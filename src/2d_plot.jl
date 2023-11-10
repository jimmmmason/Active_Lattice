cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using DrWatson, KernelDensity, Peaks

## running simulation 
    function _2d_new_param(DT::Float64, v0::Float64, DR::Float64, N::Int64, Nθ::Int64, Δx::Float64, Lx::Float64, Ly::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, Δt::Float64 = 0.01, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
        Dθ  = DR
        D   = DT/Lx^2
        λ  = v0/Lx
        N₁  = Int64(Lx*N ÷ 1)
        N₂  = Int64(Lx*N ÷ 1)
        Nx  = Int64((Lx/Δx) ÷ 1)
        param = uniform_initial_param(; name = name, D =D , λ =λ ,ρa = ϕa, ρp = ϕp, L=N₁, Δt = Δt, Dθ = Dθ, T=T)
        @pack! param = DT, v0, DR, N, Nθ, Δx, Lx, Ly, ϕa, ϕp, δt, δ, T, name, Nx, N₁, N₂, save_interval, save_on
        return param
    end
#

## proccessing funcitons 
    function rho_to_rgb(f)
        Nx, Ny, k = size(f)
        rgb_image = ones(Ny,Nx,3)
    
        rgb_image[:,:,3] = -(f[:,Ny:-1:1,1]' +f[:,Ny:-1:1,2]').^2 .+1
        rgb_image[:,:,1] = -(f[:,Ny:-1:1,3]') .+1
        rgb_image[:,:,2] = -(f[:,Ny:-1:1,1]' + f[:,Ny:-1:1,2]'  -f[:,Ny:-1:1,3]').^2  .+1
    
        return rgb_image
    end
#

## plotting funcitons
    function plot_eta!(ax, param::Dict{String,Any}, t::Float64, η::Array{Float64,3}; title = true, ϵ = 0.1)
        @unpack name, N, λ, ρa, ρp, Lx, N₁ = param
        ax.clear()
        passive, active, directions = extract_points(η; N₁=N₁, N=N, ϵ=ϵ)
        dx = cos.(directions)
        dy = sin.(directions)
        t = round(t; digits=2)
        loc_den = local_density(η; N₁=N₁, N=N, ϵ=ϵ)
        # Plot points
        loc_den = loc_den[:,N₁:-1:1]
        xs = []
        ys = []
        for x₁ in 1:N₁, x₂ in 1:N₁
            push!(xs, (x₁-0.5)/N)
            push!(ys, (x₂-0.5)/N)
        end
        colmap = bone()
        norm1 = matplotlib.colors.Normalize(vmin= 0, vmax= 1) 
        im1 = ax.matshow(-loc_den'.+1; norm = norm1,  cmap = colmap, extent = [0,Lx,0,Lx], alpha = 0.5)
            
        ax.quiver(active[1,:],active[2,:], dx,dy, cmap = hsv(),
                directions,
                scale_units = "x",
                pivot = "mid",
                minlength = 0.1,
                minshaft = 1,
                width =1/N,
                headlength = 5,
                headaxislength=5,
                scale = N,
            )
        ax.errorbar(passive[1,:],passive[2,:], 
            markersize = 400/N, 
            fmt= "o", 
            color = "black",
            alpha=0.8,
        )

        #figure configuration
        ax.xaxis.set_ticks(0.:1.0:Lx)
            ax.yaxis.set_ticks(0.:1.0:Lx)
            ax.tick_params()
            ax.axis([0., Lx, 0., Lx])
            ax.set_aspect("equal")
            if title
                @unpack v0, ϕa, ϕp, Lx = param
                latex_title = latexstring("\$ \\phi = $(ϕa+ϕp), \\mathrm{v0} = $(v0), L_x = $(Lx), t = $(t)\$")
                ax.set_title(latex_title,fontsize=15)
            else
                #ax.set_title("Φ = $(Φ)")
            end
            ax.tick_params(direction = "in")
            ax.xaxis.tick_bottom()
            ax.set_xlabel(L"x",fontsize = 15)
            ax.set_ylabel(L"y",fontsize = 15)
    end 

    function plot_pde_rho(ax, cbar_ax, param::Dict{String,Any}, f; cmin = -1, cmax = 1, cbar = true)
        @unpack Nx, Nθ, ρa, ρp, Lx= param
        #collect data
            #m = reshape(mag_f(f; Nx=Nx, Nθ=Nθ),2,Nx,Nx)
            ρ = ρ_f(f; Nx=Nx, Nθ=Nθ)
            #absmag  = sqrt.(m[2,:,:].^2+m[1,:,:].^2)
            dx = 1/Nx
        #plot data
            colmap = PyPlot.plt.cm.magma
            norm1 = matplotlib.colors.Normalize(vmin= cmin, vmax= cmax) 
            im1 = ax.matshow(ρ'; norm = norm1,  cmap = colmap, extent = [0,Lx,0,Lx], interpolation = "bilinear")

            if cbar
                fig.colorbar(im1, norm=norm1, cmap = colmap, cax = cbar_ax)
                cbar_ax.set_ylabel(L"\rho",fontsize = 15)
                cbar_ax.tick_params(direction = "in")
                cbar_ax.yaxis.set_ticks(cmin:0.2:cmax)
            end
        #figure configuration
            ax.xaxis.set_ticks(0:1:Lx)
            ax.yaxis.set_ticks(0:1:Lx)
            ax.tick_params(direction = "in",labelbottom = true, labeltop = false)
            ax.axis([0, Lx, 0, Lx])
            ax.set_aspect("equal")
            ax.set_xlabel(L"x",fontsize = 15)
            ax.set_ylabel(L"y",fontsize = 15)
    end

    function plot_pde_mag(ax, cbar_ax, param::Dict{String,Any}, f::Array{Float64,3}; cmin = -1, cmax = 1, cbar = true, density = 1.0)
        @unpack Nx, Nθ, ρa, ρp, Lx= param
        #collect data
            m = reshape(mag_f(f; Nx=Nx, Nθ=Nθ),2,Nx,Nx)
            #ρ = ρ_f(f; Nx=Nx, Nθ=Nθ)
            absmag  = sqrt.(m[2,:,:].^2+m[1,:,:].^2)
            dx = Lx/Nx
        #figure configuration
            ax.xaxis.set_ticks(0:1:Lx)
            ax.yaxis.set_ticks(0:1:Lx)
            ax.tick_params(direction = "in", labelbottom = true)
            ax.axis([0, Lx, 0, Lx])
            ax.set_aspect("equal")
            ax.set_xlabel(L"x",fontsize = 15)
            ax.set_ylabel(L"y",fontsize = 15)
        #grid points
            x = dx:dx:Lx
            y = dx:dx:Lx
            xx = [x̃ for x̃ ∈ x, ỹ ∈ y]'
            yy = [ỹ for x̃ ∈ x, ỹ ∈ y]'
        #plot data
            colmap = PyPlot.plt.cm.viridis
            norm1 = matplotlib.colors.Normalize(vmin=cmin, vmax= cmax);
            im1 = ax.streamplot(xx, yy, m[1,:,:]', m[2,:,:]', color = absmag', cmap = colmap, norm = norm1, density = density)#2.5
            
            if cbar
                fig.colorbar(im1.lines, cax = cbar_ax)
                cbar_ax.set_ylabel(L"{\bf m}",fontsize = 15)
                cbar_ax.tick_params(direction = "in")
                cbar_ax.yaxis.set_ticks(cmin:0.1:cmax)
            end
    end

    function plot_sim_rho(ax, cbar_ax, param::Dict{String,Any}, η::Array{Float64,3}; cmin = -1, cmax = 1, cbar = true, ϵ = 0.1)
        @unpack N₁, ρa, ρp, Lx, N= param
        #collect data
            # imag_m = local_polarisation(η; N₁ = N₁, N = N, ϵ = ϵ)
            # m = zeros(2,N₁,N₁)
            # m[1,:,:] = real.(imag_m)
            # m[2,:,:] = imag.(imag_m)
            ρ = local_density(η; N₁ = N₁, N = N, ϵ = ϵ)
            #absmag  = sqrt.(m[2,:,:].^2+m[1,:,:].^2)
            dx = Lx/N₁
        #plot data
            colmap = PyPlot.plt.cm.magma
            norm1 = matplotlib.colors.Normalize(vmin= cmin, vmax= cmax) 
            im1 = ax.matshow(ρ'; norm = norm1,  cmap = colmap, extent = [0,Lx,0,Lx], interpolation = "bilinear")

            if cbar
                fig.colorbar(im1, norm=norm1, cmap = colmap, cax = cbar_ax)
                cbar_ax.set_ylabel(L"\rho",fontsize = 15)
                cbar_ax.tick_params(direction = "in")
                cbar_ax.yaxis.set_ticks(cmin:0.2:cmax)
            end
        #figure configuration
            ax.xaxis.set_ticks(0:1:Lx)
            ax.yaxis.set_ticks(0:1:Lx)
            ax.tick_params(direction = "in",labelbottom = true, labeltop = false)
            ax.axis([0, Lx, 0, Lx])
            ax.set_aspect("equal")
            ax.set_xlabel(L"x",fontsize = 15)
            ax.set_ylabel(L"y",fontsize = 15)
    end

    function plot_sim_mag(ax, cbar_ax, param::Dict{String,Any}, η; cmin = -1, cmax = 1, cbar = true, density = 1.0, ϵ = 0.1)
        @unpack N₁, ρa, ρp, Lx, N= param
        #collect data
            imag_m = local_polarisation(η; N₁ = N₁, N = N, ϵ = ϵ)
            m = zeros(2,N₁,N₁)
            m[1,:,:] = real.(imag_m)
            m[2,:,:] = imag.(imag_m)
            #ρ = local_density(η; N₁ = N₁, N = N, ϵ = ϵ)
            absmag  = sqrt.(m[2,:,:].^2+m[1,:,:].^2)
            dx = Lx/N₁
        #figure configuration
            ax.xaxis.set_ticks(0:1:Lx)
            ax.yaxis.set_ticks(0:1:Lx)
            ax.tick_params(direction = "in", labelbottom = true)
            ax.axis([0, Lx, 0, Lx])
            ax.set_aspect("equal")
            ax.set_xlabel(L"x",fontsize = 15)
            ax.set_ylabel(L"y",fontsize = 15)
        #grid points
            x = dx:dx:Lx
            y = dx:dx:Lx
            xx = [x̃ for x̃ ∈ x, ỹ ∈ y]'
            yy = [ỹ for x̃ ∈ x, ỹ ∈ y]'
        #plot data
            colmap = PyPlot.plt.cm.viridis
            norm1 = matplotlib.colors.Normalize(vmin=cmin, vmax= cmax);
            im1 = ax.streamplot(xx, yy, m[1,:,:]', m[2,:,:]', color = absmag', cmap = colmap, norm = norm1, density = density)#2.5
            
            if cbar
                fig.colorbar(im1.lines, cax = cbar_ax)
                cbar_ax.set_ylabel(L"{\bf m}",fontsize = 15)
                cbar_ax.tick_params(direction = "in")
                cbar_ax.yaxis.set_ticks(cmin:0.1:cmax)
            end
    end
#


println("v2.1")