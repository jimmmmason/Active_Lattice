cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using DrWatson, KernelDensity, Peaks

## running simulation
    function new_param(DT::Float64, v0::Float64, DR::Float64, N::Int64, Δx::Float64, Lx::Float64, Ly::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, name::String = "test", pert::String = "lin", save_interval::Float64 = 0.001, save_on::Bool = false)
        param::Dict{String, Any} = Dict{String,Any}()
        N₁::Int64,N₂::Int64 = Int64(Lx*N ÷ 1), Int64(Ly*N ÷ 1)
        Nx::Int64 = Int64(Lx/Δx ÷ 1)
        @pack! param = DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ, T, name, Nx, N₁, N₂, save_interval, save_on, pert
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

    function find_xpeak_ft(ts , ft; time_length = 0.1)
        index_end   = length(ts)
        indebottom_gap_2 = index_end - length([ t for t in ts if t>ts[end]-time_length])
        
        av_rho = sum( ft[indebottom_gap_2:1:index_end,:,:]; dims = (1,3) )[1,:,1]
        N = length(av_rho)
        M = 3*N
        av_rho = reduce(vcat,[av_rho,av_rho,av_rho])
    
        smooth_rho = KernelDensitySJ.smooth(collect(1.0:M),av_rho,Int64(N/10 ÷1), collect(1.0:M))
        pks, vals = findmaxima(smooth_rho)
        pks = [x for x in pks if (x>M/6)&(x<5*M/6)]
        pks, proms = peakproms(pks, smooth_rho)
        pk = pks[argmax(proms)]
        return pk
    end

    function t_dff(ts , ft; N=100, gap = 1)
        ts = ts[gap:gap:end]
        Nt = length(ts)

        f_dt = zeros(Nt)

        for i in 2:Nt
            local fdiff, tdiff
            tdiff = ts[i] - ts[i-1]
            fdiff = ft[i*gap,:,:] - ft[(i-1)*gap,:,:]
            f_dt[i]    = norm(fdiff /tdiff )/sqrt(N)
        end
        return ts, f_dt
    end
#

println("v2.1")