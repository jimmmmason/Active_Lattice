cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using DrWatson

## running simulation 
    function new_param(DT::Float64, v0::Float64, DR::Float64, N::Int64, Δx::Float64, Lx::Float64, Ly::Float64, ϕa::Float64, ϕp::Float64, δt::Float64, δ::Float64; T::Float64 = 0.001, name::String = "test", save_interval::Float64 = 0.001, save_on::Bool = false)
        param::Dict{String, Any} = Dict{String,Any}()
        N₁::Int64,N₂::Int64 = Int64(Lx*N ÷ 1), Int64(Ly*N ÷ 1)
        Nx::Int64 = Int64(Lx/Δx ÷ 1)
        @pack! param = DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ, T, name, Nx, N₁, N₂, save_interval, save_on
        return param
    end
#

## proccessing funcitons 
    function rho_to_rgb(f)
        Nx, Ny, k = size(f)
        rgb_image = ones(Ny,Nx,3)
    
        rgb_image[:,:,3] = -f[:,Ny:-1:1,1]' -f[:,Ny:-1:1,2]' .+1
        rgb_image[:,:,1] = -f[:,Ny:-1:1,3]' .+1
        rgb_image[:,:,2] = -f[:,Ny:-1:1,1]' -f[:,Ny:-1:1,2]'  -f[:,Ny:-1:1,3]'  .+1
    
        return rgb_image
    end
#

println("v1.0")