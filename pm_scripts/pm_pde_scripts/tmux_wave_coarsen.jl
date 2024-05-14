cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
# Load relevant fuctions
# include("/home/jm2386/Active_Lattice/src/pm_pde_functions.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");


function perturb_safe_pde!(f::Matrix{Float64}, param::Dict{String, Any}; wave_num = 1)
    @unpack DT, v0, DR, Δx, Lx, ϕa, ϕp, T , name, Nx, save_interval, save_on, δt, δ, pert = param
    
    if pert == "rand"
        pertf = 2*rand(Nx,3) .-1
    elseif pert == "double"
        ω, value, vector = lin_pert_values(param;wave_num = wave_num,wave_choice=3)

        wave   = exp.((1:Nx)*(wave_num*im*2*π/Nx))
        pertf  = zeros(Nx,3)

        pertf[:,1] = real.( wave*(vector[2]- vector[3])/2 )
        pertf[:,2] = real.( wave*(vector[2]+ vector[3])/2 ) 
        pertf[:,3] = real.( wave*(vector[1]-vector[2]) )

        # ω, value, vector = lin_pert_values(param;wave_num = wave_num,wave_choice=2)

        # wave   = exp.(-(1:Nx)*(wave_num*im*2*π/Nx))
        # pertf  = zeros(Nx,3)

        # pertf[:,1] = real.( wave*(vector[1]- vector[2])/2 )
        # pertf[:,2] = real.( wave*(vector[1]+ vector[2])/2 ) 
        # pertf[:,3] = real.( wave*(vector[3]) )
        
        pertf += pertf[end:-1:1,[2,1,3]]

        print("sym pert: $(node_sym(pertf))")
    else
        ω, value, vector = lin_pert_values(param;wave_num = wave_num)

        wave   = exp.((1:Nx)*(wave_num*im*2*π/Nx))
        pertf    = zeros(Nx,3)

        pertf[:,1] = real.( wave*(vector[2]- vector[3])/2 )
        pertf[:,2] = real.( wave*(vector[2]+ vector[3])/2 ) 
        pertf[:,3] = real.( wave*(vector[1]-vector[2]) )

        println("max pert")
    end

    c = norm(pertf)/sqrt(Nx)

    max_pert = maximum(δ*sum(pertf;dims=2)/c)
    max_f = maximum(sum(f;dims=2))
    min_f = minimum(sum(f;dims=2))
    ϵ = min(1e-5,(1-max_f)/2)
    if max_pert>(1-(max_f))
        sf = (1-(max_f)-ϵ)/max_pert
    elseif max_pert>(min_f)
        sf = ((min_f-ϵ))/max_pert
    else
        sf = 1
    end

    f += sf*δ*pertf/c

    return f
end

# Load initial wave
    param = get_grid_param(12,4)
    param["ϕa"] = 0.5
    param["ϕp"] = 0.1
    param["Lx"] = 20.0
    param["v0"] = 20.0
    param["δ"] = 0.01
    param["Δx"] = param["Lx"]/param["Nx"]
    param["save_interval"] = 0.01
    param["T"] = 40.0
    t_saves, f_saves = load_compress_pde(param);
#

f = f_saves[end]
param["Δx"] = param["Lx"]/param["Nx"]
param["δ"] = 0.01
param["T"] = 40.0
println(maximum(sum(f;dims=2)), minimum(sum(f;dims=2)))
f = perturb_safe_pde!(f, param; wave_num = 1);
println(maximum(sum(f;dims=2)), minimum(sum(f;dims=2)))

param["name"] = "coarsening_test"
run_current_pde(param,40.0, f,0.0);
ts,fs = load_compress_pde(param);


# export PATH="/home/jm2386/.local/ /bin:${PATH}"
# nice -n 19 julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/tmux_wave_coarsen.jl