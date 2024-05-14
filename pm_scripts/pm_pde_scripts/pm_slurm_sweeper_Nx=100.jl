#Set up
cd("/home/jm2386/Active_Lattice/");
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_sims.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");

# Access the command-line argument 'i'
input = parse(Int64, ARGS[1]);


# # v push 
# ϕa = d2.(collect(0.2:0.01:0.75))[input] #56
# ϕp = 0.25

# Lx, Nx, ϕa, ϕp, v0 = 100.0, 3200, ϕa, ϕp, 7.5
# # ϕp_sweep = d2.(collect(0.3:(-0.01):0.0))
# ϕp_sweep = d2.(collect(0.3:(0.01):0.4))

# f,u,c = load_full(Lx, Nx, 0.5, 0.3, v0)
# no_load = true
# j = length(ϕp_sweep)
# ϵ = 1e-4
# while (no_load)&(j>0)
#     global f,u,c,ϕp,Lx,Nx,ϕa,ϕp,v0,no_load,j,dj
#     try
#         ϕp = ϕp_sweep[j]
#         f,u,c = load_full(Lx, Nx, ϕa, ϕp, v0)
#         normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#         if normf > ϵ
#             no_load = false
#         else
#             j -= 1
#         end
#     catch
#         j -=1
#     end
# end

# # ϕp_sweep = d2.(collect(ϕp:(-0.01):0.0))
# ϕp_sweep = d2.(collect(ϕp:(0.01):0.4))
# ϕa_sweep = fill(ϕa,length(ϕp_sweep))

# sweep = (Lx, Nx, ϕa, ϕp, v0)
# param_sweep = [(Lx, Nx, d2(ϕa), d2(ϕp), v0) for (ϕa,ϕp) in zip(ϕa_sweep,ϕp_sweep)];

# h push 
# ϕa = 0.27
ϕa = 0.27
ϕp = d2.(collect(0.4:-0.01:0.0))[input] #41

Lx, Nx, ϕa, ϕp, v0 = 100.0, 3200, ϕa, ϕp, 7.5
ϕa_sweep = d2.(collect(ϕa:(-0.01):0.0))
# ϕa_sweep = d2.(collect(ϕa:(0.01):(1-ϕp)))


f,u,c = load_full(Lx, Nx, 0.5, 0.3, v0)
no_load = true
j = length(ϕa_sweep)
ϵ = 1e-4
while (no_load)&(j>0)
    global f,u,c,ϕp,Lx,Nx,ϕa,ϕp,v0,no_load,j,dj
    try
        ϕa = ϕa_sweep[j]
        f,u,c = load_full(Lx, Nx, ϕa, ϕp, v0)
        normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
        if normf > ϵ
            no_load = false
        else
            j -= 1
        end
    catch
        j -=1
    end
end

ϕa_sweep = d2.(collect(ϕa:-0.01:0.0))
# ϕa_sweep = d2.(collect(ϕa:0.01:(1-ϕp)))
ϕp_sweep = fill(ϕp,length(ϕa_sweep))

sweep = (Lx, Nx, ϕa, ϕp, v0)
param_sweep = [(Lx, Nx, d2(ϕa), d2(ϕp), v0) for (ϕa,ϕp) in zip(ϕa_sweep,ϕp_sweep)];
f,u,c = load_full(sweep...)

# sweep = (Lx, Nx, ϕa, ϕp, v0)
# param, ps = get_param_full(sweep...)
# err, erri, avmag, c = check_u_full(u,ps);
# if err>1e-6
#     global u,c
#     c = c*Lx
#     u[end] = u[end]*Lx
# end

for sweep in param_sweep
    local Lx, Nx, ϕa, ϕp, v0, param
    global f,u,c
    try
        Lx, Nx, ϕa, ϕp, v0 = sweep
        g,uu,c = load_full(sweep...)
        ϵ = 1e-4
        normf = sqrt(sum( (g[:,1] .- ϕa/2).^2 + (g[:,2] .- ϕa/2).^2 + (g[:,3] .- ϕp).^2)/Nx)
        if normf > ϵ
            f = g
            u = uu
        else
            f,u,c = solve_full(Lx,Nx,ϕa,ϕp,v0,u)
            normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
            if normf < ϵ
                break
            end
        end
    catch
        Lx, Nx, ϕa, ϕp, v0 = sweep
        f,u,c = solve_full(Lx,Nx,ϕa,ϕp,v0,u)
    end
end


# ϕa = 0.5
# ϕp = 0.3
# ϕa_end = [0.5, 0.0, 0.69][input]
# ϕp_end = [0.0, 0.3, 0.30][input]
# Lx     = [50.0][1] #[input]

# Lx, Nx, ϕa, ϕp, v0 = Lx, 3200, ϕa, ϕp, 7.5

# steps = Int64(max(abs(ϕa_end-ϕa)÷0.005,abs(ϕp_end-ϕp)÷0.005))

# Δa = (ϕa_end-ϕa)/steps
# Δp = (ϕp_end-ϕp)/steps

# if Δa == 0.0
#     ϕa_sweep = fill(ϕa,steps)
# else
#     ϕa_sweep = (ϕa+Δa):Δa:ϕa_end;
# end

# if Δp == 0.0
#     ϕp_sweep = fill(ϕp,steps)
# else
#     ϕp_sweep = (ϕp+Δp):Δp:ϕp_end;
# end


# sweep = (Lx, Nx, ϕa, ϕp, v0)
# param_sweep = [(Lx, Nx, d3(ϕa), d3(ϕp), v0) for (ϕa,ϕp) in zip(ϕa_sweep,ϕp_sweep)];

# f,u,c = load_full(sweep...)

# for sweep in param_sweep
#     local Lx, Nx, ϕa, ϕp, v0, param
#     global f,u,c
#         Lx, Nx, ϕa, ϕp, v0 = sweep
#         param, ps = get_param_full(sweep...)
#         f,u,c = solve_full(Lx,Nx,ϕa,ϕp,v0,u)
#         err, avmag, c = check_u_full(u,ps);
#         println("ϕa,ϕp=($(ϕa),$(ϕp)): c=$(c), avmag=$(avmag), err=$(err)")
# end

# ϕa = collect(0.29:0.01:0.53)[input]; # 25 long

# ϕp_sweep = 0.3:(0.01):0.4;
# # ϕp_sweep = 0.3:(-0.01):0.0;

# param_sweep = [(100.0, 1024, ϕa, ϕp, 7.5) for ϕp in ϕp_sweep];
# sweep = param_sweep[1]
# f,u,c = load_out_2(sweep...)

# for sweep in param_sweep
#     local Lx, Nx, ϕa, ϕp, v0
#     global f,u,c
#     try
#         Lx, Nx, ϕa, ϕp, v0 = sweep
#         param, p2 = get_outer_param_2(sweep...)
#         f,u,c = load_out_2(sweep...)
#         mag_av = Lx*sum(f[:,2]-f[:,1])/Nx
#         err = check_F(u,p2)
#         println("solved: Lx, Nx, ϕa, ϕp, v0 = $(sweep) | c=$(c), m_av =$(mag_av), err = $(err) ")
#     catch
#         Lx, Nx, ϕa, ϕp, v0 = sweep
#         param, p2 = get_outer_param_2(sweep...)
#         f,u,c = solve_out_2(Lx,Nx,ϕa,ϕp,v0,u)
#         mag_av = Lx*sum(f[:,2]-f[:,1])/Nx
#         err = check_F(u,p2)
#         println("solved: Lx, Nx, ϕa, ϕp, v0 = $(sweep) | c=$(c), m_av =$(mag_av), err = $(err) ")
#     end
# end

# load initial solution 
# Load initial wave
# Lx = 100.0
# Nx = 3200
# param = get_dense_param(Lx,Lx/Nx)
#     filename = steady_save_name(param)
#     data = load(filename)
#     @unpack f, c = data
# #
# print("loaded sol")
# for i in 1:2
#     global param, f, c
#     # Densify
#     f, param = double_sol(param,f)

#     # Set vars
#         @unpack DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = param
#         ps = (DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp)
#         u0 = get_u(f,c)
#     #

#     # Set problem
#         using NonlinearSolve, DifferentialEquations
#         prob = NonlinearProblem(ff,u0, ps; abstol = 1e-8, reltol =  1e-8);
#     #
#     print("solving problem")
#     # Solve 
#     sol  = solve(prob)

#     # Save
#     u = sol.u
#         f = get_f(u)
#         c = u[end]
#         filename    = steady_save_name(param)
#         data        = Dict("f" => f, "c" => c)
#         safesave(filename,data)
#     #
#     print("saved")
# end

# Lenghten 
# for Lx in 80.0:5.0:100.0
#     global param, f, c
#     stretch_param(param,Lx);

#     # Set vars
#         @unpack DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = param
#         ps = (DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp)
#         u0 = get_u(f,c)
#     #

#     # Set problem
#         using NonlinearSolve, DifferentialEquations
#         prob = NonlinearProblem(ff,u0, ps; abstol = 1e-8, reltol =  1e-8);
#     #

#     # Solve 
#     sol  = solve(prob)

#     # Save
#     u = sol.u
#         f = get_f(u)
#         c = u[end]
#         filename    = steady_save_name(param)
#         data        = Dict("f" => f, "c" => c)
#         safesave(filename,data)
#     #
# end


# using NonlinearSolve, DifferentialEquations
# I = 6
# for i in 1:I
#     global param, f, c
#         f, param = double_sol(param,f)

#         u0 = get_u(f,c)
#         @unpack DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp = param
#         ps = (DT, v0, DR, Δx, Nx, Lx, ϕa, ϕp)
#         #
#         if input ==1
#             using Symbolics
#             F0 = copy(u0)
#             jac_sparsity = Symbolics.jacobian_sparsity((F, u) -> ff(F, u, ps), F0, u0)
#             f_sparse = NonlinearFunction(ff; sparsity = jac_sparsity)
#             prob = NonlinearProblem(ff_sparse,u0, ps; abstol = 1e-8, reltol =  1e-8);
#         else
#             prob = NonlinearProblem(ff,u0, ps; abstol = 1e-8, reltol =  1e-8);
#         end
#         #
#         sol  = solve(prob)
#         u = sol.u
#         f = get_f(u)
#         c = u[end]
#         filename    = steady_save_name(param)
#         data        = Dict("f" => f, "c" => c)
#         safesave(filename,data)
#         print(param["Δx"])
# end

# i_max = 31
# i = (input %i_max) +1
# j = (input ÷ i_max)+1

# Lx = collect(10:1:80)[input] #71 long
# Δx = 0.01
# save_interval = 5.0

# param = get_dense_param(Lx, Δx; save_interval = save_interval)
# loaded, f, t = quiet_load_last_pde(param)
# save_interval = 10.0
# @pack! param = save_interval
# relax_sol(param,f,t; threshold = 1e-8)

# densify(Lx, Δx; save_interval = save_interval, threshold = 1e-6)

# run solution i j
# param = get_grid_param_wide(i,j)
# param["T"] = 800.0

# params = [ get_grid_param_wide(i,j) for i in 1:31 for j in 1:20];
# cut_off = 200.0
# loads = [quiet_load_last_pde(param) for param in params];
# no_loads = [param for (param,load) in zip(params,loads) if (load[3]<cut_off)&(param["ϕa"]+param["ϕp"] < 0.95)]

# param = no_loads[input]

# function get_stretch_param(Lx)
#     param = get_grid_param(21,11)
#     @unpack Nx = param
#     param["save_interval"] = 100.0
#     param["T"] = 2000.0
#     param["name"] = "soliton_stretch"
#     param["Lx"] = Float64(Lx)
#     param["Δx"] = Float64(Lx/Nx)
#     return param
# end
# Lx = Int64(input)
# param = get_stretch_param(Lx)
# load_and_run_pde(param)

# if param["ϕa"]+param["ϕp"] < 0.95
#     print("running solution $(param["ϕa"]) $(param["ϕp"])")
#     run_new_pde(param)
#     # load_and_run_pde(param)
# end


# loaded, f, t = load_last_pde(param)

# if !loaded
#     global t, f, param
#     param["name"] = "soliton_vertical_sweep"
#     loaded, f, t = load_last_pde(param)
#     param["name"] = "soliton_grid"
# end
# param["save_interval"] = 100.0
# if loaded
#     global t, f
#     dt = max(1000.0 - t,0)
#     print("running solution $(param["ϕa"]) $(param["ϕp"])")
#     t, f = run_current_pde(param,dt, f,t)
# else
#     print("invalid load $(param["ϕa"]) $(param["ϕp"]) ")
# end
