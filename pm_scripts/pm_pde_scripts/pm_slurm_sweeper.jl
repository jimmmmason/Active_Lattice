#Set up
cd("/home/jm2386/Active_Lattice/");
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_sims.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");
include("/home/jm2386/Active_Lattice/src/Hetrocline.jl");

# Access the command-line argument 'i'
input = parse(Int64, ARGS[1]); # 1445, 200, 346 missing 

function tp_param(i)
    DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 25.0, 0.5, 0.1,0.1,1e-5, 0.1);
    ϕa,ϕp = [(x,y) for x in 0.1:0.02:1.0, y in 0.0:0.01:0.4 if x+y<1][i]
    T, save_interval, param_name, pert = (2000.0, 10.0, "tp_experiment", "single")
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
    return param
end
function cp_param_coord(ϕa,ϕp)
    DT, v0, DR, N, Δx, Lx, Ly, _, _, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 25.0, 0.5, 0.1,0.1,1e-5, 0.1);
    T, save_interval, param_name, pert = (2000.0, 10.0, "cp_experiment", "double")
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
    return param
end
function tp_param_coord(ϕa,ϕp)
    DT, v0, DR, N, Δx, Lx, Ly, _, _, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 25.0, 0.5, 0.1,0.1,1e-5, 0.1);
    T, save_interval, param_name, pert = (2000.0, 10.0, "tp_experiment", "single")
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
    return param
end
function cp_param(i)
    DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ = (1.0, 7.5, 1.0, 100, 0.05, 25.0, 0.5, 0.1,0.1,1e-5, 0.1);
    ϕa,ϕp = [(x,y) for x in 0.1:0.02:1.0, y in 0.0:0.01:0.4 if x+y<1][i]
    T, save_interval, param_name, pert = (2000.0, 10.0, "cp_experiment", "double")
    param = new_param(DT, v0, DR, N, Δx, Lx, Ly, ϕa, ϕp, δt, δ; T = T, name = param_name, save_interval = save_interval, save_on = true, pert= pert)
    return param
end

# i = [87, 88, 131, 132, 175, 219][input]
# param = cp_param(i)

# i = [506, 805, 835, 840, 874, 875, 904, 905, 906, 908, 934, 935, 974, 975, 1189, 1190, 1192, 1275, 1276, 1277, 1278, 1303, 1304, 1305, 1329, 1330, 1331, 1332, 1357][input]
i = [545, 583, 584, 585, 622, 835, 856, 857, 909, 962, 963, 964, 965, 996, 997, 998, 999, 1000, 1186, 1187, 1188, 1191, 1246, 1271, 1272, 1273, 1274, 1300, 1301, 1302, 1358][input]
param = tp_param(i)

# ϕa, ϕp = [(x,y) for x in 0.2:0.02:0.4, y in 0.01:0.01:0.07][input]
# param = cp_param_coord(ϕa, ϕp)

# param = [cp_param_coord(0.54, 0.05)] #cp_param_coord(0.58, 0.04),cp_param_coord(0.6, 0.03),cp_param_coord(0.64, 0.03),cp_param_coord(0.66, 0.02)][input]

# run_new_pde(param)
load_and_run_pde(param)

# ind = collect(16:8:1008)[input] #125

# function load_adjacent_6(Lx, Nx, ϕa, ϕp, v0, ind, γ) # f,u,c,ps,ind = / "fail"
#     ind_range = [(i,j) for j in [0,1,-1,2,-2], i in [0,1,-1]]
#     input_params = [(Lx,Nx,ϕa,ϕp,v0, Int64(ind+8*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
#     ϵ = 1e-4
#     for input_param in input_params
#         try
#             Lx,Nx,ϕa,ϕp,v0,ind,γ = input_param
#             param, ps = get_outer_param_6(Lx,Nx,ϕa,ϕp,v0,ind,γ)
            
#             f,u,c = load_out_6(input_param...)
#             ϕp = sum(f)/Nx-sum(f[:,1:2])/Nx
#             ϕa = sum(f[:,1:2])/Nx

#             normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#             err, erri, avmag, cep = check_u_6(u,ps);
#             rhodiff = u[ind+1] - u[ind]

#             if (normf>ϵ)&&(err<ϵ)&&(rhodiff>ϵ)
#                 print("loaded: $((Lx,Nx,d2(ϕa),d2(ϕp),v0,ind,γ))")
#                 return f,u,c,ps,ind
#             end
#         catch
#         end
#     end
#     return "fail"
# end

# if ind < 120
#     global Lx,Nx,ϕa,ϕp,v0,ind,γ
#     Lx,Nx,ϕa,ϕp,v0,ind,γ = 100.0, 1024, 0.3, 0.3, 7.5, ind, 1.8
#     ind_range = [(i,0) for i in collect(-1:-1:-14)]
# elseif ind < 328
#     global Lx,Nx,ϕa,ϕp,v0,ind,γ
#     Lx,Nx,ϕa,ϕp,v0,ind,γ = 100.0, 1024, 0.3, 0.3, 7.5, ind, 1.75
#     ind_range = [(i,0) for i in collect(-1:-1:-14)]
# elseif ind < 904
#     global Lx,Nx,ϕa,ϕp,v0,ind,γ
#     Lx,Nx,ϕa,ϕp,v0,ind,γ = 100.0, 1024, 0.3, 0.3, 7.5, ind, 1.65
# else
#     global Lx,Nx,ϕa,ϕp,v0,ind,γ
#     Lx,Nx,ϕa,ϕp,v0,ind,γ = 100.0, 1024, 0.3, 0.3, 7.5, ind, 1.65
# end

# range_max = (2.2-γ)÷0.01
# ind_range = [(0,i) for i in collect(1:1:(range_max))]
# input_params = [(Lx,Nx,ϕa,ϕp,v0, Int64(ind),d2(γ+0.01*j) ) for (i,j) in ind_range]

# # range_max = (γ-1)÷0.01
# # ind_range = [(0,i) for i in collect(-1:-1:(-range_max))]
# # input_params = [(Lx,Nx,ϕa,ϕp,v0, Int64(ind),d2(γ+0.01*j) ) for (i,j) in ind_range]

# ϵ = 1e-5;

# # load adjecnt then solve; stop if err/small 
# for input_param in input_params
#     global Lx,Nx,ϕa,ϕp,v0,ind,γ,f,u,c
#     output = load_adjacent_6(input_param...);
#     if output == "fail"
#         print("load_fail: $(input_param)")
#         break 
#     else
#         f,u,c,p2,ind2 = output;
#         Lx,Nx,ϕa,ϕp,v0,ind,γ = input_param
#         param, ps = get_outer_param_6(Lx,Nx,ϕa,ϕp,v0,ind,γ)
        
#         if ind2 > ind
#             u,ps = decrease_inds(u,p2,ind2-ind)
#         elseif ind2 < ind
#             u,ps = increase_inds(u,p2,ind-ind2)
#         end
#         f,u,c = solve_out_6(Lx,Nx,ϕa,ϕp,v0,ind,γ,u; tol = 1e-8, maxiters = 10)
#         ϕp = d2(sum(f)/Nx-sum(f[:,1:2])/Nx)
#         ϕa = d2(sum(f[:,1:2])/Nx)

#         param, ps = get_outer_param_6(Lx,Nx,ϕa,ϕp,v0,ind,γ)
#         err, erri, avmag, cep = check_u_6(u,ps)
#         rhodiff = u[ind+1] - u[ind]
#         normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)

#         if (normf < ϵ)
#             print("norm end: $(input_param)")
#             break
#         elseif (err > ϵ)
#             print("err end: $(input_param)")
#             break
#         elseif (rhodiff < ϵ)
#             print("rho end: $(input_param)")
#             break
#         elseif (c < ϵ^2)
#             print("speed end: $(input_param)")
#             break
#         end
#     end
# end





# # outer_5 starter ϕ sweep
# function load_adjacent_5(Lx, Nx, ϕa, ϕp, v0, ϕ, γ) # f,u,c = / "fail"
#     ind_range = [(i,j) for j in [0,1,-1], i in [0,1,-1,2,-2]]
#     input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
#     ϵ = 1e-4
#     for input_param in input_params
#         try
#             Lx,Nx,ϕa,ϕp,v0,ϕ,γ = input_param
#             param, ps = get_outer_param_5(input_param...)
#             f,u,c = load_out_5(input_param...)
#             ϕp = sum(f)/Nx-sum(f[:,1:2])/Nx
#             ϕa = sum(f[:,1:2])/Nx
#             normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#             err, erri, avmag, cep = check_u_5(u,ps)
#             if (normf > ϵ)&&(err<ϵ)
#                 print("loaded: $((Lx,Nx,d2(ϕa),d2(ϕp),v0,ϕ,γ))")
#                 return f,u,c
#             end
#         catch
#         end
#     end
#     return "fail"
# end
# # start point
# Lx,Nx,ϕa,ϕp,v0,ϕ,_ = 100.0, 1024, 0.3, 0.3, 7.5, 0.56, 2.2
# γ = collect(1.26:0.01:2.2)[input] # 95
# # sweep range
# ind_range = [(i,0) for i in collect(-1:-1:-18)]
# input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
# ϵ = 1e-4
# # ind_range = [(i,0) for i in collect(1:41) ] 
# # input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
# # ϵ = 1e-4
# # load adjecnt then solve; stop if err/small 
# for input_param in input_params
#     global Lx,Nx,ϕa,ϕp,v0,ϕ,γ,f,u,c
#     output = load_adjacent_5(input_param...);
#     if output == "fail"
#         print("load_fail: $(input_param)")
#         break 
#     else
#         f,u,c = output;
#         Lx,Nx,ϕa,ϕp,v0,ϕ,γ = input_param
        
#         f,u,c = solve_out_5(Lx,Nx,ϕa,ϕp,v0,ϕ,γ,u; tol = 1e-8, maxiters = 10)
#         ϕp = d2(sum(f)/Nx-sum(f[:,1:2])/Nx)
#         ϕa = d2(sum(f[:,1:2])/Nx)

#         param, ps = get_outer_param_5(Lx,Nx,ϕa,ϕp,v0,ϕ,γ)
#         err, erri, avmag, cep = check_u_5(u,ps)
#         normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)

#         if (normf < ϵ)
#             print("norm end: $(input_param)")
#             break
#         elseif (err > ϵ)
#             print("err end: $(input_param)")
#             break
#         end
#     end
# end




# # outer_3 starter v sweep
# function load_adjacent_3(Lx, Nx, ϕa, ϕp, v0) # f,u,c = / "fail"
#     ind_range = [(i,j) for j in [0,1,-1], i in [0,1,-1,2,-2]]
#     input_params = [(Lx, Nx, d2(ϕa + 0.01*i), d2(ϕp +0.01*j), v0) for (i,j) in ind_range]
#     ϵ = 1e-4
#     for input_param in input_params
#         try
#             f,u,c = load_out_3(input_param...)
#             normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#             param, p2 = get_outer_param_3(input_param...)
#             err, avmag, c = check_u_3(u,p2);
#             if (normf > ϵ)&&(err<ϵ)
#                 print("loaded: $(input_param)")
#                 return f,u,c
#             end
#         catch
#         end
#     end
#     return "fail"
# end
# # start point
# ind_range = [(i,j) for j in [0], i in collect(-10:1:10)] #21 
# ϕa, ϕp = [ (d2(0.45 + 0.01*i), d2(0.35 +0.01*j)) for (i,j) in ind_range][input]
# Lx, Nx, _, _, v0 = 100.0, 1024, 0.45, 0.35, 7.5
# # sweep range
# ind_range = [(i,i) for i in collect(-1:-1:-20)]
# input_params = [(Lx, Nx, d2(ϕa + (-0.36/25)*i), d2(ϕp +0.01*j), v0) for (i,j) in ind_range]
# ϵ = 1e-4
# # load adjecnt then solve; stop if err/small 
# for input_param in input_params
#     output = load_adjacent_3(input_param...);
#     if output == "fail"
#         print("load_fail: $(input_param)")
#         break 
#     else
#         f,u,c = output;
#         Lx, Nx, ϕa, ϕp, v0 = input_param
#         f,u,c = solve_out_3(Lx,Nx,ϕa,ϕp,v0,u; tol = 1e-8)
#         param, ps = get_outer_param_3(input_param...)
#         err, avmag, c = check_u_3(u,ps);
#         normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#         if (normf < ϵ)|(err > ϵ)
#             print("end: $(input_param)")
#             break
#         end
#     end
# end



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

# # h push 
# # ϕa = 0.27
# ϕa = 0.27
# ϕp = d2.(collect(0.4:-0.01:0.0))[input] #41

# Lx, Nx, ϕa, ϕp, v0 = 100.0, 3200, ϕa, ϕp, 7.5
# ϕa_sweep = d2.(collect(ϕa:(-0.01):0.0))
# # ϕa_sweep = d2.(collect(ϕa:(0.01):(1-ϕp)))


# f,u,c = load_full(Lx, Nx, 0.5, 0.3, v0)
# no_load = true
# j = length(ϕa_sweep)
# ϵ = 1e-4
# while (no_load)&(j>0)
#     global f,u,c,ϕp,Lx,Nx,ϕa,ϕp,v0,no_load,j,dj
#     try
#         ϕa = ϕa_sweep[j]
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

# ϕa_sweep = d2.(collect(ϕa:-0.01:0.0))
# # ϕa_sweep = d2.(collect(ϕa:0.01:(1-ϕp)))
# ϕp_sweep = fill(ϕp,length(ϕa_sweep))

# sweep = (Lx, Nx, ϕa, ϕp, v0)
# param_sweep = [(Lx, Nx, d2(ϕa), d2(ϕp), v0) for (ϕa,ϕp) in zip(ϕa_sweep,ϕp_sweep)];
# f,u,c = load_full(sweep...)

# # sweep = (Lx, Nx, ϕa, ϕp, v0)
# # param, ps = get_param_full(sweep...)
# # err, erri, avmag, c = check_u_full(u,ps);
# # if err>1e-6
# #     global u,c
# #     c = c*Lx
# #     u[end] = u[end]*Lx
# # end

# for sweep in param_sweep
#     local Lx, Nx, ϕa, ϕp, v0, param
#     global f,u,c
#     try
#         Lx, Nx, ϕa, ϕp, v0 = sweep
#         g,uu,c = load_full(sweep...)
#         ϵ = 1e-4
#         normf = sqrt(sum( (g[:,1] .- ϕa/2).^2 + (g[:,2] .- ϕa/2).^2 + (g[:,3] .- ϕp).^2)/Nx)
#         if normf > ϵ
#             f = g
#             u = uu
#         else
#             f,u,c = solve_full(Lx,Nx,ϕa,ϕp,v0,u)
#             normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
#             if normf < ϵ
#                 break
#             end
#         end
#     catch
#         Lx, Nx, ϕa, ϕp, v0 = sweep
#         f,u,c = solve_full(Lx,Nx,ϕa,ϕp,v0,u)
#     end
# end


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
