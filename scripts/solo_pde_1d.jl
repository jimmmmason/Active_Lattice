cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions_1d.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
name = "stability_1d_actpass_2"
Pe = 30.
ρa = 0.55
Dθ = 1.
ρp = 0.2
    δ = 1e-2
    δt = 1e-5
    T  = 1.0
    save_interval = 0.01
    Nx = 50
    Nθ = 20

param = pde_param_1d(; 
    name = name, ρa = ρa, ρp = ρp, Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = save_interval, δ= δ
)
fig, ax = PyPlot.subplots(figsize =(10, 10))
t_saves, fa_saves, fp_saves = load_pdes_1d(param,T; save_interval = 0.0001)
    dist_saves = time_dist_from_unif_1d(param, fa_saves, fp_saves)
    ax.plot(t_saves,dist_saves)
    ax.set_xlabel("t")
    ax.set_ylabel("‖ρ-ϕ‖₂")
    ax.set_title("Dθ = $(Dθ) ρ = $(ρa) Pe = $(Pe)")
display(fig)

lin_stab_line(0.4;Dx = Dx, Pe =Pe, Dθ = 1)

#non-linear instability? 
name = "stability_1d_5"
Pe = 3.
ρa = 0.6
Dθ = 100.
    δ = 1e-2
    δt = 1e-5
    T  = 0.4
    save_interval = 0.01
    Nx = 50
    Nθ = 20
#




####



n = length(fa_saves)
fa = fa_saves[n]
fp = fp_saves[n]

param = pde_param(; 
    name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = 0.0001, δ= δ
)

name = "stability_1d_actpass"
Pe = 10.
ρa = 0.45
ρp = 0.1
Dθ = 100.
    δ = 1e-2
    δt = 1e-5
    T  = 1.0
    save_interval = 0.01
    Nx = 50
    Nθ = 20

param = pde_param_1d(; 
    name = name, ρa = ρa, ρp = ρp, Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = save_interval, δ= δ
)


density = initialize_density_1d(param)
@unpack T, save_interval, max_steps, pert, δ = param
perturb_pde_1d!(param,density; pert = pert, δ = δ)
@unpack fa, fp, t = density
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(fa)
display(fig)

density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
perturb_pde!(param,density; pert = "n=1",δ = δ);

@time for i in 1:100 pde_stepper_1d!(param,density) end 

@unpack fa = density
maximum(fa)


pde_stepper!(param,density)

@unpack fa, fp, t = density
maximum(abs.(fa1[:,:] - fa[:,1,:]))

i = 1
fa = fa_saves[i]
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    #@unpack fa1, fp1, t = density1
    #ρ = sum(fa1; dims =2)[:,1].*(2*π/Nθ)
    #ax.plot(ρ)
    #fa1 = density1["fa"];
    ax.plot(fa)
    #ρ = sum(fa1; dims =2)[:,1].*(2*π/Nθ)
    #ax.plot(ρ)
    i+=1
display(fig)


dt = pde_stepper_1d!(param1,density1)
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    @unpack fa1, fp1, t = density1
    m = mag_1d(fa1;Nx = Nx, Nθ = Nθ)
    ax.plot(m)
display(fig)


@unpack δt, Nx, Nθ, Dθ, λ,pert = param1
perturb_pde_1d!(param1,density1; pert = pert, δ = δ)

@unpack δt, Nx, Nθ, Dθ, λ,pert = param1
@unpack fa1, fp1, t = density1
ρ =  sum(fa1; dims =2)[:,1].*(2*π/Nθ)


logmfa1 = map(x -> (x>0 ? log(x) : logtol), fa1);
logmfp1 = map(x -> (x>0 ? log(x) : logtol), fp1);
p_rho  = p.(ρ);

Ua  = -midpoint_bond_diff_θ_1d(logmfa1 .+ p_rho; Nx=Nx, Nθ=Nθ).+ λ*midpoint_bond_av_1d(coeff_mag_s_1d(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ) .+ λ*eθ 
Up  = -midpoint_bond_diff_1d(  logmfp  + p_rho; Nx=Nx       ) + λ*midpoint_bond_av_1d(coeff_mag_s_1d(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx )
Uθ  = -midpoint_Θ_diff_1d(fa; Nx=Nx, Nθ = Nθ)




Ua1,   Up,   Uθ   = U_velocities_1d(fa1,fp1,ρ; Nx=Nx, Nθ=Nθ, λ=λ)
    moba, mobp, mobθ = mob_1d(fa1,fp1,ρ)
    Fa,   Fp,   Fθ   = F_fluxes_1d(Ua1, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)

    a = maximum(abs.(Ua1));
    b = maximum(abs.(Up));
    c = maximum(abs.(Uθ));

    tempu = 1/(6*max(a*Nx, b*Nx, c*Nθ*Dθ/(2*π)))
    dt = min(δt, tempu)

    dfa1 = site_div_θ_1d(Fa; Nx=Nx, Nθ=Nθ)
    dfthe1 = site_θ_diff_1d(Fθ; Nx=Nx, Nθ=Nθ)

    fa1 -= dt*(dfa1 + Dθ*dfthe1)
    #fp1 -= dt*site_div_1d(fp1; Nx=Nx)
    ρ = fp1 + sum(fa1; dims =2)[:,1].*(2*π/Nθ)

    @pack! density1 = fa1, fp1, t

dist_from_unif_1d(param1, fa1, fp1)





t_saves, fa1_saves, fp1_saves = load_pdes_1d(param1,T; save_interval = 0.0001)
    dist_saves = time_dist_from_unif_1d(param1, fa1_saves, fp1_saves)
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    ax.plot(t_saves,dist_saves)
    ax.set_xlabel("t")
    ax.set_ylabel("‖ρ-ϕ‖₂")
    ax.set_title("Dθ = $(Dθ) ρ = $(ρa) Pe = $(Pe)")
display(fig)









density1 = initialize_density1_1d(param1)
#perturb_pde!(param1,density1; pert = "rand",δ = δ);
perturb_pde!(param1,density1; pert = "n=1",δ = δ);

pde_stepper_1d!(param1, density1)


@unpack fa1, fp1, t = density1
dist_from_unif(param1, fa1, fp1)
pde_step!(param1,density1)
@unpack fa1, fp1, t = density1
dist_from_unif(param1, fa1, fp1)



name = "linear_stability"
Pe = 1000.
    ρa = 0.7
    Dθ = 100.
    δ = 1e-4
    δt = 1e-7
    T  = 0.001
    λ = Pe*sqrt(Dθ)
Nx = 100
Nθ = 20
param1 = pde_param_1d(; 
    name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = 0.0001, δ= δ
)
density1 = initialize_density1(param1)
#perturb_pde!(param1,density1; pert = "rand",δ = δ);
perturb_pde!(param1,density1; pert = "n=1",δ = δ);
@unpack fa1, fp1, t = density1
logtol = log(1e-10);
eθ = reshape([cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)]',2,1,1,Nθ);

    ρ = fp1 + sum(fa1; dims =3)[:,:,1].*(2*π/Nθ);

    logmfa1 = map(x -> (x>0 ? log(x) : logtol), fa1);
    logmfp1 = map(x -> (x>0 ? log(x) : logtol), fp1);
    p_rho  = p.(ρ);#functon p is labelled W in the pdf

    Ua1 =  zeros(2,Nx,Nx,Nθ) .+ λ*eθ # -midpoint_bond_diff_θ(logmfa1 .+ p_rho; Nx=Nx, Nθ=Nθ) .+ λ*midpoint_bond_av(coeff_mag_s(fa1,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ) #
    Up = zeros(2,Nx,Nx) #-midpoint_bond_diff(  logmfp1  + p_rho; Nx=Nx       ) + λ*midpoint_bond_av(coeff_mag_s(fa1,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx )
    Uθ = zeros(Nx,Nx,Nθ)# -midpoint_Θ_diff(fa1; Nx=Nx, Nθ = Nθ)

    moba, mobp, mobθ = mob(fa1,fp1,ρ)
    fa1,   fp1,   Fθ   = F_apθ(Ua1, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)

    fa1 -= δt*(site_div_θ(fa1; Nx=Nx, Nθ=Nθ))# + Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ))

dist_from_unif(param1, fa1, fp1)-δ
sqrt(2*π*sum( (site_div_θ(fa1; Nx=Nx, Nθ=Nθ) ).^2)/(Nx*Nx*Nθ))


density1 = initialize_density1(param1)
#perturb_pde!(param1,density1; pert = "rand",δ = δ);
perturb_pde!(param1,density1; pert = "n=1",δ = δ);
@unpack fa1, fp1, t = density1
fa12,fp12,dt = time_stepper(fa1, fp1, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)
dist_from_unif(param1, fa12, fp12)

density1 = initialize_density1(param1)
#perturb_pde!(param1,density1; pert = "rand",δ = δ);
perturb_pde!(param1,density1; pert = "n=1",δ = δ);
pde_stepper!(param1,density1)
@unpack fa1, fp1, t = density1
dist_from_unif(param1, fa1, fp1)

maximum(midpoint_bond_av(coeff_mag_s(fa1,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ))


maximum(abs.(Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ)))



