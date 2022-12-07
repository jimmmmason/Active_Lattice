cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")
###
name = "linear_stability"
    Pe = 50.
    ρa = 0.9
    Dθ = 100.
    δ = 1e-4
    δt = 1e-7
    T  = 0.001
Nx = 50
Nθ = 20

param = pde_param(; 
    name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = 0.0001, δ= δ
)

density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
perturb_pde!(param,density; pert = "n=1",δ = δ);
pde_stepper!(param,density)
pde_stepper!(param,density)
pde_stepper!(param,density)
pde_stepper!(param,density)
pde_stepper!(param,density)
pde_stepper!(param,density)


perturb_pde_run(param)

t_saves, fa_saves, fp_saves = load_pdes(param,T; save_interval = 0.0001)
    dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
    fig, ax = PyPlot.subplots(figsize =(10, 10))
    ax.plot(t_saves,dist_saves)
    ax.set_xlabel("t")
    ax.set_ylabel("‖ρ-ϕ‖₂")
    ax.set_title("Dθ = $(Dθ) ρ = $(ρa) Pe = $(Pe)")
display(fig)





density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
δ = 0.05
δt = 1e-6
@pack! param = δt
perturb_pde!(param,density; pert = "rand",δ = δ);

@unpack fa, fp, t = density
dist_from_unif(param, fa, fp)
pde_stepper!(param,density)
@unpack fa, fp, t = density
dist_from_unif(param, fa, fp)
ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ)
maximum(ρ)

for i in 1:811 pde_stepper!(param,density) end

@unpack t,fa,fp = density
@unpack δt, Nx, Nθ, Dθ, λ = param
fa, fp , dt = time_stepper(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)

Ua, Up, Uθ = U_velocities(fa,fp,ρ; Nx=Nx, Nθ=Nθ, λ=λ);
moba, mobp, mobθ =  mob(fa,fp,ρ);

a = maximum(abs.(Ua))
b = maximum(abs.(Up))
c = maximum(abs.(Uθ));
1/(6*max(a*Nx, b*Nx, c*Nθ*Dθ/(2*π)))

fp
ds
minimum(moba)
minimum(mobp)
Fa,Fp,Fθ =  F_fluxes(Ua, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ);
maximum(abs.(Fa))*2*pi
maximum(abs.(Fp))


1/(2*pi)
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
param = pde_param(; 
    name = name, ρa = ρa, ρp = 0., Pe = Pe, T = T, Dθ = Dθ, δt = δt, Nx = Nx, Nθ = Nθ, save_interval = 0.0001, δ= δ
)
density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
perturb_pde!(param,density; pert = "n=1",δ = δ);
@unpack fa, fp, t = density
logtol = log(1e-10);
eθ = reshape([cos.((1:Nθ)*2π/Nθ) sin.((1:Nθ)*2π/Nθ)]',2,1,1,Nθ);

    ρ = fp + sum(fa; dims =3)[:,:,1].*(2*π/Nθ);

    logmfa = map(x -> (x>0 ? log(x) : logtol), fa);
    logmfp = map(x -> (x>0 ? log(x) : logtol), fp);
    p_rho  = p.(ρ);#functon p is labelled W in the pdf

    Ua =  zeros(2,Nx,Nx,Nθ) .+ λ*eθ # -midpoint_bond_diff_θ(logmfa .+ p_rho; Nx=Nx, Nθ=Nθ) .+ λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ) #
    Up = zeros(2,Nx,Nx) #-midpoint_bond_diff(  logmfp  + p_rho; Nx=Nx       ) + λ*midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx )
    Uθ = zeros(Nx,Nx,Nθ)# -midpoint_Θ_diff(fa; Nx=Nx, Nθ = Nθ)

    moba, mobp, mobθ = mob(fa,fp,ρ)
    Fa,   Fp,   Fθ   = F_apθ(Ua, Up, Uθ, moba, mobp, mobθ; Nx=Nx, Nθ=Nθ)

    fa -= δt*(site_div_θ(Fa; Nx=Nx, Nθ=Nθ))# + Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ))

dist_from_unif(param, fa, fp)-δ
sqrt(2*π*sum( (site_div_θ(Fa; Nx=Nx, Nθ=Nθ) ).^2)/(Nx*Nx*Nθ))


density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
perturb_pde!(param,density; pert = "n=1",δ = δ);
@unpack fa, fp, t = density
fa2,fp2,dt = time_stepper(fa, fp, δt; Nx=Nx, Nθ=Nθ, λ=λ, Dθ=Dθ)
dist_from_unif(param, fa2, fp2)

density = initialize_density(param)
#perturb_pde!(param,density; pert = "rand",δ = δ);
perturb_pde!(param,density; pert = "n=1",δ = δ);
pde_stepper!(param,density)
@unpack fa, fp, t = density
dist_from_unif(param, fa, fp)

maximum(midpoint_bond_av(coeff_mag_s(fa,ρ; Nθ=Nθ, Nx=Nx ); Nx =Nx ))


maximum(abs.(Dθ*site_θ_diff(Fθ; Nx=Nx, Nθ=Nθ)))



