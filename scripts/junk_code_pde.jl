cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/pde_functions.jl")
include("/home/jm2386/Active_Lattice/src/plot_functions.jl")


#Example 
#=
#Parameters
name = "random_test_fixed_2"
ρa = 0.6
λ = 10.
t = 0.701
Dθ = 10.
Nx = 50
T = 0.03
param = pde_param(;name = name, 
        λ = λ , ρa = ρa, ρp = 0., T = T, Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = 20, save_interval = 0.01, max_steps = 1e8
)

perturb_pde_run(param)

density = initialize_density(param)
@unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S = param
@unpack fa, fp, t = density
perturb_pde!(param,density; pert = "n=1")
for i in 1:100 pde_step!(param,density) end
ρ = fp + 2*π*sum(fa;dims = 3)[:,:,1]/Nθ
p.(ρ)

@unpack T, save_interval, max_steps = param
T = 0.5

density = initialize_density(param)
perturb_pde!(param,density; pert = "rand", δ = 0.1);
perturb_pde!(param,density; pert = "n=1");
run_pde_until!(param,density,T; save_on = true, max_steps = max_steps, save_interval = save_interval)

midpoint_bond_diff(fp; Nx = Nx)
midpoint_bond_diff_θ(fa;  Nx = Nx,  Nθ = Nθ)
midpoint_bond_av(coeff_mag_s(fa,ρ); Nx =Nx )

@time for i in 1:10 pde_step!(param,density) end

pde_step!(param,density)
fig, axs = plt.subplots(1, 2, figsize=(12,5))
plot_pde(fig,axs,param,density)
display(fig)


t = 0.01

@unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";;;
        data = wload(filename)
        fa = data["fa"]
        fp = data["fp"]
        #t = data["t"]
        @pack! density = fa, fp, t
    fig, axs = plt.subplots(1, 2, figsize=(12,5))
    plot_pde(fig,axs,param,density)
display(fig)
t += 0.01

density = initialize_density(param)
perturb_pde!(param,density; pert = "rand", δ = 0.1);
#perturb_pde!(param,density; pert = "n=1");

pde_step!(param,density)
@unpack fa = density
m  = mag(fa; Nθ = Nθ, Nx =Nx);
maximum(abs.(m[:,1:24,1]+m[:,49:(-1):26,1]))
maximum(abs.(sum(fa[1:24,1,:]; dims = 2)-sum(fa[49:(-1):26,1,:]; dims = 2)))

ρ =  sum(fa[:,:,:]; dims = 3)[:,:,1]
m[:,14,16]

fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(m[2,:,1])
ax.plot((1/50):(1/50):(1.),ρ[:,18,1])
display(fig)


density
m2 = zeros(Nx,Nx,Nθ)

for x in 1:Nx, θ in 1:Nθ
    m2[:,x,y] = fa[:,y,x]

end


magx = []
magy = []
for t= 0.:0.01:0.4
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";;;
    data = wload(filename)
    fa = data["fa"]
    m  = mag(fa; Nθ = Nθ, Nx =Nx)
    av_m = sum(m; dims = (2,3))/(Nx*Nx)
    push!(magx, av_m[1])
    push!(magy, av_m[2])
end
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(0.:0.01:0.4,magx)
ax.plot(0.:0.01:0.4,magy)
display(fig)


maxmag = []
for t= 0.:0.01:0.4
    filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";;;
    data = wload(filename)
    fa = data["fa"]
    m  = mag(fa; Nθ = Nθ, Nx =Nx);
    maximum(abs.(m[:,1:24,:]+m[:,49:(-1):26,:]))
    push!(magx, av_m[1])
    push!(magy, av_m[2])
end
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(0.:0.01:0.4,magx)
ax.plot(0.:0.01:0.4,magy)
display(fig)










@unpack fa, fp, t = density
dist_from_unif(param,fa,fp)

#expand variables
@unpack name, D, λ, ρa, ρp, δt, Nx, Nθ, S = param
@unpack fa, fp, t = density
#Run and save
run_pde_until!(param,density,0.002; save_on = true, max_steps = 20, save_interval = 0.001)
#for pmap
param = pde_param(λ = 0, T = 0.003, name = "test8")
perturb_pde_run(param; max_steps = 60, save_interval = 0.0005)
#Loading
t_saves, fa_saves, fp_saves = load_pdes(param,0.003; save_interval = 0.0005)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
#plotting
fig, ax = PyPlot.subplots(figsize =(10, 10))
y = self_diff.(x)
ax.plot(y)
display(fig)
density = initialize_density(param)
n = 1#length(t_saves)
t, fa, fp = t_saves[n], fa_saves[n], fp_saves[n]
@pack! density = fa, fp, t
#
pde_step!(param,density)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_pde_mass(fig,ax,param,density)
display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_pde_mag(fig,ax,param,density)
display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(t_saves,dist_saves)
display(fig)
##
pde_step!(param,density)
@unpack fa, fp, t = density
dist_from_unif(param,fa,fp)


for x₂ in 1:Nx
        local dsx, ρx, magx, y₁,y₂
        ## 2 direction

        y₂ = (x₂ +Nx)%Nx +1

        println(y₂)
end
=#


#########################
ρa = 0.05
T = 0.03
max_runs = 2
param = pde_param(; name = name, λ = λ , ρa = ρa, ρp = 0., T = T, Dθ = 100., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = max_runs, λ_step = 5., λmax = 204., λs=λs)
run_stab_search(param)


# /store/DAMTP/jm2386/Active_Lattice/data/pde_raw/
t_saves, fa_saves, fp_saves = load_pdes(param,0.2; save_interval = 0.001)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
fig, ax = PyPlot.subplots(figsize =(10, 10))
n = 4
t, fa, fp = t_saves[n], fa_saves[n], fp_saves[n]

display(fig)
#plot_pde_mag(fig,ax,param,density)
#display(fig)
density = Dict{String,Any}()
name = "high_density_stability_v5"
ρa = 0.6
λ = 50.
t = 0.701
Dθ = 100.
param = pde_param(;name = name, 
        λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8
)
@unpack name, Nx, Nθ, λ, ρa, ρp, δt = param
        filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_raw/$(name)/Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ)/time=$(round(t; digits = 2))_Nx=$(Nx)_Nθ=$(Nθ)_active=$(ρa)_passive=$(ρp)_lamb=$(λ)_dt=$(δt)_Dθ=$(Dθ).jld2";;;
        data = wload(filename)
        fa = data["fa"]
        fp = data["fp"]
        t = data["t"]
        @pack! density = fa, fp, t
        fig, ax = PyPlot.subplots(figsize =(10, 10))
        plot_pde_mass(fig,ax,param,density)
        #plot_pde_mag(fig,ax,param,density)
        display(fig)
#
fig, ax = PyPlot.subplots(figsize =(10, 10))
param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.01)
dist_saves = time_dist_from_unif(param, fa_saves, fp_saves)
fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(t_saves,dist_saves)
display(fig)
#
PyPlot.close()
fig = figure(figsize=(30,10))
i = 1 
name = "high_density_stability_v4"
for λ in [50., 40., 30., 20., 10.]
for ρa in collect(0.05:0.05:0.85)
        ax = fig[:add_subplot](5,17,i)
        param = pde_param(;name = name, λ = λ , ρa = ρa, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e8)
        t_saves, fa_saves, fp_saves = load_pdes(param,1.0; save_interval = 0.001)
        plot_error(fig,ax, param, t_saves, fa_saves, fp_saves; t_max = 1.0, y_max = 0.2)
        if λ == 10.
                ax.set_xlabel("t")
        end
        if ρa == 0.05
                ax.set_ylabel("e")
        end
        i += 1
end
end
fig.tight_layout()
display(fig)
PyPlot.savefig("/store/DAMTP/jm2386/Active_Lattice/data/time_stability_test.pdf",dpi = 300, format = "pdf")

### high density Dθ = 10
stabdata = refresh_stab_data(; ρs = 0.05:0.05:0.95,   Dθ = 10., Nx = 50, Nθ = 20, λs = 5.:5.:100., name = "high_density_stability_v4")
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.9:0.01:0.99 ,xs = collect(0.01:0.01:0.99), xtic = 0.9:0.1:1, ytic = 0:10:100, axlim = [0.9, 1., 15., 100.])
display(fig)
###
λs = 20.:10.:200.
Dθ = 100.
name = "high_density_stability_v5"
stabdata = refresh_stab_data(; ρs = 0.05:0.05:1.0,   Dθ = Dθ, Nx = 50, Nθ = 20, λs = λs, name = name)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.05:0.05:0.95 ,xs = collect(0.01:0.01:0.99), xtic = 0:0.2:1, ytic = λs , axlim = [0., 1., 30., 300.], Dθ = Dθ)
display(fig)
###
stabdata = refresh_stab_data(; ρs = [0.05],   Dθ = Dθ, Nx = 50, Nθ = 20, λs = λs, name = name)
@unpack name, Nx, Nθ = param
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ).jld2"
#wsave(filename,stabdata)
stabdata = wload(filename)
###
ρa = 0.6
param = pde_param(; name = name , ρa = ρa, ρp = 0., T = 1.0, Dθ = Dθ, δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7, max_runs = 10, λ_step = 20., λmax = 404., λs=λs)
next_param(stabdata, param; λmax = 100., λ_step = 20.)
#
ρ = 0.95
param = pde_param(; name = name, λ = λ , ρa = ρ, ρp = 0., T = 1.0, Dθ = 10., δt = 1e-5, Nx = 50, Nθ = 20, save_interval = 0.01, max_steps = 1e7)
#
     
using CurveFit

ρs = 0.05:0.05:0.95
xs = collect(0.4:0.01:0.99)
binodal_y = Array{Float64,1}([])
binodal_x = Array{Float64,1}([])
for ρ in ρs
            @unpack stable, unstable, unsure = stabdata["ρ = $(ρ)"]
        try
                local y
                y = (maximum(stable)+minimum(unstable))/2
                push!(binodal_y,y)
                push!(binodal_x,ρ)
        catch
        end
end 

function lineared(binodal_x,binodal_y,xs)
        new_biny = Array{Float64,1}([])
        new_binx = Array{Float64,1}([])
        for i in 1:(length(binodal_x)-1)
                local m 
                m = (binodal_y[i+1]-binodal_y[i])/(binodal_x[i+1]-binodal_x[i])
                for x in xs
                        if (x≥ binodal_x[i])&(x<binodal_x[i+1])
                                local y
                                y = binodal_y[i] + m*(x-binodal_x[i])
                                push!(new_binx,x)
                                push!(new_biny,y)
                        end 
                end
        end
        return new_binx, new_biny
end

binodal_x,binodal_y = lineared(binodal_x,binodal_y,xs)

fig, ax = PyPlot.subplots(figsize =(10, 10))
ax.plot(binodal_x,binodal_y)
display(fig)

fit = curve_fit(RationalPoly, collect(binodal_x), collect(binodal_y), 8,1)
aprox_binodal_y = fit.(binodal_x)
aprox_binodal_x = binodal_x
ax.plot(aprox_binodal_x,aprox_binodal_y)
display(fig)



x = 0.0:0.02:2.0
y0 = @. 1 + x + x*x + randn()/10
fit = curve_fit(Polynomial, collect(x), y0, 2)
y0b = fit.(x) 




unfinished, next, param = next_param(stabdata, param; λmax = 105., λ_step = 5.)


λsym(0.9; Dθ=100., Dx=1. 