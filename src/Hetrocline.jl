cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
println("Loading ...")
##
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#
#functions 
using QuadGK


rtol = 1e-14
tol = 1e-14

R_prime(ρ) = -1/self_diff(ρ)/(1-ρ)^3

function R(ρ; rtol=rtol, tol = tol)
    α::Float64= π/2 -1;
    c::Float64 =  -1/(α*(2*α-1)/(2*α+1) - α +1);
    f(x) = R_prime(x) - c/(1-x)^4;
    return quadgk(f, 0.0, ρ, rtol=rtol, atol = tol)[1] + c/3/(1-ρ)^3 - c/3
end


function R_old(ρ; rtol=rtol)
    return quadgk(x -> R_prime(x),0.0, ρ, rtol=rtol)[1]
end

# x = collect(0.9:0.001:0.999)
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# ax.plot(x,R.(x)-R_old.(x), color = "blue", label = "new")
# display(fig)

κ(ρ,Pe) = self_diff(ρ) / Pe / (1 - ρ)
Λ(ρ,Pe) = -2 * self_diff(ρ) / Pe / (1 - ρ)^2
g0(ρ;Pe = 10, γ = 1, atol = 1e-10 ) = Pe*( 1-γ*(1- ρ) )*self_diff(ρ) - 2*log(1 - ρ)/Pe
dg0(ρ;Pe = 10, γ = 1 ) = Pe*(γ )*self_diff(ρ) +Pe*( 1-γ*(1- ρ) )*self_diff_prime(ρ) + 2/(1 - ρ)/Pe

dΦ_dρ(ρ,Pe,γ) = g0(ρ;Pe = Pe,γ = γ)*R_prime(ρ)

function ΦoR(x; Pe = 10, γ = 3.4,rtol=rtol, tol = tol) # Φ as a funciton of ρ ie composed with R
    α::Float64= π/2 -1;
    c::Float64 =  2/Pe/(α*(2*α-1)/(2*α+1) - α +1);
    f(x) = dΦ_dρ(x,Pe,γ) - c*log(1 - x)/(1-x)^4;
    return quadgk(f, 0.95, x, rtol=rtol, atol = tol)[1]+c*(1+3*log(1-x))/9/(1-x)^3-c*(1+3*log(1-0.95))/9/(1-0.95)^3
end

# function ΦoR_old(x; Pe = 10,γ = 3.4, rtol=rtol) # Φ as a funciton of ρ ie composed with R
#     f(x) = dΦ_dρ(x,Pe,γ);
#     return quadgk(f, 0.95, x, rtol=rtol)[1]
# end

# # prev h alt
# function h0_alt(x; Pe = 10,γ = 3.4, rtol=rtol, atol = atol) # Φ as a funciton of ρ ie composed with R
#     f(x) = dg0(x;Pe = Pe, γ = γ )*R(x; tol = atol, rtol = rtol);
#     return quadgk(f, 0.95, x, rtol=rtol, atol = atol)[1]
# end


# x = 0.8
# Pe = 10
# γ = 3.5
# f(x) = dg0(x;Pe = Pe, γ = γ )*R(x; tol = atol, rtol = rtol);
# x = collect(0.1:0.01:0.99)
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# ax.plot(x,dg0.(x; Pe = Pe,γ = γ), color = "blue", label = "new")
# #ax.plot(x,R.(x)/10, color = "red", label = "new")
# ax.plot(x,f.(x)/1000, color = "green", label = "new")
# # ax.plot(x,ΦoR.(x;Pe = 10,γ = 3.4)-ΦoR_old.(x;Pe = 10,γ = 3.4), color = "blue", label = "new")
# ax.axis([0.5,1,-1.0,1.0])
# display(fig)

# Pe = 7

function h0(ρ;Pe = 10,γ = 1, atol = 1e-12)
    return g0(ρ;Pe = Pe,γ = γ)*R(ρ; tol = atol)-ΦoR(ρ;Pe = Pe,γ = γ, tol = atol)
end

# function h0_old(ρ;Pe = 10,γ = 1)
#     return g0(ρ;Pe = Pe,γ = γ)*R_old(ρ)-ΦoR_old(ρ)
# end



#### h,g plot

# x = collect(0.0:0.01:0.99)
# x = append!(x,collect(0.991:0.001:0.999))
# x = append!(x,collect(0.9991:0.0001:0.9999))
# x = append!(x,collect(0.99991:0.00001:0.99999))

# Pe = 3.9
# G = g0.(x;Pe = Pe)
# H = h0.(x;Pe = Pe)

# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# ax.plot(x,10*G, color = "red", label = "G")
# ax.plot(x,H.-H[1], color = "blue", label = "H")
# ax.axis([0.0,1,-10,20])
# display(fig)


# x = collect(0.999991:0.0000001:0.9999999)
# Pe = 5
# #G = g0.(x;Pe = Pe)
# H = h0.(x;Pe = Pe)
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
#     #ax.plot(x,5000*G, color = "red", label = "G")
#     ax.plot(-log.(-x.+1),log.(abs.(H)), color = "blue", label = "H")
#     ax.plot(-log.(-x.+1),-2.9*log.(-x.+1).+log.(abs.(H))[1].+2.9*log.(-x.+1)[1], color = "black", label = "-3", linestyle = "--")
#     ax.plot(-log.(-x.+1),-3*log.(-x.+1).+log.(abs.(H))[1].+3*log.(-x.+1)[1], color = "black", label = "-3", linestyle = "--")
#     ax.plot(-log.(-x.+1),-3.1*log.(-x.+1).+log.(abs.(H))[1].+3.1*log.(-x.+1)[1], color = "black", label = "-3.1",linestyle = "--")
#     ax.legend()
#     #ax.axis([0.99,1,-10000,10000])
# display(fig)


# PyPlot.close("all")
####

using Roots

function find_G_density(ϕ1; Pe = 5, γ = 1, limits = (0.8,0.9999) , atol = atol)
    try
        f(x) = g0(x;Pe = Pe,γ = γ, atol = atol)-g0(ϕ1; Pe = Pe,γ = γ, atol = atol)
        ϕ2 = find_zero(f, limits)
        return ϕ2
    catch
        return "fail"
    end
end

function find_G_density_no_error(ϕ1; Pe = 5, γ = 1, limits = (0.8,0.9999), fail = 0., atol = 1e-12 )
    try
        f(x) = g0(x;Pe = Pe,γ = γ, atol = atol)-g0(ϕ1; Pe = Pe,γ = γ, atol = atol)
        ϕ2 = find_zero(f, limits)
        return ϕ2
    catch
        return fail
    end
end

function find_H_density_no_error(ϕ1; Pe = 5, γ = 1, limits = (0.8,0.9999), fail = 0. , atol = 1e-12)
    try
        f(x) = h0(x;Pe = Pe,γ = γ, atol = atol)-h0(ϕ1;Pe = Pe,γ = γ, atol = atol)
        ϕ2 = find_zero(f, limits)
        return ϕ2
    catch
        return fail
    end
end

function find_zero_density_difference(; Pe = 5, γ = 1, lower_limits = (0.35,0.4), upper_limits = (0.8,0.9999), atol = 1e-12)
    f(x) = find_H_density_no_error(x; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[2], atol = atol) - find_G_density_no_error(x; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[1], atol = atol)
    ϕ1 = find_zero(f, lower_limits)
    ϕ2 = find_G_density(ϕ1; Pe = Pe, γ = γ, limits = upper_limits)
    return ϕ1, ϕ2
end

function find_zero_density_difference_2(; Pe = 5, γ = 1, lower_limits = (0.35,0.4), upper_limits = (0.8,0.9999))
    f(x) = find_H_density_no_error(x; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[2]) - find_G_density_no_error(x; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[1])
    ϕ2 = find_zero(f, upper_limits)
    ϕ1 = find_G_density(ϕ1; Pe = Pe, γ = γ, limits = upper_limits)
    return ϕ1, ϕ2
end


function make_x(;npoints=100, start_pt = 0., end_pt = 1.)
    length = (end_pt-start_pt)
    Δx = length/npoints
    x1 = collect(start_pt:Δx:end_pt)
    x2 = -length*exp.(-collect(1:npoints)).+end_pt
    x = append!(x1,x2)
    return sort(x)
end

function g_tunring_points(;Pe = 5, γ = 3.52, initial_Δ = 1e-8, atol = atol)
    x = 0.:(initial_Δ):(1-initial_Δ)
    gmin, gmax = (0., 1.)
    DG = dg0.(x; Pe = Pe, γ= γ)
    DG_min = minimum(DG)
    gmin = 1.0
    gmax = 0.0
    if DG_min> 0
        return println("no turn in g0")
    else
        i = argmin(DG)
        x0 = x[i]
        f(x) = dg0(x; Pe = Pe, γ= γ)
        gmax = find_zero(f, (0.,x0))
        gmin = find_zero(f, (x0,1.))
    end
    return gmin,gmax
end

function h_maximum(;Pe = 5., γ = 1., initial_Δ = initial_Δ, atol = atol)
    x = 0.:(initial_Δ):(1-initial_Δ)
    return x[argmax(h0.(x; Pe = Pe, γ= γ, atol = atol))]
end

function h_minimum(hmax;Pe = 5., γ = 1., initial_Δ = initial_Δ, atol = atol)
    X = 0.:(initial_Δ):(1-initial_Δ)
    x = [x for x in X if x < hmax]
    return x[argmin(h0.(x; Pe = Pe, γ= γ, atol = atol))]
end

function initial_intervals(;Pe = 5, γ = 1, rho_max = (1-10e-16), initial_Δ = 1e-5, atol = atol)
    local lower_limits, upper_limits
    gmin, gmax = g_tunring_points(;Pe = Pe, γ = γ, initial_Δ = initial_Δ, atol = atol)
    lower_1 = find_G_density_no_error(gmin; Pe = Pe, γ = γ, limits = (0.0,gmax), fail = 0., atol = atol )
    upper_2 = find_G_density_no_error(gmax; Pe = Pe, γ = γ, limits = (gmin,rho_max), fail = 1.0, atol = atol )

    h_max = h_maximum(;Pe = Pe, γ = γ, initial_Δ = initial_Δ, atol = atol)

    lower_limits = (lower_1, gmax)
    upper_limits = (max(gmin,h_max), upper_2)
    return lower_limits, upper_limits
end

function restrict_intervals_h(lower_limits, upper_limits; Pe = 5, γ = 1, atol = atol)
    upper_1 = find_H_density_no_error(lower_limits[1]; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[1], atol = atol )
    upper_2 = find_H_density_no_error(lower_limits[2]; Pe = Pe, γ = γ, limits = upper_limits, fail = upper_limits[2], atol = atol )
    upper_limits = (min(upper_1,upper_2), max(upper_1,upper_2))
    return lower_limits, upper_limits
end

function restrict_intervals_g(lower_limits, upper_limits; Pe = 5, γ = 1, atol = atol)
    lower_1 = find_G_density_no_error(upper_limits[1]; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[1], atol = atol )
    lower_2 = find_G_density_no_error(upper_limits[2]; Pe = Pe, γ = γ, limits = lower_limits, fail = lower_limits[2], atol = atol )
    lower_limits = (min(lower_1,lower_2), max(lower_1,lower_2))
    return lower_limits, upper_limits
end

function interval_size(lower_limits, upper_limits)
    return abs(lower_limits[2]-lower_limits[1])+abs(upper_limits[2]-upper_limits[1])
end

function colapse_sol_interval(;Pe = 5, γ = 2.9, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 20, tol = 1e-8, atol = 1e-12)
    lower_limits, upper_limits = initial_intervals(;Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ, atol = atol)
    precision = interval_size(lower_limits, upper_limits)
    i = 0
    while (i<max_iter)&(precision>tol)
        lower_limits, upper_limits = restrict_intervals_h(lower_limits, upper_limits;Pe = Pe, γ = γ, atol = atol)
        lower_limits, upper_limits = restrict_intervals_g(lower_limits, upper_limits;Pe = Pe, γ = γ, atol = atol)
        precision = interval_size(lower_limits, upper_limits)
        i+=1
    end
    if precision ≤ tol
        find_sol = true
    else
        find_sol = false
    end

    return find_sol, lower_limits, upper_limits
end

function colapse_sol_approx(;Pe = 5, γ = 2.9, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 20, tol = 1e-8)
    hmax = h_maximum(;Pe = Pe, γ = γ, initial_Δ = initial_Δ)
    hmin = h_minimum(hmax;Pe = Pe, γ = γ, initial_Δ = initial_Δ)

    upper = find_H_density_no_error(hmin; Pe = Pe, γ = γ, limits = (hmax,rho_max), fail = hmax)
    gmin, lower = g_tunring_points(;Pe = Pe, γ = γ, initial_Δ = initial_Δ)
    try 
        lower = find_G_density_no_error(hmin; Pe = Pe, γ = γ, limits = (0.1,lower), fail = lower)
    catch
    end
    return lower, upper
end



# Pe = 10
# γ = 3.5
# x = make_x(npoints=1000, start_pt = 0., end_pt = 0.999)
# G0 = g0.(x; Pe = Pe, γ= γ)
# mG = minimum(G0)
# H0 = h_alt.(x; Pe = Pe, γ= γ, atol = 1e-14).-h_alt.(0.1; Pe = Pe, γ= γ, atol = 1e-14)
# H0 = h0.(x; Pe = Pe, γ= γ, atol = 1e-14).-h0.(0.1; Pe = Pe, γ= γ, atol = 1e-14)
# mH = maximum(H0)

# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# ax.plot(x,-G0/mG/0.1 .-0.4, color = "red", label = "g0")
# ax.plot(x,-H0/mH/100, color = "blue", label = "h0")
# ax.axis([0,1,-2,2])
# display(fig)




# Pe = 10
# γ = 3.4
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# dG0 = dg0.(x; Pe = Pe, γ= γ)
# ax.plot(x,dG0, color = "red", label = "g0")
# ax.axis([0.6,1,-0.1,0.1])
# display(fig)



Pe = 10
γ = 3.54
find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 10, tol = 1e-4, atol = 1e-12)


find_zero_density_difference(; Pe = Pe, γ = γ, lower_limits = lower_limits, upper_limits = upper_limits, atol = 1e-12)



#rtol = 1e-12
# function ΦoR(x; Pe = 10, γ = 1,rtol=rtol) # Φ as a funciton of ρ ie composed with R
#     α::Float64= π/2 -1;
#     c::Float64 =  2/Pe/(α*(2*α-1)/(2*α+1) - α +1);
#     f(x) = dΦ_dρ(x,Pe,γ) - c*log(1 - x)/(1-x)^4;
#     return quadgk(f, 0.95, x, rtol=rtol)[1]+c*(1+3*log(1-x))/9/(1-x)^3-c*(1+3*log(1-0.95))/9/(1-0.95)^3
# end

# function limit_finder(;Pe = 5, γ = 1,  rho_max = (1-10e-16), x = append!(collect(0.01:0.01:0.95),collect(0.951:0.00001:0.9999)))
#     ### g limits
#     # find high stationary point
#         x_high = [x for x in x if x>0.5]
#         index_min = argmin( g0.(x_high;Pe = Pe, γ = γ) )
#         if index_min ==1
#             println("error Pe too small -> no solution")
#         end
#     gx_min = x_high[index_min]

#     gϕ_min = find_G_density(gx_min; Pe = Pe, γ = γ, limits = (0.0,0.5) )
#     if gϕ_min == "fail"
#         gϕ_min = 0.
#     end



#     # find low stationary point
#         x_low = [x for x in x if x<0.5]
#         index_max = argmax( g0.(x_low;Pe = Pe, γ = γ) )
#     gx_max = x_low[index_max]

#     if gx_max>0.
#         gϕ_max = find_G_density(gx_max; Pe = Pe, γ = γ, limits = (gx_min,rho_max) )
#         if gϕ_max-gx_min == "fail"
#             println("need higher rho cur off")
#             gϕ_max = 1.
#         end
#     else
#         gϕ_max = 1.
#         gx_max = gx_min
#     end
#     #

#     ### h limits
#     #
#     # find high stationary point
#         x_high = [x for x in x if x>0.5]
#         index_max = argmax( h0.(x_high;Pe = Pe, γ = γ) )
#     hx_max = x_high[index_max]
#     if hx_max == maximum(x)
#         println("need higher rho cur off")
#         hx_max = gx_min
#     end
#     # find low stationary point
#         x_low = [x for x in x if x<hx_max]
#         index_min = argmin( h0.(x_low;Pe = Pe, γ = γ) )
#     hx_min = x_low[index_min]

#     hϕ_max = find_H_density(hx_min; Pe = Pe, γ = γ, limits = (hx_max,rho_max) )
#     if hϕ_max == "fail"
#         hϕ_max = 1.0
#     end
#     #
#     lower_limits = (gϕ_min, gx_max )
#     #upper_limits = (max(gx_min,hx_max), min(hϕ_max,gϕ_max) )
#     upper_limits = (max(gx_min,hx_max), rho_max )
#     #
#     # extra bonus
#     #upper_high_x = upper_limits[2]
#     #lower_high_x = find_G_density(upper_high_x-1e-16; Pe = Pe, γ = γ, limits = lower_limits )
#     #
#     #lower_limits = (gϕ_min, min(lower_limits[2],lower_high_x) )
#     return lower_limits, upper_limits
# end

# function solution_check(;Pe = 5, γ = 1,  rho_max = (1-10e-16), x = append!(collect(0.01:0.01:0.95),collect(0.951:0.00001:0.9999)))
#     solution = true
#     ### h limits
#     #
#     # find high stationary point
#         x_high = [x for x in x if x>0.5]
#         index_max = argmax( h0.(x_high;Pe = Pe, γ = γ) )
#     hx_max = x_high[index_max]
    
#     # find low stationary point
#         x_low = [x for x in x if x<hx_max]
#         index_min = argmin( h0.(x_low;Pe = Pe, γ = γ) )
#     hx_min = x_low[index_min]

#     hϕ_max = find_H_density(hx_min; Pe = Pe, γ = γ, limits = (hx_max,rho_max) )
#     if hϕ_max == "fail"
#         hϕ_max = 1.0
#         solution = false
#     end
#     #

#     ### g limits
#     # find high stationary point
#     x_high = [x for x in x if x>0.5]
#     index_min = argmin( g0.(x_high;Pe = Pe, γ = γ) )
#     gx_min = x_high[index_min]

#     gϕ_min = find_G_density(gx_min; Pe = Pe, γ = γ, limits = (0.0,0.5) )
#     if gϕ_min == "fail"
#         gϕ_min = 0.
#         solution = false
#     end
#         # find low stationary point
#         x_low = [x for x in x if x<0.5]
#         index_max = argmax( g0.(x_low;Pe = Pe, γ = γ) )
#         gx_max = x_low[index_max]

#         if gx_max>0.
#             gϕ_max = find_G_density(gx_max; Pe = Pe, γ = γ, limits = (gx_min,rho_max) )
#             if gϕ_max == "fail"
#                 println("need higher rho cut off")
#                 gϕ_max = 1.
#             end
#         else
#             gϕ_max = 1.
#             gx_max = gx_min
#             solution = false
#         end

#     if g0(gx_max;Pe = Pe, γ = γ) < g0(hx_max;Pe = Pe, γ = γ)
#         solution = false
#     end


#     if g0(gϕ_min;Pe = Pe, γ = γ) > g0(hx_max;Pe = Pe, γ = γ)
#         solution = false
#     end

#     lower_limits = (gϕ_min, gx_max )
#     upper_limits = (max(gx_min,hx_max), rho_max )

#     return solution, lower_limits, upper_limits
# end

#### Calc binodal 

# Pe  = 5 
# γ = 1
# lower_limits, upper_limits = limit_finder(Pe = Pe, γ = γ)
# ϕ1, ϕ2 = find_zero_density_difference(;Pe = Pe, γ = γ, lower_limits = lower_limits, upper_limits = upper_limits)

###check solution
# h0(ϕ1; Pe = Pe)
# h0(ϕ2; Pe = Pe)

# g0(ϕ1; Pe = Pe)
# g0(ϕ2; Pe = Pe)


### Stream plot

# ϵ = 1e-10
# Δx = 0.0001
# Δy = 0.0001

# x = collect((ϕ1-0.1+ϵ):Δx:(ϕ1+0.6))
# y = collect((-0.2):Δy:0.2)

# xx = [x̃ for x̃ ∈ x, ỹ ∈ y]
#     yy = [ỹ for x̃ ∈ x, ỹ ∈ y]
#     Nx = length(x)
#     Ny = length(y)
#     v = zeros(Nx,Ny,2)
#     for i in 1:Nx, j in 1:Ny
#         ρ = x[i]
#         m = y[j]
#         v[i,j,1] = Pe*(1-ρ)*m
#         v[i,j,2] = -Pe*m^2 + Pe*(ρ*self_diff(ρ) -ϕ1*self_diff(ϕ1) )/self_diff(ρ) -(2/Pe)*log( (1-ρ)/(1-ϕ1) )/self_diff(ρ)
#     end
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
#     max_value = 1
#     absmag  = min.(sqrt.(v[:,:,1].^2+v[:,:,2].^2), max_value)
#     colmap = PyPlot.plt.cm.viridis
#     streamlines = ax.streamplot(xx', yy', v[:,:,1]', v[:,:,2]', color = absmag', cmap = colmap, density = 4)#2.5
#     norm1 = matplotlib.colors.Normalize(vmin=0., vmax= max_value);
#     fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm1, cmap = colmap), ax = ax, fraction = 0.0455)
# pertubation = [ϕ1+ϵ 0.0]
# ax.streamplot(xx', yy', v[:,:,1]', v[:,:,2]', color="red", start_points= pertubation)
# # fig options
# rc("text", usetex=true)
#     #ax.yaxis.set_ticks(ytic)
#     axlim = [0.3, 0.95, -0.2, 0.2]
#     ax.xaxis.set_tick_params(labelsize=15)
#     ax.xaxis.tick_bottom()
#     ax.yaxis.set_tick_params(labelsize=15)
#     ax.axis(axlim)
#     #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
#     ax.set_xlabel(L"\rho",fontsize=20)
#     ax.set_ylabel(L"m", fontsize=20)
#     #ax.legend(loc = "upper left", fontsize=20)
#     #ax.set_aspect(0.25*Δρ/ΔPe)
#     #title = latexstring("\$ \\ell = $(round(1/sqrt(Dθ); digits = 2)), \\chi = $(χ) \$")
#     #ax.set_title(title,fontsize=20)
#     fig.tight_layout()
# display(fig)
# #save fig
# name = "Streamplot_Pe=$(Pe)"
#     pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)";
#     mkpath(pathname)
#     filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)/Streamplot.pdf";
#     PyPlot.savefig(filename,dpi = 100, format = "pdf") #bbox_extra_artists=( ldg,)
# ##



###

using DifferentialEquations

function f(du,u,parameters,t)
    Pe = parameters[1]
    γ = parameters[2]
    ϕ1 = parameters[3]
    du[1] = Pe*(1-u[1])*u[2]
    du[2] = -Pe*u[2]^2 + Pe*( (1-γ*(1-u[1]))*self_diff(u[1]) -(1-γ*(1-ϕ1))*self_diff(ϕ1) )/self_diff(u[1]) -(2/Pe)*log( (1-u[1])/(1-ϕ1) )/self_diff(u[1])
    return du
end

function f_jac(J,u,parameters,t)
    Pe = parameters[1]
    γ = parameters[2]
    ϕ1 = parameters[3]
    J[1,1] = -Pe*u[2]
    J[1,2] =  Pe*(1-u[1])
    J[2,1] =  self_diff_prime(u[1])*            (2/Pe)*log( (1-u[1])/(1-ϕ1) )/self_diff(u[1])^2              + (2/Pe)/(1-u[1])/self_diff(u[1])
    J[2,1] += self_diff_prime(u[1])*Pe*((1-γ*(1-u[1]))*self_diff(u[1]) -(1-γ*(1-ϕ1))*self_diff(ϕ1) )/self_diff(u[1])^2  + Pe*(γ*self_diff(u[1]) +(1-γ*(1-u[1]))*self_diff_prime(u[1]) )/self_diff(u[1])
    J[2,2] = -2*Pe*u[2]
    return J
end

# parameters = (Pe, γ = γ, ϕ1, ϕ2)
# ff = ODEFunction(f;jac=f_jac)
# ϵ = 1e-10
# initial_position = [ϕ1+ϵ, 0.0]
# time_interval = (0.0, 20.0)

# ff = ODEFunction(f;jac=f_jac)
# prob = ODEProblem(ff,initial_position,time_interval, parameters)

# sol = DifferentialEquations.solve(prob,abstol = 1e-12, reltol = 1e-12)


# PyPlot.close("all")
# fig, axs = plt.subplots(2, 1, figsize=(8,10))
# xlabels = [L"\rho", L"m"]

# t_mid_arg = argmax(sol[2,:])
# t_middle = sol.t[t_mid_arg]
# t_max = maximum(sol.t[:].-t_middle)
# t_min = minimum(sol.t[:].-t_middle)
# t_lim = round(min(t_max, - t_min))

# for i in 1:2
#     ax = axs[i]
#     ax.plot(sol.t[:].-t_middle, sol[i,:], color = "black")
#     axlim = [ -t_lim, t_lim, minimum(sol[i,:])-0.01, maximum(sol[i,:])+0.01]
#     ax.xaxis.set_tick_params(labelsize=15)
#     ax.xaxis.tick_bottom()
#     ax.yaxis.set_tick_params(labelsize=15)
#     ax.axis(axlim)
#     #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
#     ax.set_ylabel(xlabels[i],fontsize=20)
#     ax.set_xlabel(L"x",fontsize=20)
# end
# fig_title = latexstring("\$ \\mathrm{Pe} = $(Pe) \$")
# axs[1].set_title(fig_title,fontsize=20)
# display(fig)
#save fig
# name = "Sol_Pe=$(Pe)"
#     pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)";
#     mkpath(pathname)
#     filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)/Solution.pdf";
#     PyPlot.savefig(filename,dpi = 100, format = "pdf") #bbox_extra_artists=( ldg,)
# ##

# initial_position = [ϕ2-ϵ, 0.0]
# time_interval = (0.0, 20.0)

# ff = ODEFunction(f;jac=f_jac)
# prob = ODEProblem(ff,initial_position,time_interval, parameters)

# sol = DifferentialEquations.solve(prob,abstol = 1e-12, reltol = 1e-12)

# PyPlot.close("all")
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# ax.plot(sol.t[:], sol[1,:])
# ax.plot(sol.t[:], sol[2,:])
# display(fig)

###


Pe = 10.0
γ = 4.0

fig, ax = plt.subplots(1, 1, figsize=(10,10))
G0 = g0.(x; Pe = Pe, γ= γ).-g0.(0.1; Pe = Pe, γ= γ)
mG = minimum(G0)
H0 = h0.(x; Pe = Pe, γ= γ).-h0.(0.1; Pe = Pe, γ= γ)
mH = maximum(H0)
ax.plot(x,-G0/mG/10, color = "red", label = "g0")
ax.plot(x,H0/mH, color = "blue", label = "h0")
ax.axis([0,1,-2,2])
display(fig)
x[argmax(H0)]
x[argmin(G0)]

fig, ax = plt.subplots(1, 1, figsize=(10,10))
ax.plot(x,f.(x), color = "blue", label = "f")
display(fig)

function colapse_sol_interval_high(;Pe = 5, γ = 2.9, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 20, tol = 1e-8)
    lower_limits, upper_limits = initial_intervals_high(;Pe = Pe, γ = γ, rho_max = rho_max, initial_Δ = initial_Δ)
    precision = interval_size(lower_limits, upper_limits)
    i = 0
    while (i<max_iter)&(precision>tol)
        lower_limits, upper_limits = restrict_intervals_h(lower_limits, upper_limits;Pe = Pe, γ = γ)
        lower_limits, upper_limits = restrict_intervals_g(lower_limits, upper_limits;Pe = Pe, γ = γ)
        precision = interval_size(lower_limits, upper_limits)
        i+=1
    end
    if precision ≤ tol
        find_sol = true
    else
        find_sol = false
    end

    return find_sol, lower_limits, upper_limits
end

function initial_intervals_high(;Pe = 5, γ = 1, rho_max = (1-10e-16), initial_Δ = 1e-5)
    local lower_limits, upper_limits
    gmin, gmax = g_tunring_points(;Pe = Pe, γ = γ, initial_Δ = initial_Δ)
    lower_1 = find_G_density_no_error(gmin; Pe = Pe, γ = γ, limits = (gmax, gmin), fail = 0. )
    upper_2 = find_G_density_no_error(gmax; Pe = Pe, γ = γ, limits = (gmin,rho_max), fail = 1.0 )

    h_max = h_maximum(;Pe = Pe, γ = γ, initial_Δ = initial_Δ)

    lower_limits = (gmax,lower_1 )
    upper_limits = (max(gmin,h_max), upper_2)
    return lower_limits, upper_limits
end


Pe = 10.0
γ = 3.5222333
find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 10, tol = 1e-2)
        

Pe = 10.0
γ = 3.522
colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 10, tol = 1e-2)
        
gamma_converter(γ, 0.91)

γ = 3.5222333

Γ = collect(1:0.1:3.5)
Γ = append!(Γ,collect(3.51:0.01:3.522))
Γ = append!(Γ,collect(3.522:0.001:3.5222333))

γ = 3.5222333
for γ in Γ
    find_sol = false
    try 
        find_sol, lower_limits, upper_limits = colapse_sol_interval(;Pe = Pe, γ = γ, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 10, tol = 1e-2, atol = 1e-12)
        if find_sol
            push!(ϕ1s,lower_limits[1])
            push!(ϕ2s,upper_limits[1])
            push!(γs, γ)
        else
            # lower, upper = colapse_sol_approx(;Pe = Pe, γ = γ, rho_max = (1-10e-16), initial_Δ = 1e-4, max_iter = 40, tol = 1e-8)
            # push!(ϕ1s,lower)
            # push!(ϕ2s,upper)
            # push!(γs, γ)
            # println("approx solution Pe=$(Pe), γ=$(γ)")
        end
    catch
            println("no solution Pe=$(Pe), γ=$(γ)")
            push!(errors,γ)
    end
end

Pes = []
γs = []
ϕ1s = []
ϕ2s = []
approx = []
errors = []


average_ϕs = (ϕ1s+ ϕ2s)./2
# subset = 10:10:200
# average_ϕs = average_ϕs[subset]
# γ_subset = γs

function chi_converter(γ, ϕ)
    ϕa = 1 - γ*(1 - ϕ)
    return ϕa/ϕ
end

function chis_converter(γs, ϕs)
    map(γs,ϕs) do γ, ϕ
        chi_converter(γ, ϕ)
    end
end

χs = chis_converter(γs, average_ϕs)


PyPlot.close("all")

using JLD2
data = Dict{String,Any}()
@pack! data = Pe, γs, ϕ1s, ϕ2s
filename = "/store/DAMTP/jm2386/Active_Lattice/data/binodal/Pe=$(Pe).jld2"
wsave(filename,data)

data = wload(filename)
@unpack Pe, γs, ϕ1s, ϕ2s = data

Pe = 2.5

γs[25:25:250]
fig, ax = plt.subplots(1, 1, figsize=(10,10))
ax.plot(ϕ1s,γs, color = "red", label = "Binodal")
ax.plot(ϕ2s,γs, color = "red", label = "Bindoal")
display(fig)

function gamma_converter(γ, ϕ)
    ϕa = 1 - γ*(1 - ϕ)
    ϕp = ϕ-ϕa
    return ϕa, ϕp
end

function gammas_converter_a(γs, ϕs)
    map(γs,ϕs) do γ, ϕ
        gamma_converter(γ, ϕ)[1]
    end
end

function gammas_converter_p(γs, ϕs)
    map(γs,ϕs) do γ, ϕ
        gamma_converter(γ, ϕ)[2]
    end
end

fig, ax = plt.subplots(1, 1, figsize=(10,10))
binod = ax.plot(gammas_converter_a(γs, ϕ1s), gammas_converter_p(γs, ϕ1s), color = "red", label = "Binodal")
ax.plot(gammas_converter_a(γs, ϕ2s), gammas_converter_p(γs, ϕ2s), color = "red", label = "_Bindoal")
ax.plot(0.:0.1:1., 1.:-0.1:0., color = "black", label = "_Full")
#
x = -collect(ϕ1s[2]:0.001:ϕ2s[2]).+1
γ = γs[2]
ax.plot(-γ*x.+1, (γ-1)*x, color = "gray", label = "Cloud line")
for i in 4:2:24
    x = -collect(ϕ1s[i]:0.001:ϕ2s[i]).+1
    γ = γs[i]
    ax.plot(-γ*x.+1, (γ-1)*x, color = "gray", label = "_cloud line")
end
#fig, ax = plt.subplots(1, 1, figsize=(10,10))
ϕas = collect(0.0001:0.0001:0.99)
    ϕa_sols, low_ϕps, high_ϕps = return_complex_boundary_outer(ϕas; Pe = Pe)
    n = length([x for x in ϕa_sols if x < 0.42])
    N = length(ϕa_sols)
    ax.plot(ϕa_sols[1:n], low_ϕps[1:n], color = "green", label = "Complex")
    ax.plot(ϕa_sols[1:n], high_ϕps[1:n], color = "green", label = "_Complex")
    ax.plot(ϕa_sols[(n+1):1:N], low_ϕps[(n+1):1:N], color = "green", label = "_Complex")
    ax.plot(ϕa_sols[(n+1):1:N], high_ϕps[(n+1):1:N], color = "green", label = "_Complex")


    ϕa_sols, low_ϕps1, high_ϕps1, low_ϕps2, high_ϕps2 = return_complex_boundary_inner(ϕas; Pe = Pe)
    ax.plot(ϕa_sols, low_ϕps1, color = "green", label = "_Complex")
    ax.plot(ϕa_sols, high_ϕps1, color = "green", label = "_Complex")
    ax.plot(ϕa_sols, low_ϕps2, color = "green", label = "_Complex")
    ax.plot(ϕa_sols, high_ϕps2, color = "green", label = "_Complex")
    
    ϕa_low, ϕa_high, low_ϕps2, high_ϕps2 = return_stable_boundary(ϕas; Pe = Pe)
    ax.plot(ϕa_low, low_ϕps2, color = "blue", label = "Spinodal")
    ax.plot(ϕa_high, high_ϕps2, color = "blue", label = "_Spindoal")

    ϕa_low, ϕa_mid, ϕa_high, low_ϕps, mid_ϕps, high_ϕps = return_stable_boundary_extra(ϕas; Pe = Pe)
    ax.plot(ϕa_low, low_ϕps, color = "blue", label = "_Spinodal")
    ax.plot(ϕa_high, high_ϕps, color = "blue", label = "_Spindoal")
    ax.plot(ϕa_mid, mid_ϕps, color = "blue", label = "_Spindoal")

    # ax.axis([0,1, 0,1])
    # display(fig)
# # fig options
rc("text", usetex=true)
    #ax.yaxis.set_ticks(ytic)
    axlim = [0, 1, 0, 1]
    #axlim = [0.4, 0.45, 0, 0.1]
    ax.xaxis.set_tick_params(labelsize=15)
    ax.xaxis.tick_bottom()
    ax.yaxis.set_tick_params(labelsize=15)
    ax.axis(axlim)
    #ax.set_title(L"\Re{ \lambda_n^\mathrm{max}} = 0",fontsize=20)
    ax.set_xlabel(L"\phi^a",fontsize=20)
    ax.set_ylabel(L"\phi^p", fontsize=20)
    ax.legend(loc = "upper right", fontsize=20)
    #ax.set_aspect(0.25*Δρ/ΔPe)
    title = latexstring("\$ \\mathrm{Pe} = $(Pe)\$")
    ax.set_title(title,fontsize=20)
    fig.tight_layout()
display(fig)
#save fig
name = "TheoreticalPhasePlot_Pe=$(Pe)"
    pathname = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)";
    mkpath(pathname)
    filename = "/store/DAMTP/jm2386/Active_Lattice/plots/active_passive_1d/$(name)/Streamplot.pdf";
    PyPlot.savefig(filename,dpi = 100, format = "pdf") #bbox_extra_artists=( ldg,)
##

low_ϕps2[1543:10:1583]

ϕa_low[1563]



function α(ϕa, ϕp)
    ϕ  = ϕa + ϕp
    ϕ0 = 1- ϕ
    ds = self_diff(ϕ)
    dsp = self_diff_prime(ϕ)
    return ϕa*ϕ0*(ds+dsp*ϕ)+ds^2*ϕp
end

function β(ϕa, ϕp)
    ϕ  = ϕa + ϕp
    ϕ0 = 1- ϕ
    ds = self_diff(ϕ)
    dsp = self_diff_prime(ϕ)
    return 2*(1-ds)*ϕ
end


return_complex_boundary_pt(ϕa; Pe = 10.)
return_complex_boundary([ϕa]; Pe = Pe)
return_complex_boundary([0.1, ϕa]; Pe = Pe)
return_complex_boundary(collect(0.1:0.1:0.2); Pe = Pe)
return_complex_boundary(collect(0.1:0.01:0.2); Pe = Pe)


is_complex_value(0.2, 0.57; Pe = 10.)
is_complex_value(0.2, 0.59; Pe = 10.)
is_complex_value(0.2, 0.0; Pe = 10.)


is_stable_value(0.2, 0.57; Pe = 10.)
is_stable_value(0.99, 0.0; Pe = 10.)

return_stable_boundary_pt(ϕa; Pe = 10.)

return_stable_boundary_pt(0.4366; Pe = 10.)


return_stable_boundary_pt(0.438; Pe = 10.)

ϕa = 0.99999
ϕp = 0.0

function is_complex_value(ϕa, ϕp; Pe = 10.)
    ϕ  = ϕa + ϕp
    ϕ0 = 1- ϕ
    ds = self_diff(ϕ)
    dsp = self_diff_prime(ϕ)
    expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
    return -(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds^2*Pe^2)
end

function is_stable_value(ϕa, ϕp; Pe = 10.)
    ϕ  = ϕa + ϕp
    ϕ0 = 1- ϕ
    ds = self_diff(ϕ)
    dsp = self_diff_prime(ϕ)
    expr1 = (β(ϕa, ϕp) - α(ϕa, ϕp)*Pe^2)
    expr2 = (β(ϕa, ϕp) + α(ϕa, ϕp)*Pe^2)
    return expr1 -4*ϕ + abs(real( sqrt(expr2^2 - 4*ϕp*β(ϕa, ϕp)*ds^2*Pe^2 +0*im)))
end

function return_complex_boundary_pt(ϕa; Pe = 10.)
    f(x) = is_complex_value(ϕa, x; Pe = Pe)
    return find_zeros(f,(0,1-ϕa-1e-8))
end

function return_complex_boundary_outer(ϕas; Pe = 10.)
    high_ϕps = []
    low_ϕps = []
    ϕa_sols = []
    for ϕa in  ϕas
        ϕp = return_complex_boundary_pt(ϕa; Pe = 10.)
        if length(ϕp)==2
                push!(high_ϕps,ϕp[2])
                push!(low_ϕps,ϕp[1])
                push!(ϕa_sols,ϕa)
        end
    end
    return ϕa_sols, low_ϕps, high_ϕps
end

function return_complex_boundary_inner(ϕas; Pe = 10.)
    high_ϕps1 = []
    high_ϕps2 = []
    low_ϕps1 = []
    low_ϕps2 = []
    ϕa_sols = []
    for ϕa in  ϕas
        ϕp = return_complex_boundary_pt(ϕa; Pe = 10.)
        if length(ϕp)==4
                push!(high_ϕps1,ϕp[2])
                push!(low_ϕps1,ϕp[1])
                push!(high_ϕps2,ϕp[4])
                push!(low_ϕps2,ϕp[3])
                push!(ϕa_sols,ϕa)
        end
    end
    return ϕa_sols, low_ϕps1, high_ϕps1, low_ϕps2, high_ϕps2
end

function return_stable_boundary_pt(ϕa; Pe = 10.)
    f(x) = is_stable_value(ϕa, x; Pe = Pe)
    return find_zeros(f,(0,1-ϕa-1e-8))
end


# return_stable_boundary_pt(0.6; Pe = 10.)


function return_stable_boundary_outer(ϕas; Pe = 10.)
    high_ϕps = []
    low_ϕps = []
    ϕa_low = []
    ϕa_high = []
    for ϕa in  ϕas
        ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
        if length(ϕp)>1
            if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)<0)
                push!(low_ϕps,minimum(ϕp))
                push!(ϕa_low,ϕa)
            end
            if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)<0)
                push!(high_ϕps,maximum(ϕp))
                push!(ϕa_high,ϕa)
            end
        elseif length(ϕp)>0
            if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)<0)
                push!(high_ϕps,ϕp[1])
                push!(ϕa_high,ϕa)
            end
        end
    end
    return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary(ϕas; Pe = 10.)
    high_ϕps = []
    low_ϕps = []
    ϕa_low = []
    ϕa_high = []
    for ϕa in  ϕas
        ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
        if length(ϕp)==2
                push!(low_ϕps,minimum(ϕp))
                push!(ϕa_low,ϕa)
                push!(high_ϕps,maximum(ϕp))
                push!(ϕa_high,ϕa)
        elseif length(ϕp)==1
                push!(high_ϕps,ϕp[1])
                push!(ϕa_high,ϕa)
        end
    end
    return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_outer(ϕas; Pe = 10.)
    high_ϕps = []
    low_ϕps = []
    ϕa_low = []
    ϕa_high = []
    for ϕa in  ϕas
        ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
        if length(ϕp)>1
            if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)<0)
                push!(low_ϕps,minimum(ϕp))
                push!(ϕa_low,ϕa)
            end
            if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)<0)
                push!(high_ϕps,maximum(ϕp))
                push!(ϕa_high,ϕa)
            end
        elseif length(ϕp)>0
            if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)<0)
                push!(high_ϕps,ϕp[1])
                push!(ϕa_high,ϕa)
            end
        end
    end
    return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_inner(ϕas; Pe = 10.)
    high_ϕps = []
    low_ϕps = []
    ϕa_low = []
    ϕa_high = []
    for ϕa in  ϕas
        ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
        if length(ϕp)==2
            #if (is_complex_value(ϕa, minimum(ϕp); Pe = Pe)>0)
                push!(low_ϕps,minimum(ϕp))
                push!(ϕa_low,ϕa)
            #end
            #if (is_complex_value(ϕa, maximum(ϕp); Pe = Pe)>0)
                push!(high_ϕps,maximum(ϕp))
                push!(ϕa_high,ϕa)
            #end
        elseif length(ϕp)>0
            #if (is_complex_value(ϕa, ϕp[1]; Pe = Pe)>0)
                push!(high_ϕps,ϕp[1])
                push!(ϕa_high,ϕa)
            #end
        end
    end
    return ϕa_low, ϕa_high, low_ϕps, high_ϕps
end

function return_stable_boundary_extra(ϕas; Pe = 10.)
    high_ϕps = []
    mid_ϕps = []
    low_ϕps = []
    ϕa_low = []
    ϕa_mid = []
    ϕa_high = []
    for ϕa in  ϕas
        ϕp = return_stable_boundary_pt(ϕa; Pe = Pe)
        if length(ϕp)==3
            push!(low_ϕps,minimum(ϕp))
            push!(ϕa_low,ϕa)
            push!(high_ϕps,maximum(ϕp))
            push!(ϕa_high,ϕa)
            push!(mid_ϕps,ϕp[2])
            push!(ϕa_mid,ϕa)
        end
    end
    return ϕa_low, ϕa_mid, ϕa_high, low_ϕps, mid_ϕps, high_ϕps
end



# ϕas = collect(0.01:0.01:0.99)
# ϕa_sols, low_ϕps, high_ϕps = return_complex_boundary(ϕas; Pe = 10.)

# ϕas = collect(0.01:0.01:0.99)
# ϕa_sols2, low_ϕps2, high_ϕps2 = return_stable_boundary(ϕas; Pe = 10.)

