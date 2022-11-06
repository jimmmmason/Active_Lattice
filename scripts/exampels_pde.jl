### Dθ = 10 stability 
name = "high_density_stability_v4"
Nx = 50
Nθ = 20
filename = "/store/DAMTP/jm2386/Active_Lattice/data/pde_pro/$(name)/stability_Nx=$(Nx)_Nθ=$(Nθ).jld2"
#wsave(filename,stabdata)
stabdata = wload(filename)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.05:0.05:0.95 ,xs = collect(0.01:0.01:0.99), xtic = 0.0:0.1:1, ytic = 0:10:100, axlim = [0., 1., 15., 100.])
display(fig)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.9:0.01:0.99 ,xs = collect(0.01:0.01:0.99), xtic = 0.9:0.1:1, ytic = 0:10:100, axlim = [0.9, 1., 15., 100.])
display(fig)
### Dθ = 100 stability 
λs = 10.:10.:200.
Dθ = 100.
name = "high_density_stability_v5"
stabdata = refresh_stab_data(; ρs = 0.05:0.05:1.0,   Dθ = Dθ, Nx = 50, Nθ = 20, λs = λs, name = name)
fig, ax = PyPlot.subplots(figsize =(10, 10))
plot_stab(fig, ax, stabdata; ρs = 0.05:0.05:0.95 ,xs = collect(0.01:0.01:0.99), xtic = 0:0.2:1, ytic = λs , axlim = [0., 1., 30., 200.], Dθ = Dθ)
display(fig)

refresh_stab_data(; ρs = [0.6], Dθ = Dθ, Nx = 50, Nθ = 20, λs = λs, name = name)



using SymPy
const sympy_parsing_mathematica = SymPy.PyCall.pyimport("sympy.parsing.mathematica")
mathematica2julia(s::AbstractString, substitutions::Pair{<:AbstractString,<:AbstractString}...) = SymPy.lambdify(sympy_parsing_mathematica["mathematica"](s, Dict(substitutions...)))

SymPy.lambdify()

f = sympy_parsing_mathematica["mathematica"](" Sqrt[((-2 + pi)/((-1 + pi) (26 + (-11 + pi) pi)))]x")
g = lambdify(f)
g(1)

x = 1
g.([1])

mathematica2julia("Log[Log[x]]")
mathematica2julia(" Sqrt[((-2 + pi)/((-1 + pi) (26 + (-11 + pi) pi)))]x")