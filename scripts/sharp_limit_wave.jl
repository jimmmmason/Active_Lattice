cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
include("/home/jm2386/Active_Lattice/src/article_src.jl")
#



#
using TensorKit, LinearAlgebra
function tensor_equations(m_max,k_max ; param = param, wave_speed =0.)
    @unpack ds, ϕₚ, ϕₐ, DD, Pe, alpha, beta, gamma = param
    tensor = zeros(m_max,k_max,m_max+2,k_max+1)
    c = wave_speed
    for m in 1:m_max
        #
            # Equation 1 coefficients
            #
            tensor[m, 1, m+2, 1] = sqrt((m+2)*(m+1))/2 * (ds + ϕₚ * DD)
            tensor[m, 1, m+2, 2] = sqrt((m+2)*(m+1))/2 * ϕₚ * DD      
            tensor[m, 1, m, 1] = (2m+1) * (ds + ϕₚ * DD)
            tensor[m, 1, m, 2] = (2m+1) * ϕₚ * DD
            #
            tensor[m, 1, m+1, 1] = sqrt((m+1)/2) * c
            tensor[m, 1, m+1, 3] = sqrt((m+1)/2) * ϕₚ * beta
            #
            # Equation 2 coefficients
            #
            tensor[m, 2, m+2, 1] = sqrt((m+2)*(m+1))/2 * (ϕₐ * DD)
            tensor[m, 2, m+2, 2] = sqrt((m+2)*(m+1))/2 * (ds + ϕₐ * DD)
            tensor[m, 2, m, 1] = (2m+1) * (ϕₐ * DD)
            tensor[m, 2, m, 2] = (2m+1) * (ds + ϕₐ * DD)
            #
            tensor[m, 2, m+1, 2] = sqrt((m+1)/2) * c
            tensor[m, 2, m+1, 3] = sqrt((m+1)/2) * (ϕₐ * beta + alpha)
            #
            # Equation 3 coefficients
            #
            tensor[m, 3, m+2, 3] = sqrt((m+2)*(m+1))/2 * ds
            tensor[m, 3, m, 3] = (2m+1) * ds -1

            #
            tensor[m, 3, m+1, 1] = sqrt((m+1)/2) * (gamma)
            tensor[m, 3, m+1, 2] = sqrt((m+1)/2) * (gamma + 2*alpha)
            tensor[m, 3, m+1, 3] = sqrt((m+1)/2) * c
            tensor[m, 3, m+1, 3] = sqrt((m+1)/2) * (alpha)      
            #

        if m >2 
            tensor[m, 1, m-2, 1] = sqrt(m*(m-1))/2 * (ds + ϕₚ * DD)
            tensor[m, 1, m-2, 2] = sqrt(m*(m-1))/2 * (ds + ϕₚ * DD)

            tensor[m, 2, m-2, 1] = sqrt(m*(m-1))/2 * (ϕₐ * DD)
            tensor[m, 2, m-2, 2] = sqrt(m*(m-1))/2 * (ds + ϕₐ * DD)

            tensor[m, 3, m-2, 3] = sqrt(m*(m-1))/2 * ds
        end

        if m >1
            tensor[m, 1, m-1, 1] = sqrt(m/2) * c
            tensor[m, 1, m-1, 3] = sqrt(m/2) * ϕₚ * beta

            tensor[m, 2, m-1, 2] = sqrt(m/2) * c
            tensor[m, 2, m-1, 3] = sqrt((m+1)/2) * (ϕₐ * beta + alpha)

            tensor[m, 3, m-1, 1] = sqrt((m)/2) * (gamma)
            tensor[m, 3, m-1, 2] = sqrt((m)/2) * (gamma + 2*alpha)
            tensor[m, 3, m-1, 3] = sqrt((m)/2) * c
            tensor[m, 3, m-1, 3] = sqrt((m)/2) * (alpha)
        end


        #
        for k in 1:k_max
            #
            # Equation k coefficients
            #
            if k >= 4
                #
                tensor[m, k, m+2, k] = sqrt((m+2)*(m+1))/2 * ds
                tensor[m, k, m, k] = (2m+1) * ds - (k-2)^2
                if m >2
                    tensor[m, k, m-2, k] = sqrt(m*(m-1))/2 * ds
                end
                #
                tensor[m, k, m+1, k-1] = sqrt((m+1)/2) * (alpha)
                tensor[m, k, m+1, k] = sqrt((m+1)/2) * c
                tensor[m, k, m+1, k+1] = sqrt((m+1)/2) * (alpha)
                if m >1
                    tensor[m, k, m-1, k-1] = sqrt((m)/2) * (alpha)
                    tensor[m, k, m-1, k] = sqrt((m)/2) * c
                    tensor[m, k, m-1, k+1] = sqrt((m)/2) * (alpha)
                end
            end
        end
    end
    return tensor[:,:,1:m_max,1:k_max]
end
function add_parameters(param)
    @unpack ρa, ρp, Pe = param
    ϕ = ρa+ρp
    ds = self_diff(ϕ)
    dsp = self_diff_prime(ϕ)
    ϕₚ  = ρp
    ϕₐ = ρa
    DD = (1-ds)/ϕ
    s = DD- 1
    alpha = Pe*ds/2
    beta = -Pe*s/2
    gamma = -Pe*ϕₐ*dsp
    @pack! param = ds, ϕₚ, ϕₐ, DD, Pe, alpha, beta, gamma
    return param
end
function tri_inv( m_max ; )
    matrix = zeros(m_max,m_max)
    for m in 2:m_max
            matrix[m-1,m] = sqrt(m)/2
            matrix[m,m-1] = sqrt(m-1)/2
    end
    return inv(matrix)
end
function tensor_eigen_problem(m_max, k_max; param = param)
    equations = tensor_equations(m_max, k_max; param = param)
    matrix_inv = tri_inv(m_max)
    tensor = zeros(m_max, k_max, m_max, k_max)
    for k in 1:k_max
        for l in 1:k_max
        tensor[:,k,:,l] = matrix_inv*equations[:,k,:,l]
        end
    end
    return tensor
end
#

#
params = []
    χ = 0.1
    Dθ = 1000.0
    Dx = 1.
    ρ_start = 0.4
    ρ_end = 1.0
    Pe_end = 50
    Pe_start = 0
    xs = collect(0.4:0.001:0.999)

    using Roots
    k= 5
    X=[]
    Y=[]
    for x in xs
            try
                local f
                f(y) = lin_stab_line_fraction(x,χ; Dx =Dx ,Pe = y, Dθ = Dθ, k = k)
                Pe = find_zero(f, (0.,  100.))
                push!(Y,Pe)
                push!(X,x)
            catch
            end
    end

    Pemin = minimum(Y)
    i = argmin(Y)
    ϕmin = X[i]

    params = []
    pert = "n=1"
    ρ = ϕmin
    Pe = Pemin
    
    # χ = 0.1
    # Dθ = 4.0
    # ϕmin = 0.929
    # Pemin = 41.1445070772628
    # lin_imaginary_fraction(ρ,χ; Dx =Dx ,Pe = Pe, Dθ = Dθ, k = k)
    # lin_stab_line_fraction(ρ,χ; Dx =Dx ,Pe = Pe, Dθ = Dθ, k = k)

    T  = 24.0
    save_interval = 0.01
    δ  = 1e-3
    k = 20
    Nx = 128
    Nθ = 64
    name = "wave_pert+rand_2d_δ=$(δ)"

    param = pde_param_fraction(; name = name, 
                    ρ = ρ, Pe = Pe, χ = χ, T = T, 
                    Dθ = Dθ, δt = 1e-5, Nx = Nx, Nθ = Nθ, 
                    save_interval = save_interval, max_steps = 1e8,
                    pert = pert, δ = δ, k = k, video_length = 80000.,
                    cbar_max = 1.0, cbar_min = 0.8
    )
    params = [param]
#


#
tparam = add_parameters(param)
m = 1000
k = 10
m_max = m
k_max = k
tensor = tensor_eigen_problem(m, k; param = tparam)
tensor_map = TensorMap(tensor, (ℂ^m ⊗ ℂ^k) , (ℂ^m ⊗ ℂ^k))

D, V = eig(tensor_map)

DA = convert(Array,D)
DA = diag(DA)

i = argmin(abs.(DA))
DA[i]



VA = convert(Array,V)
VA[1:1000,:,i]

#i = argmin(sum( abs.( VA[1000,:,:]) ; dims = 1  )[1,:]   )

j = argmax(abs.(VA[:,:,i]))
j = argmax(abs.(VA[:,1,i]))



using FastTransforms
using ApproxFun

using Pkg
Pkg.add("FastTransforms")



f = Fun(Hermite(), [0.,1.])

x = -1000:1:1000
y = f.(x)
fig , ax = plt.subplots(1, 1, figsize=(10,10))
ax.plot(x,y)
dsipaly(fig)


using UpdateJulia
