cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using LinearAlgebra

function ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k::Int64 = 10, ω = 2*π,γ= 0.0)
    ρ = ρa + ρp
    v0 = Pe*sqrt(Dθ)
    #
    ds = self_diff(ρ;γ= γ) 
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    #=
    ds = self_diff(ρ)+γ
    dp = self_diff_prime(ρ)
    DD = (1+2*γ-ds)/ρ
    =#
    #=
    # oiginal paramters: 
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    =#
    s = DD - 1
    p = -Dx*ds*ω^2
    q = -v0*im*ω*ds/2 #should have /2 ??? 
    matrix = Complex.(zeros(k+1, k+1))
    for u in 1:(k+1)
        for v in 1:(k+1)
            if abs(u - v) == 1
                matrix[u, v] = q
            elseif u == v
                matrix[u, v] = p - Dθ*(u-2)^2 
            end
        end
    end
    matrix[1, 1] = - Dx*(ω^2)*(ds+ρp*DD)
    matrix[1, 2] = - 2*π*Dx*(ω^2)*ρp*DD
    matrix[1, 3] = - π*v0*im*ω*ρp*s
    matrix[2, 1] = - Dx*(ω^2)*(ρa/(2*π))*DD
    matrix[2, 2] = - Dx*(ω^2)*(ds+ρa*DD)
    matrix[2, 3] = - v0*im*ω*(ρa*s+ds)/2
    matrix[3, 1] = - v0*im*ω*(ρa/(2*π))*dp
    matrix[3, 2] = - v0*im*ω*(ρa*dp+ds)
    return matrix
end

function ap_MathieuEigen(matrix)
    Egs = eigen(matrix)
    return Egs.values, Egs.vectors
end 

function ap_MathieuEigen_lite(matrix; k=20)
    return  eigvals(matrix)[k+1]
end

function ap_lin_stab_line(ρa,ρp; Dx =1. ,Pe = 50., Dθ = 100., k = 40 )
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return real(amin)
end

function ap_lin_stab_imaginary(ρa,ρp; Dx =1. ,Pe = 1., Dθ = 100., k = 40 )
    matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end

###
###


function a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k::Int64 = 10, ω = 2*π,γ= 0.0)
    ρ = ρa + ρp
    v0 = Pe*sqrt(Dθ)
    #
    ds = self_diff(ρ;γ= γ) 
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    #=
    ds = self_diff(ρ)+γ
    dp = self_diff_prime(ρ)
    DD = (1+2*γ-ds)/ρ
    =#
    #=
    # oiginal paramters: 
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    =#
    s = DD - 1
    p = -Dx*ds*ω^2
    q = -v0*im*ω*ds/2 #should have /2 ??? 
    matrix = Complex.(zeros(k, k))
    for u in 1:(k)
        for v in 1:(k)
            if abs(u - v) == 1
                matrix[u, v] = q
            elseif u == v
                matrix[u, v] = p - Dθ*(u-1)^2 
            end
        end
    end

    matrix[1, 1] = - Dx*(ω^2)*(ds + ρa*DD)
    matrix[1, 2] = - v0*im*ω*(ρa*s + ds)/2
    matrix[2, 1] = - v0*im*ω*(ρa*dp+ds)
    return matrix
end

function a_MathieuEigen(matrix)
    Egs = eigen(matrix)
    return Egs.values, Egs.vectors
end 

function a_MathieuEigen_lite(matrix; k=20)
    return  eigvals(matrix)[k]
end

function a_lin_stab_line(ρa,ρp; Dx =1. ,Pe = 50., Dθ = 100., k = 40 )
    matrix = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = a_MathieuEigen_lite(matrix; k = k)  
    return real(amin)
end

function a_lin_stab_imaginary(ρa,ρp; Dx =1. ,Pe = 1., Dθ = 100., k = 40 )
    matrix = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = a_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end







###
###

function ap_MathieuMatrix_checker(ρa,ρp,Dx,Pe,Dθ; k::Int64 = 10, ω = [2π, 0],γ= 0.0)
    ρ = ρa + ρp
    v0 = Pe*sqrt(Dθ)
    #
    ds = self_diff(ρ;γ= γ)
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    #=
    ds = self_diff(ρ)+γ
    dp = self_diff_prime(ρ)
    DD = (1+2*γ-ds)/ρ
    =#
    #=
    # oiginal paramters: 
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    DD = (1-ds)/ρ
    =#
    s = DD - 1
    p = -Dx*ds*norm(ω)^2
    q = -v0*im*ω*ds/2
    α = -Dx*DD*[ρp,ρa]*norm(ω)^2
    β = -v0*im*ω*s/2
    γ = -v0*ρa*dp*im*ω
    matrix = Complex.(zeros(2*k, 2*k))
    for u in 1:(2*k)
        for v in 1:(2*k)
                if abs(u - v) == 2
                    matrix[u, v] = q[1]
                elseif v - u == 1
                    matrix[u, v] = q[2]
                elseif v - u == -1
                        matrix[u, v] = -q[2]
                elseif u == v
                    matrix[u, v] = p - Dθ*( floor(u/2))^2 
                end
        end
    end

    matrix[1, 1] = p + α[1]
    matrix[1, 2] = α[1]
    matrix[1, 3] = β[1]*ρp
    matrix[1, 4] = β[2]*ρp

    matrix[2, 1] = α[2]
    matrix[2, 2] = p + α[2]
    matrix[2, 3] = β[1]*ρa +q[1]
    matrix[2, 4] = β[2]*ρa +q[2]

    matrix[3, 1] = γ[1]
    matrix[3, 2] = γ[1] + 2*q[1]
    matrix[3, 3] = p - Dθ
    matrix[3, 4] = 0

    matrix[4, 1] = γ[2]
    matrix[4, 2] = γ[2] + 2*q[2]
    matrix[4, 3] = 0
    matrix[4, 4] = p - Dθ

    return matrix
end

function ap_MathieuEigen_checker(matrix)
    Egs = eigen(matrix)
    return Egs.values, Egs.vectors
end 

function ap_MathieuEigen_lite_checker(matrix; k=20)
    return  eigvals(matrix)[4]
end

function ap_lin_stab_line_checker(ρa,ρp; Dx =1. ,Pe = 50., Dθ = 100., k = 40 )
    matrix = ap_MathieuMatrix_checker(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite_checker(matrix; k = k)  
    return real(amin)
end

function ap_lin_stab_imaginary_checker(ρa,ρp; Dx =1. ,Pe = 1., Dθ = 100., k = 40 )
    matrix = ap_MathieuMatrix_checker(ρa,ρp,Dx,Pe,Dθ; k = k);
    amin = ap_MathieuEigen_lite_checker(matrix; k = k)  
    return imag(amin)
end

function lin_stab_line_fraction_checker(ρ,χ; Dx =1. ,Pe = 20., Dθ =100.,k=40, ω = [2*π,0] )
    ρa = χ*ρ
    ρp = (1-χ)*ρ
    matrix = ap_MathieuMatrix_checker(ρa,ρp,Dx,Pe,Dθ; k = k, ω = ω);
    amin = ap_MathieuEigen_lite_checker(matrix; k = k)  
    return real(amin)
end




###


#=
ρ = 0.7
Pe = 3.
Dθ = 100.
Dx = 1.

ρa = 0.6
ρp = 0.

m1 = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k =2)
v = eigvals(m1)[2]
_, V = eigen(m1)
v1 = V[:,2]

v1' *m1 *v1 / (v1' *v1)

v2 = push!(v1,0.)
m2 = a_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k =3)
eigvals(m2)

v2' *m2 *v2 / (v2' *v2)

@unpack Dθ, Dx, ρp, ρa, Pe = param
k = 10
matrix = ap_MathieuMatrix(ρa,ρp,Dx,Pe,Dθ; k = k, ω = -2*π);
a, A = ap_MathieuEigen(matrix)
hcat(A[:,k+1],A[:,k])

hcat(A[:,k+1]/A[2,k+1],A[:,k]/A[2,k])


matrix = [ 10 + im 0 0; 0 -1 0 ;  0 0 2]
eigvals(matrix)




#Mparams = stabparams(param)

ω2 = (2 * pi)^2
ds = self_diff(ρ)
c1 = Dx * ds
q,m1,m2,γ = mstabparams_lite(ρa,ρ,Dx,Pe,Dθ)
matrix = MathieuMatrix(q,m1,m2; k = 5)
a, A = MathieuEigen(matrix, γ) 




lin_stab_line(ρ; Dx =Dx ,Pe = Pe, Dθ = Dθ)

x = 0.01:0.01:0.99
y = 0.1:0.1:10
fig, ax = PyPlot.subplots(figsize =(10, 10))

n = length(x)
m = length(y)
z = zeros(m,n)
    for i in 1:m, j in 1:n
        z[i,j] = lin_stab_line(x[j]; Dx =Dx ,Pe = y[i], Dθ = Dθ)
    end
ax.contour(x,y,z; levels = [0])
display(fig)

Real(im)


matrix = MathieuMatrix(q,m1,m2; k = 5);
a, A = MathieuEigen(matrix, 0.5);
=#
