cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"

using LinearAlgebra

function stabparams(param)
    @unpack name, L, D, λ, γ, ρa, ρp, Dθ = param
    ρ = ρa + ρp
    v0 = Pe*sqrt(Dθ)
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    c1 = Dx * ds
    c2 = Dx * (1 - ds) / (2 * pi)
    c3 = -v0 * (1 -  ρ - ds)
    c4 = -v0 * (dp *  ρ) / 2
    c5 = -v0 * (ds) / 2
    c6 = Dθ
    ω = 2 * pi
    ω2 = (2 * pi)^2
    #α = pi / 2 - 1
    q = -im * ω * c5 / c6
    #s = -ω * c5 / c6
    m1 = c2 * ω2 / c6
    m2 = -ω * sqrt(-(c4 + c5) * (c3 + c5)) / c6
    γ = sqrt(-(c4 + c5) / (c3 + c5))*im
    Mparam = Dict{String,Any}()
    @pack! Mparam = c1, c2, c3, c4, c5, c6, γ, m1, m2, q
    return Mparam
end

function mstabparams_lite(ρa,ρ,Dx,Pe,Dθ)
    v0 = Pe*sqrt(Dθ)
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    c1 = Dx * ds
    c2 = Dx * ρa *(1 - ds) / (ρ) #not sure about 2pi 2 * pi * 
    c3 = -v0 * ρa * (1 -  ρ - ds) / ρ
    c4 = -v0 * (dp *  ρa ) / 2
    c5 = -v0 * (ds) / 2
    c6 = Dθ
    ω = 2 * pi
    ω2 = (2 * pi)^2
    #α = pi / 2 - 1
    q = -im * ω * c5 / c6
    #s = -ω * c5 / c6
    m1 = c2 * ω2 / c6
    m2 = -ω * sqrt(complex(-(c4 + c5) * (c3 + c5))) / c6
    γ = sqrt(complex(-(c4 + c5) / (c3 + c5)))*im
    return Complex(q),Complex(m1),Complex(m2),Complex(γ)
end

function MathieuMatrix(q::Complex{Float64}, m1::Complex{Float64}, m2::Complex{Float64}; k::Int64 = 10)
    matrix = Complex.(zeros(k, k))
    for u in 1:k
        for v in 1:k
            if abs(u - v) == 1
                matrix[u, v] = q
            elseif u == v
                matrix[u, v] = (u-1)^2
            end
        end
    end
    matrix[2, 1] = m2
    matrix[1, 2] = m2
    matrix[1, 1] = m1
    return matrix
end

function MathieuEigen(matrix, γ::Complex{Float64})
    Egs = eigen(matrix)
    Egs.vectors[1,:] = Egs.vectors[1,:]/γ
    return Egs.values, Egs.vectors
end 

function MathieuEigen_lite(matrix)
    return  eigvals(matrix)[1]
end

function lin_stab_line(ρ; Dx =1. ,Pe = 1., Dθ = 100. )
    ω2 = (2 * pi)^2
    ds = self_diff(ρ)
    c1 = Dx * ds
    q,m1,m2,γ = mstabparams_lite(ρ,ρ,Dx,Pe,Dθ)
    matrix = MathieuMatrix(q,m1,m2; k = 20);
    amin = MathieuEigen_lite(matrix)  
    return real(amin + c1*ω2/Dθ)
end

##


function ap_mstabparams_lite(ρa,ρp,Dx,Pe,Dθ)
    ρ = ρa + ρp
    v0 = Pe*sqrt(Dθ)
    ds = self_diff(ρ)
    dp = self_diff_prime(ρ)
    c = fill(Complex(0.),9)
    c[1] = Dx * ds
    c[2] = Dx * ρa *(1 - ds) / (ρ) #not sure about 2pi 2 * pi * 
    c[3] = -v0 * ρa * (1 -  ρ - ds) / (2* ρ)
    c[4] = -v0 * (dp *  ρa ) / 2
    c[5] = -v0 * (ds) / 2
    c[6] = Dθ
    c[7] = Dx * ds
    c[8] = Dx * ρp *(1 - ds) / (ρ)
    c[9] = -v0 * ρp * (1 -  ρ - ds) / (2* ρ)
    ω = 2 * pi
    return c, ω
end

function ap_MathieuMatrix(c, ω; k::Int64 = 10)
    matrix = Complex.(zeros(k+1, k+1))
    for u in 1:(k+1)
        for v in 1:(k+1)
            if abs(u - v) == 1
                matrix[u, v] = im * c[5] * ω
            elseif u == v
                matrix[u, v] = - c[1]*ω^2-c[6]*(u-2)^2 
            end
        end
    end
    matrix[1, 1] = -(c[7]+ (2*π)*c[8])*ω^2
    matrix[1, 2] = - c[8]*ω^2
    matrix[1, 3] = c[9]*im*ω
    matrix[2, 1] = - c[2]*ω^2/(2*π)
    matrix[2, 2] = -(c[1]+ c[2])*ω^2
    matrix[2, 3] = (c[3]+ c[5])*im*ω
    matrix[3, 1] = (c[4])*im*ω/(2*π)
    matrix[3, 2] = (c[4]+ c[5])*im*ω
    return matrix
end

function ap_MathieuEigen(matrix)
    Egs = eigen(matrix)
    return Egs.values, Egs.vectors
end 

function ap_MathieuEigen_lite(matrix; k=20)
    return  eigvals(matrix)[k+1]
end

function ap_lin_stab_line(ρa,ρp; Dx =1. ,Pe = 1., Dθ = 100., k = 20 )
    c,ω = ap_mstabparams_lite(ρa,ρp,Dx,Pe,Dθ)
    matrix = ap_MathieuMatrix(c,ω; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return real(amin)
end

function ap_lin_stab_imaginary(ρa,ρp; Dx =1. ,Pe = 1., Dθ = 100., k = 20 )
    c,ω = ap_mstabparams_lite(ρa,ρp,Dx,Pe,Dθ)
    matrix = ap_MathieuMatrix(c,ω; k = k);
    amin = ap_MathieuEigen_lite(matrix; k = k)  
    return imag(amin)
end

#=
ρ = 0.7
Pe = 3.
Dθ = 100.
Dx = 1.
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
