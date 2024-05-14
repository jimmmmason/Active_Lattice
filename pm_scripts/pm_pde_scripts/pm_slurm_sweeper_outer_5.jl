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

# outer_5 starter ϕ sweep
function load_adjacent_5(Lx, Nx, ϕa, ϕp, v0, ϕ, γ) # f,u,c = / "fail"
    ind_range = [(i,j) for j in [0,1,-1], i in [0,1,-1,2,-2]]
    input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
    ϵ = 1e-4
    for input_param in input_params
        try
            Lx,Nx,ϕa,ϕp,v0,ϕ,γ = input_param
            param, ps = get_outer_param_5(input_param...)
            f,u,c = load_out_5(input_param...)
            ϕp = sum(f)/Nx-sum(f[:,1:2])/Nx
            ϕa = sum(f[:,1:2])/Nx
            normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
            err, erri, avmag, cep = check_u_5(u,ps)
            if (normf > ϵ)&&(err<ϵ)
                print("loaded: $((Lx,Nx,d2(ϕa),d2(ϕp),v0,ϕ,γ))")
                return f,u,c
            end
        catch
        end
    end
    return "fail"
end
# start point
Lx,Nx,ϕa,ϕp,v0,ϕ,_ = 100.0, 1024, 0.3, 0.3, 7.5, 0.56, 2.2
γ = collect(1.26:0.01:2.2)[input] # 95
# sweep range
ind_range = [(i,0) for i in collect(-1:-1:-18)]
input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
ϵ = 1e-4
# load adjecnt then solve; stop if err/small 
for input_param in input_params
    output = load_adjacent_5(input_param...);
    if output == "fail"
        print("load_fail: $(input_param)")
        break 
    else
        f,u,c = output;
        Lx,Nx,ϕa,ϕp,v0,ϕ,γ = input_param
        
        f,u,c = solve_out_5(Lx,Nx,ϕa,ϕp,v0,u; tol = 1e-8, maxiters = 10)
        ϕp = d2(sum(f)/Nx-sum(f[:,1:2])/Nx)
        ϕa = d2(sum(f[:,1:2])/Nx)

        param, ps = get_outer_param_5(Lx,Nx,ϕa,ϕp,v0,ϕ,γ)
        err, erri, avmag, cep = check_u_5(u,ps)
        normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)

        if (normf < ϵ)
            print("norm end: $(input_param)")
            break
        elseif (err > ϵ)
            print("err end: $(input_param)")
            break
        end
        ax.scatter(ϕa, ϕp)
    end
end
# sweep range
ind_range = [(i,0) for i in collect(1:41) ] 
input_params = [(Lx,Nx,ϕa,ϕp,v0, d2(ϕ+0.01*i),d2(γ+0.01*j) ) for (i,j) in ind_range]
ϵ = 1e-4
# load adjecnt then solve; stop if err/small 
for input_param in input_params
    output = load_adjacent_5(input_param...);
    if output == "fail"
        print("load_fail: $(input_param)")
        break 
    else
        f,u,c = output;
        Lx,Nx,ϕa,ϕp,v0,ϕ,γ = input_param
        
        f,u,c = solve_out_5(Lx,Nx,ϕa,ϕp,v0,u; tol = 1e-8, maxiters = 10)
        ϕp = d2(sum(f)/Nx-sum(f[:,1:2])/Nx)
        ϕa = d2(sum(f[:,1:2])/Nx)

        param, ps = get_outer_param_5(Lx,Nx,ϕa,ϕp,v0,ϕ,γ)
        err, erri, avmag, cep = check_u_5(u,ps)
        normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)

        if (normf < ϵ)
            print("norm end: $(input_param)")
            break
        elseif (err > ϵ)
            print("err end: $(input_param)")
            break
        end
        ax.scatter(ϕa, ϕp)
    end
end