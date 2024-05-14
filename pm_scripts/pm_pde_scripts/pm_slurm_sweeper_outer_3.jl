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

# outer_3 starter h sweep
function load_adjacent_3(Lx, Nx, ϕa, ϕp, v0) # f,u,c = / "fail"
    ind_range = [(i,j) for j in [0,1,-1], i in [0,1,-1,2,-2]]
    input_params = [(Lx, Nx, d2(ϕa + 0.01*i), d2(ϕp +0.01*j), v0) for (i,j) in ind_range]
    ϵ = 1e-4
    for input_param in input_params
        try
            Lx, Nx, ϕa, ϕp, v0 = input_param
            f,u,c = load_out_3(input_param...)
            normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
            param, p2 = get_outer_param_3(input_param...)
            err, avmag, c = check_u_3(u,p2);
            if (normf > ϵ)&&(err<ϵ)
                print("loaded: $(input_param)")
                return f,u,c
            end
        catch
        end
    end
    return "fail"
end
# start point
# start point
ind_range = [(i,j) for j in [0], i in collect(-10:1:10)] #21 
ϕa, ϕp = [ (d2(0.45 + 0.01*i), d2(0.35 +0.01*j)) for (i,j) in ind_range][input]
Lx, Nx, _, _, v0 = 100.0, 1024, 0.45, 0.35, 7.5
# sweep range
ind_range = [(i,i) for i in vcat(collect(1:10),collect(-1:-1:-20)) ]
input_params = [(Lx, Nx, d2(ϕa + (-0.36/25)*i), d2(ϕp +0.01*j), v0) for (i,j) in ind_range]
ϵ = 1e-4
# load adjecnt then solve; stop if err/small 
for input_param in input_params
    output = load_adjacent_3(input_param...);
    if output == "fail"
        print("load_fail: $(input_param)")
        break 
    else
        f,u,c = output;
        Lx, Nx, ϕa, ϕp, v0 = input_param
        f,u,c = solve_out_3(Lx,Nx,ϕa,ϕp,v0,u; tol = 1e-8)
        param, ps = get_outer_param_3(input_param...)
        err, avmag, c = check_u_3(u,ps);
        normf = sqrt(sum( (f[:,1] .- ϕa/2).^2 + (f[:,2] .- ϕa/2).^2 + (f[:,3] .- ϕp).^2)/Nx)
        if (normf < ϵ)|(err > ϵ)
            print("end: $(input_param)")
            break
        end
    end
end