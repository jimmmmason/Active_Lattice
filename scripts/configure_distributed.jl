cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed

addprocs([("indium",:auto)])

rmprocs([5,6,7])

@distributed for i in 1:nprocs()
    host = gethostname()
    println(host)
end

Distributed.interrupt()