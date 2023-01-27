cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed

#= dumbo bash commands to check they are awake
ssh jm2386@ssh.maths.cam.ac.uk

wake cherry
wake cyllene
wake orthosie
wake refract

ssh cherry
exit
##ssh cyllene
exit
ssh orthosie
exit
ssh refract
exit
=#
#addprocs([("cherry",:auto)])
#addprocs([("cyllene",:auto)])
#addprocs([("indium",:auto)])
#addprocs([("orthosie",:auto)])
addprocs([("refract",:auto)])
addprocs([("radius",:auto)])

rmprocs()
rmprocs([8,12,14])


@distributed for i in 1:nworkers()
    host = gethostname()
    println(host)
end
nworkers()

Distributed.interrupt()

main_pool = WorkerPool(collect(2:19))