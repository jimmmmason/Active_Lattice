cd("/home/jm2386/Active_Lattice/")
using DrWatson
@quickactivate "Active_Lattice"
using Distributed

addprocs([("indium",:auto)])
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
addprocs([("cherry",:auto)])
addprocs([("cyllene",:auto)])
addprocs([("orthosie",:auto)])
addprocs([("refract",:auto)])
addprocs([("radius",:auto)])

rmprocs([4])
rmprocs([8,12,14])

@distributed for i in 1:nprocs()
    host = gethostname()
    println(host)
end

Distributed.interrupt()

main_pool = WorkerPool([18, 24, 26, 27, 25, 20, 21, 19, 22, 28, 31, 29])