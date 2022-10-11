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
ssh cyllene
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

rmprocs([2,3,4,5,6])
rmprocs([7])

@distributed for i in 1:nprocs()
    host = gethostname()
    println(host)
end

Distributed.interrupt()