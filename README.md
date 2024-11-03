# Active Lattice Gas
This repository contains the code used to produce the data associated with [Exact hydrodynamics and onset of phase separation for an active exclusion process](https://royalsocietypublishing.org/doi/full/10.1098/rspa.2023.0524)

## Particle simulations

`uniform_initial_param` will generate a dictionary of parameters compatible with `run_sim` and 
`run_sim_dump` will run and save a new particle simulation for a given set of parameters. 

## Time dependant solutions

`pde_param_k` will generate a dictionary of parameters compatible with
`perturb_pde_run` which will perturb the homogeneous steady state, then run and save a new pde solution for a given set of parameters. See `scripts/run_pdes.jl` for examples. 

## Linear stability

`lin_stab_line` will compute the real part of the largest eigenvlaue in (3.7). 

