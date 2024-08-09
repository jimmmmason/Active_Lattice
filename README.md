# Active Passive Lattice Gas
This repository contains the code use to produce the data associated with [Dynamical patterns in active-passive particle mixtures with non-reciprocal interactions: Exact hydrodynamic analysis](https://doi.org/10.48550/arXiv.2408.03932)

## Particle simulations

`new_sim_param` will generate a dictionary of parameters compatable with `run_new_sim` and 
`run_new_sim` will run and save a new particle simulation for a given set of parameters. 

## Time dependant solutions

`new_pde_param` will generate a dictionary of parameters compatable with `run_new_pde` and 
`run_new_pde` will run and save a new pde solution for a given set of parameters. 

## Travelling solutions

`solve_full` solves Equations (S40) and (S41) to compute travelling profiles in finite domains.
`solve_out_3`, `solve_out_5`, and `solve_out_6` sovle systems `outer0`, `outer1`, and `outer2` respectively, for travellign profiles in large domains. 

## Phase diagram 

`is_stable_value` returns the real part of the eigenvalue of the matrix (S23) maximised over all wavelengths $q$. 
`colapse_sol_interval` solves for the coexisting phases. It returns `find_sol, lower_limits, upper_limits`. If `find_sol` is `true`, then $\phi_v$ belongs the interval `lower_limits` and $\phi_l$ belongs the interval `upper_limits` and the intervals are smaller than the set tolerance. 
