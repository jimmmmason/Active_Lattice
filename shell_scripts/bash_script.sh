#!/bin/bash

export PATH="/home/jm2386/julia-1.8.1/bin:${PATH}"
# nice -n 19 julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/pm_slurm_sweeper.jl $SLURM_ARRAY_TASK_ID &
nice -n 19 julia /home/jm2386/Active_Lattice/pm_scripts/pm_pde_scripts/pm_slurm_test.jl $SLURM_ARRAY_TASK_ID &