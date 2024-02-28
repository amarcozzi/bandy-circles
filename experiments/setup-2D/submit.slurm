#!/bin/bash

#SBATCH --job-name=fds_setup
#SBATCH --nodes=2
#SBATCH --ntasks=72
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-11%4
#SBATCH --output=%a.out                        
#SBATCH --error=%a.err  

module gcc openmpi
srun ~/fds/Build/fds_ompi_gnu_linux "input_${SLURM_ARRAY_TASK_ID}.fds"
