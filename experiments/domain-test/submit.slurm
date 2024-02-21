#!/bin/bash

#SBATCH --job-name=fds_domain
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH --partition=cpu(all)
#SBATCH --account=umontana_fire_modeling

# Run FDS
mpirun -n 48 fds input_0.fds
