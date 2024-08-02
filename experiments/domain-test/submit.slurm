#!/bin/bash

#SBATCH --job-name=fds_domain
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=76
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition=atlas
#SBATCH --account=umontana_fire_modeling

# Run FDS
module gcc openmpi
srun ~/fds/Build/fds_ompi_gnu_linux input_0.fds