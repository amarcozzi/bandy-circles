#!/bin/bash

#SBATCH -J mirror-rad-2
#SBATCH --nodes=4
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH -A umontana_fire_modeling
#SBATCH -o a.out

# Run FDS
module load intel-oneapi-mpi
srun ~/fds/Build/fds_impi_intel_linux input_0.fds
