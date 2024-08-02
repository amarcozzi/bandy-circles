#!/bin/bash

#SBATCH -J test_domain        			# Job name
#SBATCH -o test.out
#SBATCH --nodes=2-10
#SBATCH -n 76
#SBATCH --mem-per-cpu=1G
#SBATCH -A umontana_fire_modeling

# Run FDS
module load gcc openmpi
srun ~/fds/Build/fds_ompi_gnu_linux input_0.fds