#!/bin/bash

#SBATCH --job-name=$job_name
#SBATCH --account=umontana_fire_modeling
#SBATCH --nodes=2
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --time=6:00:00
#SBATCH --output=$job_name_output                      
#SBATCH --error=$job_name_error


# Activate conda environment
source /project/umontana_fire_modeling/anthony.marcozzi/miniforge3/etc/profile.d/conda.sh
conda activate bandy-circles

# Activate MPI Module
module load intel-oneapi-mpi

# Run FDS
srun ~/fds/Build/fds_impi_intel_linux $fds_input_file_name
