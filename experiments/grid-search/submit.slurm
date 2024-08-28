#!/bin/bash

#SBATCH --job-name=grid-search
#SBATCH --partition=cpu(all)
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-59%4                           
#SBATCH --nodes=4
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --time=24:00:00
#SBATCH --output=simulations/%a.out                        
#SBATCH --error=simulations/%a.err  

echo "Starting job ${SLURM_ARRAY_TASK_ID}"

# Activate conda environment
echo "Activating conda environment"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bandy-circles

# Generate the FDS input file
echo "Generating FDS input file"
python3 ../../src/preprocess_from_csv.py $PWD $SLURM_ARRAY_TASK_ID

# Run FDS
echo "Running FDS"
cd simulations/simulation_$SLURM_ARRAY_TASK_ID
mpirun -n 54 fds "input_${SLURM_ARRAY_TASK_ID}.fds"