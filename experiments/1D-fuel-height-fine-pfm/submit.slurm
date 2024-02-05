#!/bin/bash

#SBATCH --job-name=fuel_height
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=72:00:00
#SBATCH --partition=cpu(all)
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-11%12                             
#SBATCH --output=output/%a.out                        
#SBATCH --error=output/%a.err  

echo "Starting job ${SLURM_ARRAY_TASK_ID}"

# Create a directory for the output files
mkdir output/simulation_$SLURM_ARRAY_TASK_ID

# Activate conda environment
echo "Activating conda environment"
source ~/miniconda3/etc/profile.d/conda.sh

# Generate the FDS input file
echo "Generating FDS input file"
conda activate bandy-circles
python3 ../../src/preprocess_from_csv.py $PWD $SLURM_ARRAY_TASK_ID

# Run FDS
echo "Running FDS"
cd output/simulation_$SLURM_ARRAY_TASK_ID
mpirun -n 12 fds "input_${SLURM_ARRAY_TASK_ID}.fds"