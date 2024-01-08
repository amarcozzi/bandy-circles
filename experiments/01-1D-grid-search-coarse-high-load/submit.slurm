#!/bin/bash

#SBATCH --job-name=1d_parameter_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=24:00:00
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-99%100                             
#SBATCH --output=output/%a.out                        
#SBATCH --error=output/%a.err  

echo "Starting job ${SLURM_ARRAY_TASK_ID}"

# Create a directory for the output files
cd output/
mkdir simulation_$SLURM_ARRAY_TASK_ID
cd simulation_$SLURM_ARRAY_TASK_ID

# Activate conda environment
echo "Activating conda environment"
source ~/miniconda3/etc/profile.d/conda.sh

# Generate the FDS input file
echo "Generating FDS input file"
conda activate bandy-circles
python3 ../../../src/preprocess_from_csv.py $PWD $SLURM_ARRAY_TASK_ID

# Run FDS
echo "Running FDS"
fds "input_${SLURM_ARRAY_TASK_ID}.fds"