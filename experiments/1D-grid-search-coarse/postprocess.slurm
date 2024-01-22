#!/bin/bash

#SBATCH --job-name=postprocess-1D-coarse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --partition=cpu(all)
#SBATCH --account=umontana_fire_modeling                          
#SBATCH --output=output/%a.out                        
#SBATCH --error=output/%a.err  

echo "Running postprocessing on job ${SLURM_ARRAY_TASK_ID}"


# Activate conda environment
echo "Activating conda environment"
source ~/miniconda3/etc/profile.d/conda.sh

# Generate the FDS input file
echo "Generating FDS input file"
conda activate bandy-circles
python3 ../../src/postprocess.py $PWD
