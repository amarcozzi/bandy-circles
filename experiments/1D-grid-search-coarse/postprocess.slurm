#!/bin/bash

#SBATCH --job-name=postprocess-1D-coarse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu(all)
#SBATCH --account=umontana_fire_modeling                          
#SBATCH --output=postprocess.out                        
#SBATCH --error=postprocess.err  

echo "Running postprocessing on job ${SLURM_ARRAY_TASK_ID}"


# Activate conda environment
echo "Activating conda environment"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bandy-circles

# Run postprocessing script
echo "Running postprocessing script"
python3 ../../src/postprocess.py .
