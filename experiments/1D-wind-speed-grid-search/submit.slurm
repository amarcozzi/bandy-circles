#!/bin/bash

#SBATCH --job-name=1d_parameter_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=48:00:00
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-99%50                             
#SBATCH --output=output/%a.out                        
#SBATCH --error=output/%a.err  

# ulimit -s unlimited
# export OMP_NUM_THREADS=1

echo $SLURM_ARRAY_TASK_ID