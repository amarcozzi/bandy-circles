#!/bin/bash

#SBATCH --job-name=grid-search
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-59%5
#SBATCH --nodes=4
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --open-mode=append
#SBATCH --time=8:30:00
#SBATCH --output=simulations/%a.out                        
#SBATCH --error=simulations/%a.err
#SBATCH --requeue

echo "Starting job ${SLURM_ARRAY_TASK_ID} with JOB_ID=${SLURM_ARRAY_JOB_ID} on ${SLURM_JOB_NODELIST}"

# Activate conda environment
echo "Activating conda environment"
source $PROJECT_DIR/miniforge3/etc/profile.d/conda.sh
conda activate bandy-circles

# Generate the FDS input file
echo "Generating FDS input file"
python3 $BANDY_DIR_SRC/preprocess_from_csv.py $PWD $SLURM_ARRAY_TASK_ID

# Activate MPI Module
echo "Activating MPI module"
module load intel-oneapi-mpi

# Run FDS
echo "Running FDS"
cd simulations/simulation_$SLURM_ARRAY_TASK_ID
timeout 8h srun ~/fds/Build/fds_impi_intel_linux "input_${SLURM_ARRAY_TASK_ID}.fds"

if [[ $? -eq 124 || $? -eq 137 ]]; then
    echo "Job failed. Requeueing job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
else
    echo "Job ${SLURM_ARRAY_TASK_ID} completed execution"
fi