#!/bin/bash

#SBATCH --job-name=grid-search
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-4%4                           
#SBATCH --nodes=4
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0:5:00
#SBATCH --output=simulations/%a.out                        
#SBATCH --error=simulations/%a.err
#SBATCH --requeue
#SBATCH --signal=SIGTERM@90

# Trap function to handle job termination
cleanup_and_requeue() {
    echo "Job ${SLURM_ARRAY_TASK_ID} received termination signal. Cleaning up and requeueing."
    scontrol requeue ${SLURM_JOB_ID}
    exit 0
}

# Set the trap
trap 'cleanup_and_requeue' SIGTERM

echo "Starting job ${SLURM_ARRAY_TASK_ID} with JOB_ID=${SLURM_JOB_ID} on ${SLURM_JOB_NODELIST}"

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
timeout --signal=SIGTERM 270s srun ~/fds/Build/fds_impi_intel_linux "input_${SLURM_ARRAY_TASK_ID}.fds"
FDS_EXIT_CODE=$?

if [ $FDS_EXIT_CODE -eq 124 ] || [ $FDS_EXIT_CODE -eq 137 ]; then
    echo "FDS timed out or was terminated. Triggering requeue."
    scontrol requeue $SLURM_JOB_ID
else
    echo "FDS completed with exit code $FDS_EXIT_CODE"
fi

echo "Job ${SLURM_ARRAY_TASK_ID} completed execution"