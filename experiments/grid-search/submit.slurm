#!/bin/bash

#SBATCH --job-name=grid-search
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-59%5
#SBATCH --nodes=2
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --open-mode=append
#SBATCH --time=6:00:00
#SBATCH --output=simulations/%a.out                        
#SBATCH --error=simulations/%a.err
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@30

handle_sigterm() {
    echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Received SIGTERM"
    echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} received SIGTERM. Requeueing."
    scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
    exit 0
}

trap handle_sigterm SIGTERM

echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Starting job ${SLURM_ARRAY_TASK_ID} with JOB_ID=${SLURM_ARRAY_JOB_ID} on ${SLURM_JOB_NODELIST}"
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"

# Activate conda environment
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Activating conda environment"
source $PROJECT_DIR/miniforge3/etc/profile.d/conda.sh
conda activate bandy-circles

# Generate the FDS input file
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Generating FDS input file"
python3 $BANDY_DIR_SRC/preprocess_from_csv.py $PWD $SLURM_ARRAY_TASK_ID

# Activate MPI Module
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Activating MPI module"
module load intel-oneapi-mpi

# Run FDS
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Running FDS"
cd simulations/simulation_$SLURM_ARRAY_TASK_ID
srun ~/fds/Build/fds_impi_intel_linux "input_${SLURM_ARRAY_TASK_ID}.fds"

# Postprocess outputs
echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Postprocessing outputs"
python3 $BANDY_DIR_SRC/process_binary_data_to_zarr.py $PWD $SLURM_ARRAY_TASK_ID
rm out*

echo "$(date +%Y-%m-%d' '%H:%M:%S.%N): Job ${SLURM_ARRAY_TASK_ID} completed execution"