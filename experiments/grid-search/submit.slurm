#!/bin/bash

#SBATCH --job-name=grid-search
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-4%4                           
#SBATCH --nodes=4
#SBATCH -n 88
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0:10:00
#SBATCH --output=simulations/%a.out                        
#SBATCH --error=simulations/%a.err
#SBATCH --requeue
#SBATCH --signal=SIGTERM@60

# Trap the SIGTERM signal and exit cleanly
trap 'clean_up' SIGTERM

clean_up() {
    echo "Job ${SLURM_ARRAY_TASK_ID} received SIGTERM signal. Cleaning up and exiting for requeue."
    rm -rf ../simulation_${SLURM_ARRAY_TASK_ID}
    exit 0  # Exit code 0 with --requeue flag will cause Slurm to requeue the job
}

echo "Starting job ${SLURM_ARRAY_TASK_ID}"

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
srun ~/fds/Build/fds_impi_intel_linux "input_${SLURM_ARRAY_TASK_ID}.fds"
