echo "Running experiment: grid-search-fine"

# Remove the output directory if it exists
rm -rf output

# Create the output directory
mkdir output

# Prep environment
source $PROJECT_DIR/miniforge3/etc/profile.d/conda.sh
conda activate bandy-circles
module load intel-oneapi-mpi

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
