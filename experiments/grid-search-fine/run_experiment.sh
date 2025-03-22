echo "Running experiment: 1D-wind-speed-grid-search"

# Remove the output directory if it exists
rm -rf output

# Create the output directory
mkdir output

# Prep environment
conda activate bandy-circles
module load intel-oneapi-mpi

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
