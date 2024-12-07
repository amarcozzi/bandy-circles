echo "Running experiment: 1D-wind-speed-grid-search"

# Remove the output directory if it exists
rm -rf output

# Create the output directory
mkdir output

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
