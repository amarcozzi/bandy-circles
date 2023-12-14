echo "Running experiment: 1D-wind-speed-grid-search"

# Check if the output directory exists
if [ -d "output" ]; then
  # Remove the output directory
  echo "Removing output directory"
  rm -rf output
fi
mkdir output

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
