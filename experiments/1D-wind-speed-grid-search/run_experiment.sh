echo "Running experiment: 1D-wind-speed-grid-search"

# Check if the output directory exists
echo "Removing old output files"
rm out_*

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
