echo "Running experiment: 1D-wind-speed-grid-search"

# Check if the output directory exists
echo "Removing old simulation files"
rm out_*
rm input_*

# Submit the job to slurm
sbatch submit.slurm
squeue -u $USER
