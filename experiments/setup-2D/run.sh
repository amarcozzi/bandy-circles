echo "Running 2D setup experiment"
rm *.out
sbatch submit.slurm
squeue -u $USER
