#!/bin/bash

#SBATCH --job-name=1d_parameter_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=48:00:00
#SBATCH --account=umontana_fire_modeling
#SBATCH --array=0-99%50                             
#SBATCH --output=output/%a.out                        
#SBATCH --error=output/%a.err  

echo "Starting job ${SLURM_ARRAY_TASK_ID}"
echo $SLURM_ARRAY_TASK_ID
echo $PWD

# Activate conda environment
echo "Activating conda environment"
conda ()
{ 
    \local cmd="${1-__missing__}";
    case "$cmd" in 
        activate | deactivate)
            __conda_activate "$@"
        ;;
        install | update | upgrade | remove | uninstall)
            __conda_exe "$@" || \return;
            __conda_reactivate
        ;;
        *)
            __conda_exe "$@"
        ;;
    esac
}

# Generate the FDS input file
echo "Generating FDS input file"
conda activate bandy-circles
python3 ../../src/preprocessing.py $PWD $SLURM_ARRAY_TASK_ID

# Run FDS
echo "Running FDS"
cd output
fds "input_${SLURM_ARRAY_TASK_ID}.fds"