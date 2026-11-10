#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

python /users/PAS3252/ayang1720/HONORS_THESIS/creates_the_3_by_2_plots/create_pickle_file_of_global_arrays_850hPa.py 1 > creating_pickle_file_using_hours_log