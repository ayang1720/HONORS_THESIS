#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

python create_pickle_file_of_global_arrays.py 1 > creating_pickle_file_using_hours_log