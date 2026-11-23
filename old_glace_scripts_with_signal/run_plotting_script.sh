#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

python old_plot_glace_using_pickle_files.py > old_plot_glace_using_pickle_files_log