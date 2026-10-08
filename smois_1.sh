#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=4:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles1
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

export times="20210714000000 20210715000000 20210716000000 20210717000000 20210718000000 20210719000000 \
             20210720000000"

for time in ${times}
do
    echo "Starting ${time}"
    python plot_glace_ONLY.py ${time} > ${time}_log 2>&1
    echo "Finished ${time}"
done