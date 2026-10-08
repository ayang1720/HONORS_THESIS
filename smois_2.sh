#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles2
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

export times="20210721000000 20210722000000 20210723000000 20210724000000 20210725000000 \
             20210726000000"

for time in ${times}
do
    echo "Starting ${time}"
    python plot_glace_ONLY.py ${time} > ${time}_log 2>&1
    echo "Finished ${time}"
done