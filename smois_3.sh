#!/bin/bash
#SBATCH --cluster=pitzer
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=pickles3
#SBATCH --account=PAS2635
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yang.6726@osu.edu

export times="20210727000000 20210728000000 20210729000000 20210730000000 20210731000000 \
             20210801000000"

for time in ${times}
do
    echo "Starting ${time}"
    python plot_glace_ONLY.py ${time} > ${time}_log 2>&1
    echo "Finished ${time}"
done