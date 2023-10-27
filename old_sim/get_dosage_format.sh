#!/bin/bash
#$ -q b.q
#$ -N dosagefmt
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 3
#$ -j y

# qsub /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/get_dosage_format.sh

# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/get_dosage_format.sh

module load tabix/0.2.6
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/:$PYTHONPATH

python /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/get_dosage_format.py

conda deactivate
module unload tabix/0.2.6