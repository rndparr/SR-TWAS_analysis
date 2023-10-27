#!/bin/bash
#
#
#$ -q b.q
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/logs
#$ -j y

#####
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_Rscript.sh


# qsub -N pred_results_setup -pe smp 10 -v script=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/pred_results_setup.R /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_Rscript.sh


# qsub -N pred_results -pe smp 10 -v script=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/pred_results.R /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_Rscript.sh

#####

Rscript ${script}

