#!/bin/bash
#$ -q b.q
#$ -N bgzip
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -j y

# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/bgzip.sh
# qsub /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/bgzip.sh

module load tabix/0.2.6

dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train/sims

for i in `seq 1 1000000`; do

	bgzip ${dir}/GTEx_EN_weights_${i}.txt
	bgzip ${dir}/ROSMAP_weights_${i}.txt

	if [[ $(( $i % 2000 )) == 0 ]]; then
		echo ${i}
	fi

done

module unload tabix/0.2.6

