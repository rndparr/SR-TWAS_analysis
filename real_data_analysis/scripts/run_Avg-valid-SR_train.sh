#!/bin/bash
#$ -q b.q
#$ -S /bin/bash
#$ -N AvgvalidSR
#$ -j y

# Set up environment
module load tabix/0.2.6
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

# run training
${SR_TWAS_dir}/Avg-valid_SR.sh \
	--chr ${chr} \
	--gene_anno ${gene_anno} \
	--out_dir ${out_dir} \
	--SR_TWAS_dir ${SR_TWAS_dir} \
	--sub_dir 0 \
	--parallel ${NSLOTS:-1} \
	--out_prefix ${out_prefix} \
	--weights ${weight1} ${weight2} ${weight3} \
	--weights_names ${weight_name1} ${weight_name2} ${weight_name3}

# Unload modules
module unload tabix/0.2.6
conda deactivate
