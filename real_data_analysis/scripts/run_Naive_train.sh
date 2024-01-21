#!/bin/bash
#$ -q b.q
#$ -S /bin/bash
#$ -N Naive
#$ -j y

# Set up environment
module load tabix/0.2.6
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

# run training
${SR_TWAS_dir}/Naive.sh \
	--chr ${chr} \
	--cvR2 1 \
	--format GT \
	--gene_exp ${gene_exp} \
	--genofile ${genofile} \
	--genofile_type vcf \
	--out_dir ${out_dir} \
	--SR_TWAS_dir ${SR_TWAS_dir} \
	--maf_diff 0 \
	--sub_dir 0 \
	--parallel ${NSLOTS:-1} \
	--train_sampleID ${valid_sampleID} \
	--out_prefix ${out_prefix} \
	--weights ${weight1} ${weight2} ${weight3} ${weight4} ${weight5} ${weight6} \
	--weights_names ${weight_name1} ${weight_name2} ${weight_name3} ${weight_name4} ${weight_name5} ${weight_name6}

# Unload modules
module unload tabix/0.2.6
conda deactivate
