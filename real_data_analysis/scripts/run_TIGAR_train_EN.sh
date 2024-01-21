#!/bin/bash
#$ -q b.q
#$ -S /bin/bash
#$ -N TIGAREN
#$ -j y

# Set up environment
module load tabix/0.2.6
conda activate tigarenv
export PYTHONPATH=/home/rparrish/.conda/envs/tigarenv/lib/python3.5/site-packages/

# run training
${TIGAR_dir}/TIGAR_Model_Train.sh \
--model elastic_net \
--gene_exp ${gene_exp} \
--train_sampleID ${train_sampleID} \
--genofile ${genofile} \
--chr ${chr} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--thread ${NSLOTS:-1} \
--out_prefix ${out_prefix} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${out_dir}

# unload modules
module unload tabix/0.2.6
conda deactivate





