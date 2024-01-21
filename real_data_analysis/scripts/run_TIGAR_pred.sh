#!/bin/bash
#$ -q b.q
#$ -S /bin/bash
#$ -N Pred
#$ -j y

# Set up environment
module load tabix/0.2.6
conda activate tigarenv
export PYTHONPATH=/home/rparrish/.conda/envs/tigarenv/lib/python3.5/site-packages/

# run prediction
${TIGAR_dir}/TIGAR_GReX_Pred.sh \
--chr ${chr} \
--weight ${weight} \
--gene_anno ${gene_anno} \
--test_sampleID ${test_sampleID} \
--genofile ${genofile} \
--format GT \
--genofile_type vcf \
--thread ${NSLOTS:-1} \
--sub_dir 0 \
--out_prefix ${out_prefix} \
--TIGAR_dir ${TIGAR_dir} \
--maf_diff 0 \
--out_dir ${out_dir}

# Unload modules
module unload tabix/0.2.6
conda deactivate
