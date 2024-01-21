#!/bin/bash
#$ -q b.q
#$ -S /bin/bash
#$ -N TWAS
#$ -j y

# Set up environment
module load tabix/0.2.6
conda activate tigarenv
export PYTHONPATH=/home/rparrish/.conda/envs/tigarenv/lib/python3.5/site-packages/

# run TWAS
${TIGAR_dir}/TIGAR_TWAS.sh \
	--asso 2 \
	--chr ${chr} \
	--gene_anno ${gene_anno} \
	--LD ${LD} \
	--weight ${weight} \
	--Zscore ${Zscore} \
	--thread ${NSLOTS:-1} \
	--TIGAR_dir ${TIGAR_dir} \
	--sub_dir 0 \
	--out_prefix ${out_prefix} \
	--out_dir ${out_dir}

# UNLOAD MODULES
module unload tabix/0.2.6
conda deactivate
