#!/bin/bash
#
#$ -q b.q
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/logs
#$ -pe smp 16
#$ -j y

#####
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/run_TIGAR_TWAS_all_chrms_LD_multi-tissue.sh



####

## SR-TWAS
# qsub -N PDTWAS -v LD_cohort=GTEx -v model=GTEx6tissues_Naive /mnt/YangFSS/data2/rparrish/SR_TWAS/run_TIGAR_TWAS_all_chrms_LD_multi-tissue.sh
# qsub -N ADTWAS2 -v LD_cohort=GTEx -v model=GTEx6tissues_Naive /mnt/YangFSS/data2/rparrish/SR_TWAS/run_TIGAR_TWAS_all_chrms_LD_multi-tissue.sh

# twas=PDTWAS

##########

### job info
# cohort=${model%%_*}
cohort=${model}
train_dat=${model#*_}
twas=${JOB_NAME:-ADTWAS2}

### directories
ROSMAP_dir=/mnt/YangFSS/data2/rparrish/ROSMAP_WGS
# TIGAR_dir=/home/rparrish/github/TIGAR
TIGAR_dir=/home/rparrish/TIGAR

## out_dir based on twas
if [[ "$twas"x == "ADTWAS"x ]]; then
	## AD TWAS
	out_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/AD_TWAS

elif [[ "$twas"x == "ADTWAS2"x ]]; then
	## AD TWAS2
	out_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/AD_TWAS2

elif [[ "$twas"x == "PDTWAS"x ]]; then
	## PD TWAS
	out_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/PD_TWAS
else
	echo "Error: Please specify twas."
	exit 1
fi

# out_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/scrap

# load modules
module load tabix/0.2.6
# conda activate myenv
# export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/

conda activate tigarenv
export PYTHONPATH=/home/rparrish/.conda/envs/tigarenv/lib/python3.5/site-packages/

echo -E "Starting TWAS for ${cohort} ${train_dat} with ${LD_cohort} LD data"

for chr in `seq 1 22`; do

# chr=21

	if [[ "$twas"x == "ADTWAS"x ]]; then
		## AD TWAS - iGAP GWAS
		Zscore=${out_dir}/IGAP_Zscore/CHR${chr}_Zscore.txt.gz	

	elif [[ "$twas"x == "ADTWAS2"x ]]; then
		## AD TWAS - GWAS; other dataset
		Zscore=/mnt/YangFSS/data2/thu/TIGAR/data/AD_sumstas/CHR${chr}_Zscore_GRCh38.txt.gz

	elif [[ "$twas"x == "PDTWAS"x ]]; then
		## PD TWAS - GWAS
		Zscore=${out_dir}/PD_GWAS/Zscore/CHR${chr}_nallsEtAl2019_excluding23andMe_allVariants_GRCh38_Zscore.txt.gz
	fi

	# LD file

	if [[ "$LD_cohort"x == "ROSMAP"x ]]; then
		LD=${ROSMAP_dir}/RefCovFiles/CHR${chr}_reference_cov.txt.gz
	elif [[ "$LD_cohort"x == "GTEx"x ]]; then
		LD=/mnt/YangFSS/data/rparrish/GTEx_V8/RefCovFiles/CHR${chr}_reference_cov.txt.gz
	else
		echo "Error: Please specify LD_cohort."
		exit 1
	fi

	## SR_TWAS infiles
	train_dat=${train_dat:-GTEx}

	in_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/GTEx_6tissue

	# weight=${in_dir}/CHR${chr}_SR_train_GTEx_eQTLweights.txt.gz
	# gene_anno=${in_dir}/CHR${chr}_SR_train_GTEx_GeneInfo.txt
	# out_prefix=${cohort}_GTEx_LD_${LD_cohort}_CHR${chr}_TWAS

	weight=${in_dir}/CHR${chr}_Naive_train_GTEx_eQTLweights.txt.gz
	gene_anno=${in_dir}/CHR${chr}_Naive_train_GTEx_GeneInfo.txt
	out_prefix=${cohort}_GTEx_LD_${LD_cohort}_CHR${chr}_TWAS

	echo -E "Starting CHR ${chr}."

	${TIGAR_dir}/TIGAR_TWAS.sh \
	--asso 2 \
	--gene_anno ${gene_anno} \
	--Zscore ${Zscore} \
	--weight ${weight} \
	--LD ${LD} \
	--chr ${chr} \
	--thread ${NSLOTS:-1} \
	--TIGAR_dir ${TIGAR_dir} \
	--sub_dir 0 \
	--out_prefix ${out_prefix} \
	--out_dir ${out_dir}

	echo "Done."

done

# UNLOAD MODULES
module unload tabix/0.2.6
conda deactivate


