#$ -q b.q
#$ -S /bin/bash
#$ -N NaiveTWAS
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/logs
#$ -j y
#$ -pe smp 16

###
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/run_Naive_train_multi-tissue.sh
# qsub /mnt/YangFSS/data2/rparrish/SR_TWAS/run_Naive_train_multi-tissue.sh
###

# arguments
SR_TWAS_dir=/home/rparrish/github/SR-TWAS
ROSMAP_dir=/mnt/YangFSS/data2/rparrish/ROSMAP_WGS

# Set up environment
module load tabix/0.2.6
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

for chr in `seq 1 22`; do

	genofile=/mnt/YangFSS/data/rparrish/GTEx_V8/GenotypeFiles/CHR${chr}_GTEx_WGS.vcf.gz

	# validation data
	valid_sampleID=/mnt/YangFSS/data/rparrish/GTEx_V8/SubjectIDFiles/Brain_Substantia_nigra_subjid.txt
	gene_exp=/mnt/YangFSS/data/rparrish/GTEx_V8/ExpressionFiles/Brain_Substantia_nigra_GTEx_Exp.txt


	out_prefix=CHR${chr}_Naive_train_GTEx

	# weight files
	weight1=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Anterior_cingulate_cortex_BA24/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name1=Brain_Anterior_cingulate_cortex_BA24

	weight2=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Caudate_basal_ganglia/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name2=Brain_Caudate_basal_ganglia

	weight3=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Cortex/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name3=Brain_Cortex

	weight4=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Nucleus_accumbens_basal_ganglia/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name4=Brain_Nucleus_accumbens_basal_ganglia

	weight5=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Brain_Putamen_basal_ganglia/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name5=Brain_Putamen_basal_ganglia

	weight6=/mnt/YangFSS/data/rparrish/GTEx_V8/TIGAR_train/Whole_Blood/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt.gz
	weight_name6=Whole_Blood


	# mkdir -p /mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/GTEx_6tissue/
	out_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/SR_TWAS_train/GTEx_6tissue/

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
		--thread ${NSLOTS:-1} \
		--train_sampleID ${valid_sampleID} \
		--out_prefix ${out_prefix} \
		--weights ${weight1} ${weight2} ${weight3} ${weight4} ${weight5} ${weight6} \
		--weights_names ${weight_name1} ${weight_name2} ${weight_name3} ${weight_name4} ${weight_name5} ${weight_name6}

done


# Unload modules
module unload tabix/0.2.6
conda deactivate
