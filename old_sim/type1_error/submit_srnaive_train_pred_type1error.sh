#!/bin/bash
#$ -q b.q
#$ -N simtype1error
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 16
#$ -j y

############
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/submit_srnaive_train_pred_type1error.sh

# qsub -v scenario=overlap_0.1_0.1_0.5_0.1_0.1_508 -v GTEx_model=_EN /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/submit_srnaive_train_pred_type1error.sh

############
# load tabix
module load tabix/0.2.6

# load conda command (testing only)
# source /sw/hgcc/Pkgs/Anaconda3/4.2.0/etc/profile.d/conda.sh 

# set directories
sim_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim

TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim
# expr_dir=${sim_dir}/expression/raw_dosage

out_dir_train=${sim_dir}/type1_error/${scenario}/train/
out_dir_pred=${sim_dir}/type1_error/${scenario}/pred/

############
# NAIVE/SR TRAIN
############
# load SR python env
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

# input sr training files
weight1=${out_dir_train}/sims/GTEx${GTEx_model}_weights
weight2=${out_dir_train}/sims/ROSMAP_weights


# sampleid and genotype/expression files
valid_sampleID=${sim_dir}/sampleid/ROSMAP_valid_400_sampleid.txt
genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz


gene_exp=${sim_dir}/type1_error/${scenario}/ROSMAP_expr.txt

annot=${sim_dir}/type1_error/TargetIDs.txt

# gene_exp=${sim_dir}/bleh.txt
out_weight_file=SR_weights${GTEx_model}
out_info_file=SR_info${GTEx_model}.txt
out_naive_weight_file=Naive_weights${GTEx_model}
out_naive_info_file=Naive_info${GTEx_model}.txt
log_file=SR_Naive_log${GTEx_model}.txt

python ${TIGAR_dir}/SR_TWAS_Naive_type1error.py \
--chr 19 \
--cvR2 1 \
--format 'DS' \
--annot ${annot} \
--gene_exp ${gene_exp} \
--genofile ${genofile} \
--genofile_type 'dosage' \
--out_weight_file ${out_weight_file} \
--out_info_file ${out_info_file} \
--out_naive_weight_file ${out_naive_weight_file} \
--out_naive_info_file ${out_naive_info_file} \
--log_file ${log_file} \
--SR_TWAS_dir ${TIGAR_dir} \
--maf_diff 0 \
--hwe 0.00001 \
--missing_rate 0.2 \
--naive 1 \
--sub_dir 0 \
--thread ${NSLOTS:-1} \
--train_sampleID ${valid_sampleID} \
--weights ${weight1} ${weight2} \
--out_dir ${out_dir_train}/sims/

# move info files up 1 directory
mv ${out_dir_train}/sims/${out_naive_info_file} ${out_dir_train}/${out_naive_info_file}
mv ${out_dir_train}/sims/${out_info_file} ${out_dir_train}/${out_info_file}

zcat ${out_dir_train}/sims/SR_weights_EN_2.txt.gz | head -n1 > ${out_dir_train}/sims/${out_weight_file}_header.txt
bgzip -f ${out_dir_train}/sims/${out_weight_file}_header.txt

zcat ${out_dir_train}/sims/Naive_weights_EN_2.txt.gz| head -n1 > ${out_dir_train}/sims/${out_naive_weight_file}_header.txt
bgzip -f ${out_dir_train}/sims/${out_naive_weight_file}_header.txt

# unload SR python env
conda deactivate


############
# PRED
############
# load TIGAR python env
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/

models=(Naive SR)
for model in ${models[@]}; do
	# from training
	# weight=${out_dir_train}/${model}_train_weight.txt.gz
	weight=${out_dir_train}/sims/${model}_weights${GTEx_model}
	pred_gene_anno=${out_dir_train}/${model}_info${GTEx_model}.txt

	# test data
	test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
	pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

	# set names of output files
	out_pred_file=${model}_pred${GTEx_model}.txt
	log_file=${model}_pred_log${GTEx_model}.txt

	python ${TIGAR_dir}/Pred_py.py \
	--chr 19 \
	--weight ${weight} \
	--gene_anno ${pred_gene_anno} \
	--test_sampleID ${test_sampleID} \
	--genofile ${pred_genofile} \
	--format 'DS' \
	--genofile_type 'dosage' \
	--window 1000000 \
	--thread ${NSLOTS:-1} \
	--maf_diff 0 \
	--missing_rate 0.2 \
	--out_pred_file ${out_pred_file} \
	--sub_dir 0 \
	--TIGAR_dir ${TIGAR_dir} \
	--log_file ${log_file} \
	--out_dir ${out_dir_pred}
done

# unload TIGAR python env
conda deactivate


########################
# unload tabix
module unload tabix/0.2.6


