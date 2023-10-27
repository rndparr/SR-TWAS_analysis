#!/bin/bash
#$ -q b.q
#$ -N simsrnEN
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 4
#$ -j y

############
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_srnaive_train_pred_EN.sh

# qsub -v suffix=_overlap -hold_jid 448949 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_srnaive_train_pred_EN.sh

# qsub /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_srnaive_train_pred_EN.sh

# qsub -v suffix=_120621 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_srnaive_train_pred_EN.sh
# qsub -v suffix=_121721 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_srnaive_train_pred_EN.sh


############
# load tabix
module load tabix/0.2.6

# load conda command (testing only)
# source /sw/hgcc/Pkgs/Anaconda3/4.2.0/etc/profile.d/conda.sh 

# set directories
sim_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim
expr_dir=${sim_dir}/expression/raw_dosage
out_dir_train=${sim_dir}/train/
out_dir_pred=${sim_dir}/pred/

############
# NAIVE/SR TRAIN
############
# load SR python env
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

# input sr training files
# weight1=${out_dir_train}/GTEx_train_weight.txt.gz
# weight2=${out_dir_train}/ROSMAP_train_weight.txt.gz
weight1=${out_dir_train}/sims/GTEx_train_EN_weight${suffix}
weight2=${out_dir_train}/sims/ROSMAP_train_weight${suffix}

#{suffix}

# sampleid and genotype/expression files
valid_sampleID=${sim_dir}/sampleid/ROSMAP_valid_400_sampleid.txt
genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz
# gene_exp=${expr_dir}/ROSMAP_expr_EN${suffix}.txt
gene_exp=${expr_dir}/ROSMAP_expr${suffix}.txt

# gene_exp=${sim_dir}/bleh.txt
# set names of output files
# out_weight_file=SR_train_weight.txt
# out_info_file=SR_train_info.txt
# out_naive_weight_file=Naive_train_weight.txt
# out_naive_info_file=Naive_train_info.txt
# log_file=SR_Naive_train_log.txt
out_weight_file=SR_train_weight_EN${suffix}
out_info_file=SR_train_info_EN${suffix}.txt
out_naive_weight_file=Naive_train_weight_EN${suffix}
out_naive_info_file=Naive_train_info_EN${suffix}.txt
log_file=SR_Naive_train_log_EN${suffix}.txt

python ${TIGAR_dir}/SR_TWAS_Naive_py.py \
--chr 19 \
--cvR2 1 \
--format 'DS' \
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

mv ${out_dir_train}/sims/${out_naive_info_file} ${out_dir_train}/${out_naive_info_file}
mv ${out_dir_train}/sims/${out_info_file} ${out_dir_train}/${out_info_file}

zcat ${out_dir_train}/sims/SR_train_weight_0.001_0.1_0.1_1.txt.gz | head -n1 > ${out_dir_train}/sims/${out_weight_file}_header.txt
bgzip -f ${out_dir_train}/sims/${out_weight_file}_header.txt

zcat ${out_dir_train}/sims/Naive_train_weight_0.001_0.1_0.1_1.txt.gz | head -n1 > ${out_dir_train}/sims/${out_naive_weight_file}_header.txt
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
	weight=${out_dir_train}/sims/${model}_train_weight_EN${suffix}
	pred_gene_anno=${out_dir_train}/${model}_train_info_EN${suffix}.txt

	# test data
	test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
	pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

	# set names of output files
	out_pred_file=${model}_pred_EN${suffix}.txt
	log_file=${model}_pred_log_EN${suffix}.txt

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


# ########################
# # setup pred results
# Rscript /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios2_pred_results_setup.R

# ########################
# # get pred results
# Rscript /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios2_pred_results.R



