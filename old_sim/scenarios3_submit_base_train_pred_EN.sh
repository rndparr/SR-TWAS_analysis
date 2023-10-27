#!/bin/bash
#$ -q b.q
#$ -N simbaseEN
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 8
#$ -j y

############
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh

# qsub -v dataset=GTEx -v suffix=_overlap /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh


# qsub -v dataset=GTEx -v suffix=_120621 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh
# qsub -v dataset=GTEx -v suffix=_121721 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh

# qsub -v dataset=GTEx -v suffix=_123121 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh

# qsub -v dataset=ROSMAP /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh

# qsub -v dataset=GTEx -v suffix=_zeros /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh
# qsub -v dataset=ROSMAP -v suffix=_zeros /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios3_submit_base_train_pred_EN.sh

############
# load modules
module load tabix/0.2.6
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/

# set directories
sim_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim/
expr_dir=${sim_dir}/expression/raw_dosage

############
# TRAIN
############
# simulation data
# Gene_Exp_train_file=${expr_dir}/${dataset}_expr.txt
# Gene_Exp_train_file=${expr_dir}/${dataset}_expr_EN${suffix}.txt
Gene_Exp_train_file=${expr_dir}/${dataset}_expr${suffix}.txt

train_sample_path=${sim_dir}/sampleid/${dataset}_train_465_sampleid.txt
geno_file_path=${sim_dir}/genotype/${dataset}_ABCA7_raw.dosage.gz

out_dir_train=${sim_dir}/train/
# out_weight_file=${dataset}_train_weight.txt
out_weight_file=${dataset}_train_EN_weight${suffix}
out_info_file=${dataset}_train_EN_info${suffix}.txt
log_file=${dataset}_train_EN_log${suffix}.txt

# # scrap
# out_dir_train=${sim_dir}/train/scrap/
# Gene_Exp_train_file=${sim_dir}/expression/${dataset}_expr_scrap.txt

# run training
python ${TIGAR_dir}/Train_EN_2.py \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_path} \
--genofile ${geno_file_path} \
--chr 19 \
--genofile_type 'dosage' \
--format 'DS' \
--maf 0.01 \
--hwe 0.00001 \
--cvR2 1 \
--alpha 0.5 \
--use_alpha 1 \
--thread ${NSLOTS:-1} \
--out_weight_file ${out_weight_file} \
--out_info_file ${out_info_file} \
--log_file ${log_file} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${out_dir_train}/sims

mv ${out_dir_train}/sims/${out_info_file} ${out_dir_train}/${out_info_file}

zcat ${out_dir_train}/sims/${dataset}_train_EN_weight_0.001_0.1_0.1_1.txt.gz | head -n1 > ${out_dir_train}/sims/${out_weight_file}_header.txt
bgzip -f ${out_dir_train}/sims/${out_weight_file}_header.txt


############
# PRED
############
# from training
# weight=${out_dir_train}/${out_weight_file}.gz
weight=${out_dir_train}/sims/${out_weight_file}
pred_gene_anno=${out_dir_train}/${out_info_file}

# test data
test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

# set names of output files
out_dir_pred=${sim_dir}/pred/
out_pred_file=${dataset}_EN_pred${suffix}.txt
log_file=${dataset}_EN_pred_log${suffix}.txt

# # scrap
# out_dir_pred=${sim_dir}/pred/scrap
# weight=${sim_dir}/train/scrap/${out_weight_file}.gz
# pred_gene_anno=${sim_dir}/train/scrap/${out_info_file}

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

# unload modules
module unload tabix/0.2.6
conda deactivate
