#!/bin/bash

# qsub -v dataset=GTEx -v train_model=EN -pe smp 10 train_pred_base.sh
# qsub -v dataset=ROSMAP -v train_model=DPR -pe smp 10 train_pred_base.sh

# dataset=GTEx
# dataset=ROSMAP

# train_model=DPR
# train_model=EN


while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done


suffix=${suffix:-''}

# ncores=${NSLOTS:-1}

############
# load modules
# module load tabix/0.2.6
source ${conda_path}
# conda init bash
conda activate ${tigar_env}
# export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH

# set directories
# sim_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim/
expr_dir=${sim_dir}/expression/raw_dosage


############
# TRAIN
############

# input files
Gene_Exp_train_file=${expr_dir}/${dataset}_expr${suffix}.txt
train_sample_path=${sim_dir}/sampleid/${dataset}_train_465_sampleid.txt
geno_file_path=${sim_dir}/genotype/${dataset}_ABCA7_raw.dosage.gz


# output files
out_dir_train=${sim_dir}/train/

# if [[ "$train_model"x == "DPR"x ]]; then
# 	model_str=''
# elif [[ "$train_model"x == "EN"x ]]; then
# 	model_str=_EN
# fi

out_weight_file=${dataset}_train_${train_model}_weight${suffix}
out_info_file=${dataset}_train_${train_model}_info${suffix}.txt
log_file=${dataset}_train_${train_model}_log${suffix}.txt

out_pred_file=${dataset}_${train_model}_pred${suffix}.txt
log_file=${dataset}_${train_model}_pred_log${suffix}.txt


if [[ "$train_model"x == "DPR"x ]]; then

	job_suf=${dataset}

	# run training
	python ${TIGAR_dir}/Train_py.py \
	--gene_exp ${Gene_Exp_train_file} \
	--train_sampleID ${train_sample_path} \
	--genofile ${geno_file_path} \
	--chr 19 \
	--genofile_type 'dosage' \
	--format 'DS' \
	--maf 0.01 \
	--hwe 0.00001 \
	--cvR2 1 \
	--dpr 1 \
	--ES fixed \
	--thread ${NSLOTS:-1} \
	--out_weight_file ${out_weight_file} \
	--out_info_file ${out_info_file} \
	--log_file ${log_file} \
	--job_suf ${job_suf} \
	--TIGAR_dir ${TIGAR_dir} \
	--sub_dir 0 \
	--out_dir ${out_dir_train}/sims



elif [[ "$train_model"x == "EN"x ]]; then

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

fi

mv ${out_dir_train}/sims/${out_info_file} ${out_dir_train}/${out_info_file}


############
# PRED
############
# from training
weight=${out_dir_train}/sims/${out_weight_file}
pred_gene_anno=${out_dir_train}/${out_info_file}

# test data
test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

out_dir_pred=${sim_dir}/pred/

python ${TIGAR_dir}/Pred_py.py \
	--chr 19 \
	--weight ${weight} \
	--gene_anno ${pred_gene_anno} \
	--test_sampleID ${test_sampleID} \
	--genofile ${pred_genofile} \
	--format 'DS' \
	--genofile_type 'dosage' \
	--window 1000000 \
	--thread ${ncores} \
	--maf_diff 0 \
	--missing_rate 0.2 \
	--out_pred_file ${out_pred_file} \
	--sub_dir 0 \
	--TIGAR_dir ${TIGAR_dir} \
	--log_file ${log_file} \
	--out_dir ${out_dir_pred}


# unload modules
# module unload tabix/0.2.6
conda deactivate



