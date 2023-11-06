#!/bin/bash

while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done


jobsuff=${jobsuff:-''}
suffix=${suffix:-''}
sampleid=${sampleid:-train_465}


############
# set directories
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim/
expr_dir=${sim_dir}/expression
out_dir_train=${sim_dir}/train/

############
# TRAIN
############
# input files
Gene_Exp_train_file=${expr_dir}/${dataset}_expr${suffix}.txt
train_sample_path=${sim_dir}/sampleid/${dataset}_${sampleid}_sampleid.txt
geno_file_path=${sim_dir}/genotype/${dataset}_ABCA7_raw.dosage.gz


# output files
out_weight_file=${train_model}-${dataset}${jobsuff}${suffix}_train_weight
out_info_file=${train_model}-${dataset}${jobsuff}${suffix}_train_info.txt
log_file=${train_model}-${dataset}${jobsuff}${suffix}_train_log.txt
job_suf=${dataset}_${sampleid}${suffix}


echo 'Training '${train_model}-${dataset}${jobsuff}' model.'
if [[ "$train_model"x == "TIGAR"x ]]; then
	# run training
	python ${TIGAR_dir}/Train_TIGAR.py \
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
		--sub_dir 'sims' \
		--out_dir ${out_dir_train}

elif [[ "$train_model"x == "PrediXcan"x ]]; then
	python ${TIGAR_dir}/Train_PrediXcan.py \
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
		--sub_dir 'sims' \
		--out_dir ${out_dir_train}

fi


############
# PRED
############
# from training
weight=${out_dir_train}/sims/${out_weight_file}
pred_gene_anno=${out_dir_train}/${out_info_file}

# test data
test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

# set names of output files
out_dir_pred=${sim_dir}/pred/
out_pred_file=${train_model}-${dataset}${jobsuff}${suffix}_pred.txt
log_file=${train_model}-${dataset}${jobsuff}${suffix}_pred_log.txt


echo 'Running prediction for '${train_model}-${dataset}${jobsuff}' model.'
python ${TIGAR_dir}/Pred.py \
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
	--TIGAR_dir ${TIGAR_dir} \
	--log_file ${log_file} \
	--out_dir ${out_dir_pred}

