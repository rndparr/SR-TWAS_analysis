#!/bin/bash


while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
      v="${1/--/}"
      declare $v="$2"
    fi
    shift
done


suffix=${suffix:-''}

################
# set directories
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim
expr_dir=${sim_dir}/expression
out_dir_train=${sim_dir}/train/
out_dir_pred=${sim_dir}/pred/


############
# NAIVE/SR TRAIN
############
# input sr training files
weight1=${out_dir_train}/sims/${w1}${suffix}_train_weight
weight2=${out_dir_train}/sims/${w2}${suffix}_train_weight

weight1_name=${w1}${suffix}
weight2_name=${w2}${suffix}

# sampleid and genotype/expression files
valid_sampleID=${sim_dir}/sampleid/ROSMAP_valid_400_sampleid.txt
genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz
gene_exp=${expr_dir}/ROSMAP_expr${suffix}.txt

# output files
out_weight_file=SR-${w1}_${w2}${suffix}_train_weight
out_info_file=SR-${w1}_${w2}${suffix}_train_info.txt
out_naive_weight_file=Naive-${w1}_${w2}${suffix}_train_weight
out_naive_info_file=Naive-${w1}_${w2}${suffix}_train_info.txt
log_file=SR_Naive-${w1}_${w2}${suffix}_train_log.txt

echo 'Training 'SR-${w1}_${w2}' & 'Naive-${w1}_${w2}' models.'
python ${TIGAR_dir}/SR-TWAS_Naive.py \
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
	--sub_dir 'sims' \
	--thread ${ncores} \
	--train_sampleID ${valid_sampleID} \
	--weights ${weight1} ${weight2} \
	--out_dir ${out_dir_train}/

############
# PRED
############
models=(Naive SR)
for model in ${models[@]}; do

	# from training
	weight=${out_dir_train}/sims/${model}-${w1}_${w2}${suffix}_train_weight
	pred_gene_anno=${out_dir_train}/${model}-${w1}_${w2}${suffix}_train_info.txt

	# test data
	test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
	pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

	# set names of output files
	out_pred_file=${model}-${w1}_${w2}${suffix}_pred.txt
	log_file=${model}-${w1}_${w2}${suffix}_pred_log.txt

	echo 'Running prediction for '${model}-${w1}_${w2}' model.'
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

done
