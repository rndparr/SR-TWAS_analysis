#!/bin/bash


while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
      v="${1/--/}"
      declare $v="$2"
    fi
    shift
done


suffix=${suffix:-''}

############
# set directories
TIGAR_dir=${TIGAR_dir:-${sim_dir}/scripts/TIGAR_SR_sim}
expr_dir=${expr_dir:-${sim_dir}/expression}
out_dir_train=${out_dir_train:-${sim_dir}/train/}
out_dir_pred=${out_dir_pred:-${sim_dir}/pred/}

############
# AVG-NOVALID TRAIN
############
# INPUT
weight1=${out_dir_train}/sims/${w1}${suffix}_train_weight
weight2=${out_dir_train}/sims/${w2}${suffix}_train_weight

weight1_name=${w1}${suffix}
weight2_name=${w2}${suffix}

gene_anno=${expr_dir}/ROSMAP_expr${suffix}.txt

## OUTPUT
out_weight_file=Avg-SRbasevalid${suffix}_train_weight
out_info_file=Avg-SRbasevalid${suffix}_train_info.txt

log_file=Avg-SRbasevalid${suffix}_train_log.txt

echo 'Training Avg-SRbasevalid model.'
python ${TIGAR_dir}/Naive_no_validation_data.py \
	--chr 19 \
	--gene_anno ${gene_anno} \
	--out_weight_file ${out_weight_file} \
	--out_info_file ${out_info_file} \
	--log_file ${log_file} \
	--SR_TWAS_dir ${TIGAR_dir} \
	--sub_dir 'sims' \
	--thread ${ncores} \
	--weights ${weight1} ${weight2} \
	--out_dir ${out_dir_train}


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
out_pred_file=Avg-SRbasevalid${suffix}_pred.txt
log_file=Avg-SRbasevalid${suffix}_pred_log.txt

echo 'Running prediction for Avg-SRbasevalid model.'

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



