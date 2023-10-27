#!/bin/bash



# qsub -v weight1=GTEx_train_EN_weight -v weight2=ROSMAP_train_DPR_weight train_pred_srnaive.sh


#weight1=GTEx_train_EN_weight${suffix}
#weight2=ROSMAP_train_weight${suffix}

# train_pred_srnaive.sh \
# 	--weight1 GTEx_train_EN_weight \
# 	--weight2 ROSMAP_train_DPR_weight \
# 	--ncores ${ncores} \
#   --suffix ${suffix}

# weights=( )
# weights_names=( )

while [ $# -gt 0 ]; do
    if [[ $1 == *"--"* ]]; then
        v="${1/--/}"
        # if [[ $v == "weights" ]]; then
        #     while [[ $2 != *"--"* ]]; do
        #         if [[ $2 == "" ]]; then 
        #             break;
        #         fi
        #         weights+=("$2")
        #         shift
        #     done
        # elif [[ $v == "weights_names" ]]; then
        #     while [[ $2 != *"--"* ]]; do
        #         if [[ $2 == "" ]]; then 
        #             break;
        #         fi
        #         weights_names+=("$2")
        #         shift
        #     done        
        # else
            declare $v="$2"
        # fi
    fi
    shift
done


suffix=${suffix:-''}



################

# ncores=${NSLOTS:-1}

############
# load tabix
# module load tabix/0.2.6

# set directories
# sim_dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim
TIGAR_dir=${sim_dir}/scripts/TIGAR_SR_sim
expr_dir=${sim_dir}/expression/raw_dosage
out_dir_train=${sim_dir}/train/
out_dir_pred=${sim_dir}/pred/


############
# NAIVE/SR TRAIN
############
# load SR python env
source ${conda_path}
# conda init bash
conda activate ${sr_env}
# export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/
# set the PYTHONPATH
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.6/site-packages/:$PYTHONPATH


weight1_file=${out_dir_train}/sims/${weight1}
weight2_file=${out_dir_train}/sims/${weight2}

out_weight_file=SR_train_weight${suffix}
out_info_file=SR_train_info${suffix}.txt
out_naive_weight_file=Naive_train_weight${suffix}
out_naive_info_file=Naive_train_info${suffix}.txt
log_file=SR_Naive_train_log${suffix}.txt

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
--thread ${ncores} \
--train_sampleID ${valid_sampleID} \
--weights ${weight1_file} ${weight2_file} \
--out_dir ${out_dir_train}/sims/


mv ${out_dir_train}/sims/${out_naive_info_file} ${out_dir_train}/${out_naive_info_file}
mv ${out_dir_train}/sims/${out_info_file} ${out_dir_train}/${out_info_file}


# unload SR python env
conda deactivate


############
# PRED
############
# load TIGAR python env
conda activate ${tigar_env}
# export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH


models=(Naive SR)
for model in ${models[@]}; do

	# from training
	weight=${out_dir_train}/sims/${model}_train_weight${suffix}
	pred_gene_anno=${out_dir_train}/${model}_train_info${suffix}.txt


	# test data
	test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
	pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

	# set names of output files
	out_pred_file=${model}_pred${suffix}.txt
	log_file=${model}_pred_log${suffix}.txt


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

done



# unload TIGAR python env
conda deactivate

########################
# unload tabix
# module unload tabix/0.2.6



# ########################
# # setup pred results
# Rscript /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios2_pred_results_setup.R

# ########################
# # get pred results
# Rscript /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/scripts/scenarios2_pred_results.R


