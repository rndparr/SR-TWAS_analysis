#!/bin/bash
#$ -q b.q
#$ -N simtype1error
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 10
#$ -j y

############
# chmod 755 /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_srnaive_pred_type1error.sh

# qsub -v scenario=overlap_0.1_0.1_0.5_0.1_0.1_508 -v SR_model=_EN -v models="Naive SR" /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_srnaive_pred_type1error.sh

# qsub -v scenario=overlap_0.1_0.1_0.5_0.1_0.1_508 -v SR_model="" -v models="GTEx_EN ROSMAP" /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/submit_srnaive_pred_type1error.sh

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
# PRED
############
# load TIGAR python env
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/

# models=(Naive SR)
for model in ${models[@]}; do

	# from training
	# weight=${out_dir_train}/${model}_train_weight.txt.gz
	weight=${out_dir_train}/sims/${model}_weights${SR_model}
	pred_gene_anno=${out_dir_train}/${model}_info${SR_model}.txt

	# test data
	test_sampleID=${sim_dir}/sampleid/ROSMAP_test_800_sampleid.txt
	pred_genofile=${sim_dir}/genotype/ROSMAP_ABCA7_raw.dosage.gz

	# set names of output files
	out_pred_file=${model}_pred${SR_model}.txt
	log_file=${model}_pred_log${SR_model}.txt

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


