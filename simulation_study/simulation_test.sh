#!/bin/bash


############
## SET UP
############
sim_dir='/home/rparrish/github/SR-TWAS_analysis/simulation_study'
cd ${sim_dir}/scripts

# set up directories
mkdir -p ${sim_dir}/{pred/logs,power/data,plot,train/{sims,logs}}

## number of cores for scripts that can run in parallel
ncores=1

## name for batch of files
suffix=''

## number of simulations to do (fixed at simulate expression step)
N_sim=2

## ggplot2 v3.3.0 library location; plotting step only
ggplot2_lib='/home/rparrish/lib/R/ggplot2_3.3.0'

## variables for generating jobs_list (see section for details)
prop_causal_snps=(0.05)
expr_hes=(0.2)
overlap=1
gtex_prop_causal_factor=1
gtex_expr_he_factor=1

## load modules/python
module load tabix/0.2.6
conda activate py36
export PYTHONPATH=/home/rparrish/.conda/envs/py36/lib/python3.6/site-packages/

############
## generate jobs list
############
## for plotting purposes, all jobs in a "batch" of jobs is assumed to have: 
	#* constant % of overlapped causal SNPs
	#* constant proportion (ROSMAP proportion of causal SNPs):(GTEx proportion of causal SNPs)
	#* constant proportion (ROSMAP heritability):(GTEx heritability) 
	#* all combinations of (ROSMAP prop_causal)*(ROSMAP expression heritability)
jobs_list=''
calc() { awk "BEGIN{print $*}"; }
for prop_causal_snp in ${prop_causal_snps[@]}; do
	gtex_prop_causal_snp=`calc ${prop_causal_snp}*${gtex_prop_causal_factor}`
	str_1=${prop_causal_snp}_${gtex_prop_causal_snp}_${overlap}
	for expr_he in ${expr_hes[@]}; do
		gtex_expr_he=`calc ${expr_he}*${gtex_expr_he_factor}`
		str_2=${expr_he}_${gtex_expr_he}
		jobs_list=`echo ${jobs_list},${str_1}_${str_2}`
	done
done
# remove extra comma at begginning of string
jobs_list="${jobs_list:1}"
echo $jobs_list

############
## simulate expression
############
# R packages: doFuture, foreach
Rscript ./simulate_expression.R ${sim_dir}/ ${jobs_list} ${suffix} ${ncores} ${N_sim}

############
## train/pred base models
############
# TIGAR-ROSMAP
./train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--train_model 'TIGAR' \
	--dataset 'ROSMAP' \
	--sampleid 'train_465' \
	--ncores ${ncores} \
	--suffix ${suffix}

# PrediXcan-GTEx
./train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--train_model 'PrediXcan' \
	--dataset 'GTEx' \
	--sampleid 'train_465' \
	--ncores ${ncores} \
	--suffix ${suffix}
	
# TIGAR-ROSMAP_valid
./train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--train_model 'TIGAR' \
	--dataset 'ROSMAP' \
	--sampleid 'valid_400' \
	--jobsuff '_valid' \
	--ncores ${ncores} \
	--suffix ${suffix}


############
## train/pred sr-naive models
############
# Naive-TIGAR-ROSMAP_PrediXcan-GTEx
# SR-TIGAR-ROSMAP_PrediXcan-GTEx
./train_pred_srnaive.sh \
	--sim_dir ${sim_dir} \
	--w1 'TIGAR-ROSMAP' \
	--w2 'PrediXcan-GTEx' \
	--ncores ${ncores} \
	--suffix ${suffix}


############
## train/pred avg-novalid models
############
# Avg-SRbasevalid
./train_pred_avg-novalid.sh \
	--sim_dir ${sim_dir} \
	--w1 'TIGAR-ROSMAP_valid' \
	--w2 'SR-TIGAR-ROSMAP_PrediXcan-GTEx' \
	--ncores ${ncores} \
	--suffix ${suffix}


############
## combine pred results
############
datasets=(TIGAR-ROSMAP TIGAR-ROSMAP_valid PrediXcan-GTEx Naive-TIGAR-ROSMAP_PrediXcan-GTEx SR-TIGAR-ROSMAP_PrediXcan-GTEx Avg-SRbasevalid)

# add header to output file
head -n1 ${sim_dir}/pred/${datasets[0]}_pred.txt | \
	awk 'BEGIN { FS = OFS = "\t" } {print "dataset\t"$0'} \
	> ${sim_dir}/pred/all_pred_results${suffix}.txt

# add data to output file
for dataset in ${datasets[@]}; do
	path=${sim_dir}/pred/${dataset}${suffix}_pred.txt
	tail -n+2 ${path} | \
		awk -v dataset=${dataset} 'BEGIN { FS = OFS = "\t" } {print dataset"\t"$0}' \
		>> ${sim_dir}/pred/all_pred_results${suffix}.txt
done


############
## results
############
# R packages: doFuture, foreach, reshape2
Rscript ./results_setup.R ${sim_dir}/ ${suffix} ${ncores}


############
## plot results
############
# reqires packages: ggsci, grid, gtable, gridExtra, paletteer, reshape2, ggplot2 v3.3.0
Rscript ./plots.R ${sim_dir}/ ${suffix} ${ggplot2_lib}


########################
module unload tabix/0.2.6
conda deactivate


