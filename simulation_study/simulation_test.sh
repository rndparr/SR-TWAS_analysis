


# weights=( )
# weights_names=( )

# while [ $# -gt 0 ]; do
#     if [[ $1 == *"--"* ]]; then
#         v="${1/--/}"
#         if [[ $v == "weights" ]]; then
#             while [[ $2 != *"--"* ]]; do
#                 if [[ $2 == "" ]]; then 
#                     break;
#                 fi
#                 weights+=("$2")
#                 shift
#             done
#         elif [[ $v == "weights_names" ]]; then
#             while [[ $2 != *"--"* ]]; do
#                 if [[ $2 == "" ]]; then 
#                     break;
#                 fi
#                 weights_names+=("$2")
#                 shift
#             done        
#         else
#             declare $v="$2"
#         fi
#     fi
#     shift
# done

# suffix=_test



## load modules
# module load tabix/0.2.6

############
# setup
############

sim_dir=/Users/randyparr/github/SR-TWAS_analysis/simulation_study
scripts_dir=${sim_dir}/scripts
# cd ${sim_dir}

suffix=_batch1
ncores=1
jobs='0.01_0.01_0.5_0.5_0.5'
N_sim=10



conda_path=/opt/anaconda3/etc/profile.d/conda.sh
tigar_env=tigarenv
sr_env=py36

############
## make executable
############
chmod 755 ${scripts_dir}/TIGAR_SR_sim/DPR
chmod 755 ${scripts_dir}/train_pred_base.sh
chmod 755 ${scripts_dir}/train_pred_srnaive.sh
chmod 755 ${scripts_dir}/simulate_expression.R


############
## Get expression
############
# # genotype only
# paste0(sim_dir, 'genotype/', dataset, '_ABCA7_raw.dosage.gz')
RScript ${scripts_dir}/simulate_expression.R ${sim_dir}/ ${jobs} ${N_sim} ${suffix}


############
## train/pred base models
############
# ROSMAP_DPR
# ${scripts_dir}/train_pred_base.sh \
# 	--sim_dir ${sim_dir} \
# 	--conda_path ${conda_path} \
# 	--tigar_env ${tigar_env} \
# 	--dataset ROSMAP \
# 	--train_model DPR \
# 	--ncores ${ncores} \
# 	--suffix ${suffix}

${scripts_dir}/train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--conda_path ${conda_path} \
	--tigar_env ${tigar_env} \
	--dataset ROSMAP \
	--train_model EN \
	--ncores ${ncores} \
	--suffix ${suffix}


# GTEx_EN
${scripts_dir}/train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--conda_path ${conda_path} \
	--tigar_env ${tigar_env} \
	--dataset GTEx \
	--train_model EN \
	--ncores ${ncores} \
	--suffix ${suffix}

############
## train/pred sr-naive models
############
${scripts_dir}/train_pred_srnaive.sh \
	--sim_dir ${sim_dir} \
	--conda_path ${conda_path} \
	--sr_env ${sr_env} \
	--tigar_env ${tigar_env} \
	--weight1 GTEx_train_EN_weight \
	--weight2 ROSMAP_train_EN_weight \
	--ncores ${ncores} \
	--suffix ${suffix}

############
## pred results
############
combine_pred.sh

pred_results_setup.R

# power
pred_results.R 




############
## plot results
############
plot_setup.R

plot_both.R



module unload tabix/0.2.6


