


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
module load tabix/0.2.6

sim_dir=


## Get expression


# # genotype onli
# paste0(sim_dir, 'genotype/', dataset, '_ABCA7_raw.dosage.gz')


RScript simulate_expression.R ${sim_dir} ${jobs} ${N_sim} ${suffix}


############
## train/pred base models
############
# ROSMAP_DPR
train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--dataset ROSMAP \
	--train_model DPR \
	--ncores ${ncores} \
	--suffix ${suffix}

# GTEx_EN
train_pred_base.sh \
	--sim_dir ${sim_dir} \
	--dataset GTEx \
	--train_model EN \
	--ncores ${ncores} \
	--suffix ${suffix}

############
## train/pred sr-naive models
############
train_pred_srnaive.sh \
	--sim_dir ${sim_dir} \
	--weight1 GTEx_train_EN_weight \
	--weight2 ROSMAP_train_DPR_weight \
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


