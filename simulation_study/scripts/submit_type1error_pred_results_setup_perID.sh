#!/bin/bash
#$ -q b.q
#$ -j y

Rscript ${sim_dir}/scripts/type1error_pred_results_setup_perID.R ${SGE_TASK_ID:-1} ${scenario_sim_i_suffix}

