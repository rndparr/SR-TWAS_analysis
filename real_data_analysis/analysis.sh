

TIGAR_dir=/home/rparrish/github/TIGAR
SR_TWAS_dir=/home/rparrish/github/SR-TWAS
wd=/home/rparrish/github/SR-TWAS_analysis/real_data_analysis/logs

# scripts_dir=.../SR-TWAS_analysis/real_data_analysis/analysis.sh

# Most jobs were done in-parallel; SGE was the resource manager for the machine used in analysis; if using a different resource manager, will need to modify submission scripts (${scripts_dir}/run_*.sh); -pe smp specifies number of cores

## TIGAR-TRAIN: DPR
qsub -wd ${wd} -pe smp 16 \
	-v chr=${chr} \
	-v TIGAR_dir=${TIGAR_dir} \
	-v train_sampleID=${train_sampleID} \
	-v gene_exp=${gene_exp} \
	-v genofile=${genofile} \
	-v out_prefix=CHR${chr}_TIGAR_train \
	-v out_dir=${out_dir} \
	${scripts_dir}/run_TIGAR_train_DPR.sh


## TIGAR-TRAIN: PrediXcan
qsub -wd ${wd} -pe smp 8 \
	-v chr=${chr} \
	-v TIGAR_dir=${TIGAR_dir} \
	-v train_sampleID=${train_sampleID} \
	-v gene_exp=${gene_exp} \
	-v genofile=${genofile} \
	-v out_prefix=CHR${chr}_PrediXcan_train \
	-v out_dir=${out_dir} \
	${scripts_dir}/run_TIGAR_train_EN.sh \

## SR-TWAS Train
qsub -wd ${wd} -pe smp 6 \
	-v chr=${chr} \
	-v SR_TWAS_dir=${SR_TWAS_dir} \
	-v valid_sampleID=${valid_sampleID} \
	-v gene_exp=${gene_exp} \
	-v genofile=${genofile} \
	-v weight1=${weight1} \
	-v weight2=${weight2} \
	-v weight3=${weight3} \
	-v weight4=${weight4} \
	-v weight5=${weight5} \
	-v weight6=${weight6} \
	-v weight_name1=${weight_name1} \
	-v weight_name2=${weight_name2} \
	-v weight_name3=${weight_name3} \
	-v weight_name4=${weight_name4} \
	-v weight_name5=${weight_name5} \
	-v weight_name6=${weight_name6} \
	-v out_dir=${out_dir} \
	-v out_prefix=CHR${chr}_SR-TWAS_train \
	${scripts_dir}/run_SR-TWAS_train.sh

## Naive Train
qsub -wd ${wd} -pe smp 8 \
	-v SR_TWAS_dir=${SR_TWAS_dir} \
	-v chr=${chr} \
	-v valid_sampleID=${valid_sampleID} \
	-v valid_sampleID=${valid_sampleID} \
	-v gene_exp=${gene_exp} \
	-v genofile=${genofile} \
	-v weight1=${weight1} \
	-v weight2=${weight2} \
	-v weight3=${weight3} \
	-v weight4=${weight4} \
	-v weight5=${weight5} \
	-v weight6=${weight6} \
	-v weight_name1=${weight_name1} \
	-v weight_name2=${weight_name2} \
	-v weight_name3=${weight_name3} \
	-v weight_name4=${weight_name4} \
	-v weight_name5=${weight_name5} \
	-v weight_name6=${weight_name6} \
	-v out_dir=${out_dir} \
	-v out_prefix=CHR${chr}_Naive_train \
	${scripts_dir}/run_Naive_train.sh

## Avg-valid+SR Train
qsub -wd ${wd} -pe smp 6 \
	-v SR_TWAS_dir=${SR_TWAS_dir} \
	-v chr=${chr} \
	-v gene_anno=${gene_anno} \
	-v weight1=${weight1} \
	-v weight2=${weight2} \
	-v weight3=${weight3} \
	-v weight_name1=${weight_name1} \
	-v weight_name2=${weight_name2} \
	-v weight_name3=${weight_name3} \
	-v out_dir=${out_dir} \
	-v out_prefix=CHR${chr}_Avg-valid-SR_train \
	${scripts_dir}/run_Avg-valid-SR_train.sh

## TIGAR pred
qsub -wd ${wd} -pe smp 4 \
	-v TIGAR_dir=${TIGAR_dir} \
	-v chr=${chr} \
	-v gene_anno=${gene_anno} \
	-v genofile=${genofile} \
	-v test_sampleID=${test_sampleID} \
	-v weight=${weight} \
	-v out_dir=${out_dir} \
	-v out_prefix=CHR${chr}_pred \
	${scripts_dir}/run_TIGAR_pred.sh

## TIGAR TWAS

# for example; using 
LD_str=reference_cov.txt.gz

qsub -wd ${wd} -pe smp 6 \
	-v TIGAR_dir=${TIGAR_dir} \
	-v chr=${chr} \
	-v gene_anno=${gene_anno} \
	-v LD_cohort=${LD_cohort} \
	-v weight=${weight} \
	-v Zscore=${Zscore} \
	-v out_prefix=CHR${chr}_TWAS \
	-v out_dir=${out_dir} \
	${scripts_dir}/run_TIGAR_TWAS.sh


