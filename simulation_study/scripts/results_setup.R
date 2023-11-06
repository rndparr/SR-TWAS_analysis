#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE)

# load libraries
library(doFuture)
library(foreach)
library(reshape2)

# load grguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==2) {
	sim_dir = as.character(args[[1]])
	suffix = ''
	ncores = as.integer(args[[2]])
} else {
	sim_dir = as.character(args[[1]])
	suffix = as.character(args[[2]])
	ncores = as.integer(args[[3]])
}

# set up parallel environment
registerDoFuture()
if (ncores > 1){
	pln <- plan(multicore, workers=ncores)
	print(paste0('Starting job with ', ncores, ' cores'))
} else {
	plan(sequential)
}

######################################################
# directories for simulation files
plot_dir <- paste0(sim_dir, 'plot/')

# sampleid
test_sampleid_path <- paste0(sim_dir, 'sampleid/', 'ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleid <- read.table(test_sampleid_path, header = FALSE)$V1

# parse arguments
datasets <- c('TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'PrediXcan-GTEx', 'Naive', 'SR', 'Avg')

datasets_filename_str <- setNames(c('TIGAR-ROSMAP', 'TIGAR-ROSMAP_valid', 'PrediXcan-GTEx', 'Naive-TIGAR-ROSMAP_PrediXcan-GTEx', 'SR-TIGAR-ROSMAP_PrediXcan-GTEx', 'Avg-SRbasevalid'), datasets)

# phenotype list
pheno_h2_list <- c(0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)

# output paths for results
all_out_pred_path <- paste0(sim_dir, 'power/', 'all_pred_results', suffix, '.txt')
all_out_power_path <- paste0(sim_dir, 'power/', 'all_power_results', suffix, '.txt')

# output columns for results
out_pred_cols <- c('job', 'i', 'TargetID', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2', 'cohort', 'expr_R2', 'pheno_R2', 'pheno_pval')
out_power_cols <- c('job', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2', 'cohort', 'power')

# set up output files; add header
cat(paste(out_pred_cols, collapse='\t', sep=''), file=all_out_pred_path, sep='\n')
cat(paste(out_power_cols, collapse='\t', sep=''), file=all_out_power_path, sep='\n')

########## 
# FUNCTIONS

# function to reshape dataframes
shapedat <- function(x, colname=NULL){
	if(is.null(colname)){ 
		colname <- deparse(substitute(x))
	}
	if ('job' %in% colnames(x)){ 
		idvars = c('TargetID', 'job')
	} else {
		idvars = 'TargetID'
	}
	return(melt(x, id.vars = idvars, variable.name = 'sample_id', value.name = colname))
}

# function to get fit and results safely in case of target failed in one of the datasets
get_fit <- function(x, dat){
	out_final <- tryCatch(
		expr = {
			expr_r2 <- cor(dat$sim_expr, dat[[paste0(x, '_pred')]])^2
			fit <- summary(lm(sim_pheno ~ get(paste0(x, '_pred')), data = dat))
			out <- data.frame('cohort'=x, 'expr_R2'=expr_r2, 
				'pheno_R2'=fit$r.squared * 100, 'pheno_pval'=fit$coefficients[2, 'Pr(>|t|)'])
			return(out)
		},
		error = function(e){
			# print(e)
			out <- data.frame('cohort'=x, 'expr_R2'=0, 'pheno_R2'=0, 'pheno_pval'=1)
			return(out)
		})
	return(out_final)
}

# function for calculating power for a cohort
pwr <- function(x, data){
	dat <- subset(data, cohort == x)
	power <- (length(which(dat[['pheno_pval']] < 2.5e-6)) / nrow(dat)) * 100
	## long
	out <- data.frame('cohort' = x, 'power' = power)
	## wide
	# out <- setNames(power, paste0(x, '_power'))
	return(out)
}

# get power values for a dataframe for a single Hp2
hp2_pwr <- function(data) {
	pheno <- unique(data[,'Hp2'])
	
	return(cbind(Hp2=pheno, t(sapply(
		datasets, pwr, 
		data=data, USE.NAMES=FALSE))))
}

# get power for job
job_pwr <- function(data) {
	res <- do.call(rbind.data.frame, by(data, data[,'Hp2'], hp2_pwr))
	rownames(res) <- NULL
	return(res)
}

######################################################
### SETUP
######################################################

## READ IN PREDICTION RESULTS
pred_path <- paste0(sim_dir, 'pred/all_pred_results', suffix, '.txt')
pred_cols <- colnames(read.table(pred_path, header=TRUE, check.names=FALSE, nrow=1))
pred_colClasses <- c('character', rep('NULL', 3), 'character', 'NULL', rep('numeric', length(pred_cols) - 6))
pred <- read.table(pred_path, header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE)
pred$job <- sub('_[^_]+$', '', pred$TargetID)

# read in prediction data, shape
get_pred_job <- function(dataset, job){
	ret_pred <- pred[which((pred$dataset == dataset) & (pred$job == job)), ]
	ret_pred <- ret_pred[, which(!(colnames(ret_pred)) %in% c('dataset'))]
	return(shapedat(ret_pred, paste0(dataset, '_pred')))
}

## READ IN SIMULATED EXPRESSION
# path to simulated expression
sim_expr_path <- paste0(sim_dir, 'expression/ROSMAP_expr', suffix, '.txt')
# read in number of columns, set column classes
sim_expr_ncols <- length(colnames(read.table(sim_expr_path, header = TRUE, check.names = FALSE, nrow=1)))
sim_expr_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', sim_expr_ncols - 5))
# read-in
sim_expr <- read.table(sim_expr_path, header=TRUE, check.names=FALSE, colClasses=sim_expr_colClasses, skipNul=TRUE)
sampleid <- colnames(sim_expr)[-1]
# center expression data
sim_expr[, sampleid] <- sim_expr[,sampleid] - rowMeans(sim_expr[,sampleid], na.rm=TRUE)
rm(list=c('sim_expr_ncols', 'pred_colClasses'))

# get jobs
sim_expr$job <- sub('_[^_]+$', '', sim_expr$TargetID)
jobs <- unique(sim_expr$job)
save(jobs, file=paste0(sim_dir, 'power/data/jobs', suffix, '.Rdata'))
job_df <- do.call(rbind, strsplit(sim_expr$job, '_'))

## get job info columns
sim_expr$ROSMAP_causal_prop <- job_df[,1]
sim_expr$GTEx_causal_prop <- job_df[,2]
sim_expr$GTEx_overlap <- job_df[,3]
sim_expr$ROSMAP_He2 <- job_df[,4]
sim_expr$GTEx_He2 <- job_df[,5]
rm(job_df)

# new variable for use in getting expression variance
sim_expr_allsamps <- sim_expr

# save unshaped dataset with only test sampleids
sim_expr <- sim_expr[, c('TargetID', 'job', test_sampleid)]
save(sim_expr, file=paste0(sim_dir, 'power/data/sim_expr_noshape_testsampleid', suffix, '.Rdata'))
rm(sim_expr)

## GET EXPRESSION VARIANCE
# reshape
sim_expr_allsamps <- shapedat(sim_expr_allsamps[,c('TargetID', 'job', test_sampleid)], 'sim_expr')
rm(list=c('test_sampleid_path', 'sim_expr_path', 'sim_expr_colClasses', 'sampleid')); # gc()

# expression variance (from ALL samples) for each simulation; used to get gamma for phenotype
expr_var <- by(sim_expr_allsamps$sim_expr, sim_expr_allsamps$TargetID, var, simplify=FALSE)

# save
save(expr_var, file=paste0(sim_dir, 'power/data/expr_var', suffix, '.Rdata'))
rm(list=c('sim_expr_allsamps', 'expr_var'))

## EXPORT DATA FOR EACH JOB
print('Getting sim_dat per job:')
for (job in jobs) {
	print(job)

	load(paste0(sim_dir, 'power/data/sim_expr_noshape_testsampleid', suffix,  '.Rdata'))

	sim_expr <- sim_expr[sim_expr$job == job, ]

	# reshape
	sim_expr <- shapedat(sim_expr)

	TIGAR_ROSMAP_pred <- get_pred_job('TIGAR-ROSMAP', job)
	TIGAR_ROSMAP_valid_pred <- get_pred_job('TIGAR-ROSMAP_valid', job)
	PrediXcan_GTEx_pred <- get_pred_job('PrediXcan-GTEx', job)
	Naive_pred <- get_pred_job('Naive-TIGAR-ROSMAP_PrediXcan-GTEx', job); colnames(Naive_pred)[4] <- 'Naive_pred'
	SR_pred <- get_pred_job('SR-TIGAR-ROSMAP_PrediXcan-GTEx', job); colnames(SR_pred)[4] <- 'SR_pred'
	Avg_pred <- get_pred_job('Avg-SRbasevalid', job); colnames(Avg_pred)[4] <- 'Avg_pred'

	sim_data <- merge(
		sim_expr,
		Reduce(function(...) merge(..., all=TRUE), list(TIGAR_ROSMAP_pred, TIGAR_ROSMAP_valid_pred, PrediXcan_GTEx_pred,	Naive_pred, SR_pred, Avg_pred)),
		all=TRUE)
	rm(list=c('TIGAR_ROSMAP_pred', 'TIGAR_ROSMAP_valid_pred', 'PrediXcan_GTEx_pred',	'Naive_pred', 'SR_pred', 'Avg_pred', 'sim_expr'))

	job_info <- strsplit(job, '_')[[1]]
	sim_data$ROSMAP_causal_prop <- job_info[1]
	sim_data$GTEx_causal_prop <- job_info[2]
	sim_data$GTEx_overlap <- job_info[3]
	sim_data$ROSMAP_He2 <- job_info[4]
	sim_data$GTEx_He2 <- job_info[5]

	sim_data_pred_cols <- paste0(datasets, '_pred')
	sim_data <- sim_data[, c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'job', 'TargetID', 'sample_id', 'sim_expr', sim_data_pred_cols)]

	save(sim_data, file=paste0(sim_dir, 'power/data/sim_data_', job, suffix, '.Rdata'))
	rm(sim_data)
}
rm(list=c('get_pred_job', 'shapedat', 'pred', 'sim_data_pred_cols'))

######################################################
### PRED RESULTS
######################################################
rm(jobs)
load(paste0(sim_dir, 'power/data/jobs', suffix, '.Rdata'))

# load expression variance
load(paste0(sim_dir, 'power/data/expr_var', suffix, '.Rdata'))

print('Prediction/Power results:')
for (job in jobs) {
	print(job)

	# get start time
	start_time <- Sys.time()

	# load job data
	load(paste0(sim_dir, 'power/data/sim_data_', job, suffix, '.Rdata'))

	# get job info
	job_info <- sim_data[1, c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2')]
	attach(job_info)

	# should be 1:1000 some jobs may have failed though, so
	i_pos <- length(strsplit(sim_data$TargetID[[1]], '_')[[1]])
	i_targets <- sort(unique(as.integer(
		do.call(rbind, strsplit(sim_data$TargetID, '_'))[, i_pos])))

	# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %:%
	foreach(i=i_targets, .combine=rbind) %dopar% {

		# get target id
		target_id <- paste0(job, '_', i)

		# subset sim_data for target
		target_data <- sim_data[sim_data$TargetID == target_id, -which(colnames(sim_data) == 'TargetID')]

		# variables for simulating phenotype
		gamma <- sqrt(pheno_h2 / expr_var[[target_id]])

		error_term <- rnorm(nrow(target_data), mean=0, sd=sqrt(1 - pheno_h2))

		# Simulate phenotype from test expression
		target_data$sim_pheno <- gamma * target_data$sim_expr + error_term


		head(target_data)
		# output for target
		out_i <- data.frame(
			'job'=job,
			'i'=i,
			'TargetID'=target_id, 
			'ROSMAP_causal_prop'=ROSMAP_causal_prop,
			'GTEx_causal_prop'=GTEx_causal_prop,
			'GTEx_overlap'=GTEx_overlap,
			'ROSMAP_He2'=ROSMAP_He2, 
			'GTEx_He2'=GTEx_He2,
			'Hp2'=pheno_h2,
			rbind(
				get_fit('TIGAR-ROSMAP', target_data),
				get_fit('TIGAR-ROSMAP_valid', target_data),
				get_fit('PrediXcan-GTEx', target_data), 
				get_fit('Naive', target_data),
				get_fit('SR', target_data),
				get_fit('Avg', target_data)))

	} %seed% 7654567

	# get power info
	out_power <- data.frame(
		'job'=job,
		'ROSMAP_causal_prop'=ROSMAP_causal_prop,
		'GTEx_causal_prop'=GTEx_causal_prop,
		'GTEx_overlap'=GTEx_overlap,
		'ROSMAP_He2'=ROSMAP_He2, 
		'GTEx_He2'=GTEx_He2,
		job_pwr(out_data))
	out_power <- apply(out_power, 2, as.character)

	# output pred results
	out_data <- out_data[, out_pred_cols]
	write.table(
		out_data,
		all_out_pred_path,
		quote=FALSE,
		append=TRUE,
		row.names=FALSE,
		col.names=FALSE,
		sep='\t')

	# output power results
	out_power <- out_power[, out_power_cols]
	write.table(
		out_power,
		all_out_power_path,
		quote=FALSE,
		row.names=FALSE,
		col.names=FALSE,
		append=TRUE,
		sep='\t')

	rm(sim_data, out_data, out_power)
	detach(job_info)

	# time elapsed
	end_time <- Sys.time()
	print(paste0('Computation time: ', end_time - start_time))
}

######################################################
### PLOT SETUP
######################################################
options(stringsAsFactors=FALSE, digits=5)
setwd(plot_dir)

# set up data frames
# sr_dat <- data.frame()
train_dat <- data.frame()

#############################
### TRAIN DAT
train_info_path <- function(x) { paste0(sim_dir, 'train/', datasets_filename_str[x], suffix, '_train_info.txt') }
# sr_rename <- function(col){ return(	sub(suffix, '', sub('Z_W', 'Z', sub('TIGAR-ROSMAP', 'W0', sub('PrediXcan-GTEx', 'W1', col)))) ) }

print('Setting up train_dat:')
for (dataset in datasets) {

	train_dat_path <- train_info_path(dataset)

	# get column names/classes for different train info filetypes
	if (startsWith(dataset, 'TIGAR')) {
		# CHROM, GeneStart, GeneEnd, TargetID, GeneName, sample_size, n_snp, n_effect_snp, CVR2, TrainPVALUE, TrainR2, CVR2_threshold
		col_classes <- rep('NULL', 12)
		col_classes[c(9,11)] <- rep('numeric', 2)
	} else if (startsWith(dataset, 'PrediXcan')) {
		# CHROM, GeneStart, GeneEnd, TargetID, GeneName, sample_size, n_snp, n_effect_snp, CVR2, TrainPVALUE, TrainR2, k-fold, alpha, Lambda, cvm, CVR2_threshold
		col_classes <- rep('NULL', 16)
		col_classes[c(9,11)] <- rep('numeric', 2)
	} else if(dataset == 'Naive') {
		# "CHROM, GeneStart, GeneEnd, TargetID, GeneName, sample_size, N_SNP, N_EFFECT_SNP, CVR2, R2, PVAL, W0_N_SNP, W0_CVR2, W0_R2, W0_PVAL, W1_N_SNP, W1_CVR2, W1_R2, W1_PVAL"
		col_classes <- rep('NULL', 19)
		col_classes[c(9:10, 13:14, 17:18)] <- rep('numeric', 6)
	} else if(dataset == 'SR') {
		# CHROM, GeneStart, GeneEnd, TargetID, GeneName, sample_size, N_SNP, N_EFFECT_SNP, CVR2, R2, PVAL, Z_TIGAR-ROSMAP, Z_PrediXcan-GTEx, TIGAR-ROSMAP_NSNP, TIGAR-ROSMAP_CVR2, TIGAR-ROSMAP_R2, TIGAR-ROSMAP_PVAL, PrediXcan-GTEx_NSNP, PrediXcan-GTEx_CVR2, PrediXcan-GTEx_R2, PrediXcan-GTEx_PVAL
		col_classes <- rep('NULL', 21)
		col_classes[c(9:10, 12:13, 15:16, 19:20)] <- rep('numeric', 8)
	} else if(dataset == 'Avg'){
		# CHROM, GeneStart, GeneEnd, TargetID, GeneName, N_SNP, N_EFFECT_SNP, TIGAR-ROSMAP_valid_NSNP, SR-TIGAR-ROSMAP_PrediXcan-GTEx_NSNP
		col_classes <- rep('NULL', 9)
		col_classes[6:9] <- rep('numeric', 4)
	}
	# TargetID
	col_classes[4] <- 'character'

	train_dat_part <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t', colClasses=col_classes, skipNul=TRUE)

	## Set up dataset pieces based on job type
	if(dataset == 'Avg'){
		train_dat_part$cohort <- dataset
		train_dat_part$CVR2 <- NA
		train_dat_part$TrainR2 <- NA		
		train_dat_part$ValidCVR2 <- NA
		train_dat_part$ValidR2 <- NA
		train_dat_part$Z0 <- NA
		train_dat_part$Z1 <- NA

	} else if(dataset == 'Naive' || dataset == 'SR') {

		# if(dataset == 'SR') {
		# 	sr_dat_part <- read.table(train_dat_path, header=TRUE, check.names=FALSE, sep='\t')[, c(1:5, 12,13)]
		# 	sr_dat_part$cohort <- dataset
		# 	sr_dat <- rbind(sr_dat_part, sr_dat)
		# } else 

		if (dataset == 'Naive') {
			train_dat_part$Z0 <- 0.5
			train_dat_part$Z1 <- 0.5
		}

		train_dat_part <- rbind(
			data.frame('TargetID'=train_dat_part$TargetID,  
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'=train_dat_part$CVR2, 
				'ValidR2'=train_dat_part$R2, 
				'cohort'=dataset,
				'Z0'=train_dat_part$Z0,
				'Z1'=train_dat_part$Z1),
			data.frame(
				'TargetID'=train_dat_part$TargetID, 
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'=train_dat_part$W0_CVR2, 
				'ValidR2'=train_dat_part$W0_R2, 
				'cohort'='PrediXcan-GTEx',
				'Z0'=NA,
				'Z1'=NA),
			data.frame(
				'TargetID'=train_dat_part$TargetID, 
				'CVR2'=NA, 'TrainR2'=NA,
				'ValidCVR2'= train_dat_part$W1_CVR2,
				'ValidR2'=train_dat_part$W1_R2, 
				'cohort'='TIGAR-ROSMAP',
				'Z0'=NA,
				'Z1'=NA))
	} else {
		colnames(train_dat_part) <- c('TargetID', 'CVR2', 'TrainR2')
		train_dat_part$cohort <- dataset
		train_dat_part$ValidCVR2 <- NA
		train_dat_part$ValidR2 <- NA
		train_dat_part$Z0 <- NA
		train_dat_part$Z1 <- NA
	}
	# merge into final file
	train_dat_part <- train_dat_part[, c('TargetID', 'CVR2', 'TrainR2', 'ValidCVR2', 'ValidR2', 'cohort', 'Z0', 'Z1')]
	train_dat <- rbind(train_dat_part, train_dat)
}

## get job info columns
train_dat$job <- sub('_[^_]+$', '', train_dat$TargetID)

job_df <- do.call(rbind, strsplit(train_dat$TargetID, '_'))
train_dat$ROSMAP_causal_prop <- job_df[,1]
train_dat$GTEx_causal_prop <- job_df[,2]
train_dat$GTEx_overlap <- job_df[,3]
train_dat$ROSMAP_He2 <- job_df[,4]
train_dat$GTEx_He2 <- job_df[,5]
train_dat$i <- job_df[,6]
rm(job_df)

train_dat$cohort <- factor(train_dat$cohort, levels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg'),
	labels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-Base+SR'))

#############################
### PRED DAT
print('Setting up pred_dat:')
pred_dat <- data.frame()
pred_dat_path <- paste0(sim_dir, 'power/all_pred_results', suffix, '.txt')
pred_dat <- read.table(pred_dat_path, header=TRUE, check.names=FALSE, sep='\t')
pred_dat <- pred_dat[, c('job', 'i', 'cohort', 'expr_R2', 'pheno_R2', 'pheno_pval', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2', 'TargetID')]
pred_dat$cohort <- factor(pred_dat$cohort, levels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg'),
	labels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-Base+SR'))

#############################
### POWER DAT
print('Setting up power_dat:')
power_dat <- data.frame()

power_dat_path <- paste0(sim_dir, 'power/all_power_results', suffix, '.txt')
power_dat <- read.table(power_dat_path, header=TRUE, check.names=FALSE, sep='\t')

# power_dat$oldjob <- power_dat$old_job
power_dat <- power_dat[, c('job', 'cohort', 'Hp2', 'power', 'ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2')]

power_dat$cohort <- factor(power_dat$cohort, levels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg'),
	labels=c('TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR-TWAS', 'Avg-Base+SR'))

#############################
# save plot data
save(train_dat, pred_dat, power_dat, file=paste0(sim_dir, 'plot/plot_data_scenario5', suffix, '.Rdata'))

