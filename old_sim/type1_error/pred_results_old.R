#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
# options(stringsAsFactors=FALSE, digits=16)
options(stringsAsFactors=FALSE)

# load libraries
library(reshape2)
library(doFuture)
library(foreach)

# set up parallel environment
registerDoFuture()
pln <- plan(multicore, workers=availableCores(methods='SGE'))
print(paste0('Starting job with ', availableCores(methods='SGE'), ' cores'))

######################################################
# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'
expr_dir <- 'raw_dosage/'

# sampleid
test_sampleid_path <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleid <- read.table(test_sampleid_path, header = FALSE)$V1




# output file paths
suffix <- '_overlap1'
all_out_pred_path <- paste0(dir, 'power/', 'all_pred_results', suffix, '.txt')
all_out_power_path <- paste0(dir, 'power/', 'all_power_results', suffix, '.txt')
# all_out_pred_path <- paste0(dir, 'power/', 'all_pred_results_test.txt')
# all_out_power_path <- paste0(dir, 'power/', 'all_power_results_test.txt')

# output columns
out_pred_cols <- c('job','i','TargetID','ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2','cohort','expr_R2','pheno_R2','pheno_pval')
out_power_cols <- c('job','ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2','cohort','power')

# # set up output files
# cat(paste(out_pred_cols, collapse='\t', sep=''), file=all_out_pred_path, sep='\n')
# cat(paste(out_power_cols, collapse='\t', sep=''), file=all_out_power_path, sep='\n')

######################################################
## FUNCTIONS

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
		c('ROSMAP','GTEx','Naive','SR'), pwr, 
		data=data, USE.NAMES=FALSE))))
}

# get power for job
job_pwr <- function(data) {
	res <- do.call(rbind.data.frame, by(data, data[,'Hp2'], hp2_pwr))
	rownames(res) <- NULL
	return(res)
}


oldjobstr <- function(job){
	x <- strsplit(job, '_')[[1]]
	if(length(x) == 3){
		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
		return(job)
	} else {
		return(paste(x[[1]], x[[5]], x[[4]], sep='_'))
	}
}

######################################################
# load expression variance
load(paste0(dir, 'power/data/expr_var', '_overlap1','.Rdata'))

# class(sim_data$GTEx_pred)
# sim_data$GTEx_pred <- as.numeric(sim_data$GTEx_pred)
# sim_data$ROSMAP_pred <- as.numeric(sim_data$ROSMAP_pred)
# sim_data$Naive_pred <- as.numeric(sim_data$Naive_pred)
# sim_data$SR_pred <- as.numeric(sim_data$SR_pred)


# sim_data[801,]

# get results for each job
for (job in jobs) {
	print(job)

	# get start time
	start_time <- Sys.time()

	# load job data
	load(paste0(dir, 'power/data/sim_data_', job, '.Rdata'))

	# get job info
	# job_vec <- as.numeric(strsplit(job, '_')[[1]])
	# causal_prop <- job_vec[1]
	# GTEx_He2 <- job_vec[2]
	# ROSMAP_He2 <- job_vec[3]
	job_info <- as.numeric(strsplit(job, '_')[[1]])
	ROSMAP_causal_prop  <- job_info[1]
	GTEx_causal_prop  <- job_info[2]
	GTEx_overlap  <- job_info[3]
	ROSMAP_He2  <- job_info[4]
	GTEx_He2  <- job_info[5]

	# should be 1:1000; included for testing with fewer
	i_pos <- length(strsplit(sim_data$TargetID[[1]], '_')[[1]])
	i_targets <- sort(unique(as.integer(
		do.call(rbind, strsplit(sim_data$TargetID, '_'))[,i_pos])))

	# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %:%
	foreach(i=i_targets, .combine=rbind) %dopar% {

		# get target id
		target_id <- paste0(job, '_', i)

		# subset sim_data for target
		target_data <- sim_data[sim_data$TargetID == target_id, -which(colnames(sim_data) == 'TargetID')]

		# variables for simulating phenotype
		##var_expr <- var(target_data$sim_expr)
		gamma <- sqrt(pheno_h2 / expr_var[[target_id]])

		error_term <- rnorm(nrow(target_data), mean=0, sd=sqrt(1 - pheno_h2))

		# Simulate phenotype from test expression
		target_data$sim_pheno <- gamma * target_data$sim_expr + error_term

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
				get_fit('ROSMAP', target_data),
				get_fit('GTEx', target_data), 
				get_fit('Naive', target_data),
				get_fit('SR', target_data)))

	} %seed% 1234567

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

	save(out_data, out_power, file=paste0(dir, 'power/data/out_power_pred_', job, '.Rdata'))
	rm(sim_data, out_data, out_power)

	# time elapsed
	end_time <- Sys.time()
	print(paste0('Computation time: ', end_time - start_time))
}
