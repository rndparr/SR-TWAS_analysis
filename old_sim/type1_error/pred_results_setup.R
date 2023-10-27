#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)

# load libraries
library(reshape2)
library(doFuture)
library(foreach)
library(iterators)

# set up parallel environment
registerDoFuture()
pln <- plan(multicore, workers=availableCores(methods='SGE'))
print(paste0('Starting job with ', availableCores(methods='SGE'), ' cores'))

# pheno_h2_list <- c(0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)

######################################################
# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'

# sampleid
test_sampleids_path <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleids <- read.table(test_sampleids_path, header = FALSE)$V1


##########
## SIM_EXPR

# path to simulated expression
sim_expr_path <- paste0(dir, 'ROSMAP_expr.txt')

# read in number of columns, set column classes
sim_expr_ncols <- length(colnames(read.table(sim_expr_path, header = TRUE, check.names = FALSE, nrow=1)))
sim_expr_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', sim_expr_ncols - 5))

# read-in
sim_expr <- read.table(sim_expr_path, header=TRUE, check.names=FALSE, colClasses=sim_expr_colClasses, skipNul=TRUE)
sampleid <- colnames(sim_expr)[-1]

# center expression data
sim_expr[,sampleid] <- sim_expr[,sampleid] - rowMeans(sim_expr[,sampleid], na.rm=TRUE)

# reduce to test sampleids and reshape
sim_expr <- sim_expr[, test_sampleids]
sim_expr <- melt(sim_expr, measure.vars=test_sampleids, variable.name='sample_id', value.name='sim_expr')

# # get variance ## NEEDS TO BE FOR ALL SAMPLES NOT JUST TEST_SAMPLEIDS
# expr_var <- var(sim_expr$sim_expr)



######################################################
pred_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/split_files/'

out_pred_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/target_files/'
datasets <- c('ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')
out_cols <- c('TargetID', 'sample_id', 'sim_expr', datasets)

# print header to its own file
cat(paste(out_cols[-1], collapse='\t', sep=''), file=paste0(out_pred_dir, 'header.txt'), sep='\n')

# read in pred data for a dataset,test_sampleid
read_pred <- function(dataset, pred_dir, test_sampleid){
	pred_path <- paste0(pred_dir, dataset, '_', test_sampleid, '.txt')
	pred_dataset <- read.table(pred_path, header=FALSE, colClasses=c('character', 'numeric'), fill=TRUE)
	colnames(pred_dataset) <- c('TargetID', dataset)
	return(pred_dataset)
}

# output to per-target files
output_pred <- function(row) {
	cat(paste(row[out_cols[-1]], collapse='\t', sep=''), file=paste0(out_pred_dir, row[['TargetID']], '.txt'), sep='\n', append=TRUE)
}

targets <- c()

ith_test_sampleid <- setNames(1:length(test_sampleids), test_sampleids)
# SM-CJFLS

# test_sampleids[402:length(test_sampleids)]

# combine dataset,test_sampleid data per test_sampleid, output by target
# for (test_sampleid in test_sampleids) {
for (test_sampleid in test_sampleids[402:length(test_sampleids)]){

	print(ith_test_sampleid[[test_sampleid]])

	sim_data <- Reduce(function(...) merge(..., all=TRUE), 
		sapply(datasets, read_pred, pred_dir=pred_dir, test_sampleid=test_sampleid, simplify=FALSE))

	sim_data$TargetID <- as.integer(sim_data$TargetID)
	sim_data <- sim_data[order(sim_data$TargetID), ]
	sim_data$sample_id <- test_sampleid
	sim_data$sim_expr <- sim_expr[sim_expr$sample_id == test_sampleid, 'sim_expr']
	sim_data <- sim_data[, out_cols]

	targets <- sort(unique(c(sim_data$TargetID, targets)))

	# output to per-target files
	# apply(sim_data, 1, output_pred)

	foreach(sim_dat=iter(sim_data, by='row'), .combine=c) %dopar% {
		output_pred(sim_dat)
	}


}

# output list of targets
cat(paste(targets, sep='\n'), file=paste0(out_pred_dir, 'targets.txt') )


# ##################################################################
# # function to get fit and results safely in case of target failed in one of the datasets
# get_fit <- function(x, dat, sim_expr_scalar){
# 	out_final <- tryCatch(
# 		expr = {
# 			expr_r2 <- cor(dat$sim_expr, dat[[x]])^2
# 			fit <- summary(lm(sim_pheno ~ get(x), data = dat))
# 			out <- data.frame('cohort'=x, 'expr_R2'=expr_r2, 
# 				'pheno_R2'=fit$r.squared * 100, 'pheno_pval'=fit$coefficients[2, 'Pr(>|t|)'])
# 			return(out)
# 		},
# 		error = function(e){
# 			print(e)
# 			out <- data.frame('cohort'=x, 'expr_R2'=0, 'pheno_R2'=0, 'pheno_pval'=1)
# 			return(out)
# 		})
# 	return(out_final)
# }

# for (target in targets){

# 	sim_data <- read.table(paste0(out_pred_dir, target, '.txt'), header=FALSE, colClasses=c('character', rep('numeric', 5)), fill=TRUE)
# 	colnames(sim_data) <- out_cols[-1]

# 	# for each phenotype
# 	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %dopar% {

# 		# variables for simulating phenotype
# 		gamma <- sqrt(pheno_h2 / expr_var)

# 		# error term
# 		error_term <- rnorm(nrow(sim_data), mean=0, sd=sqrt(1 - pheno_h2))

# 		# Simulate phenotype from test expression
# 		sim_data$sim_pheno <- gamma * sim_data$sim_expr + error_term

# 	} %seed% 1234567
# }


# get_fit('ROSMAP', sim_data)


# ############


# # function for calculating power for a cohort
# pwr <- function(x, data){
# 	dat <- subset(data, cohort == x)
# 	power <- (length(which(dat[['pheno_pval']] < 2.5e-6)) / nrow(dat)) * 100
# 	## long
# 	out <- data.frame('cohort' = x, 'power' = power)
# 	## wide
# 	# out <- setNames(power, paste0(x, '_power'))
# 	return(out)
# }

# # get power values for a dataframe for a single Hp2
# hp2_pwr <- function(data) {
# 	pheno <- unique(data[,'Hp2'])
# 	return(cbind(Hp2=pheno, t(sapply(
# 		c('ROSMAP','GTEx','Naive','SR'), pwr, 
# 		data=data, USE.NAMES=FALSE))))
# }

# # get power for job
# job_pwr <- function(data) {
# 	res <- do.call(rbind.data.frame, by(data, data[,'Hp2'], hp2_pwr))
# 	rownames(res) <- NULL
# 	return(res)
# }



# ##################################################################
# ##################################################################

# # 	# i_targets <- sort(unique(as.integer(sim_data$TargetID)))

# # 	# # for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
# # 	# pheno_h2 <- pheno_h2_list[[1]]
# # 	# out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %dopar% {
# # 	# # %:%
# # 	# # foreach(i=i_targets, .combine=rbind) %dopar% {

# # 	# 	test_sampleid_sim_expr <- sim_expr[sim_expr$sample_id == test_sampleid, 'sim_expr']

# # 	# 	# variables for simulating phenotype
# # 	# 	gamma <- sqrt(pheno_h2 / expr_var)

# # 	# 	error_term <- rnorm(nrow(sim_data), mean=0, sd=sqrt(1 - pheno_h2))

# # 	# 	# Simulate phenotype from test expression
# # 	# 	sim_data$sim_pheno <- gamma * test_sampleid_sim_expr + error_term



# # 	# } %seed% 1234567

# # # }


# # # function to get fit and results safely in case of target failed in one of the datasets
# # get_fit <- function(x, dat, sim_expr_scalar){
# # 	out_final <- tryCatch(
# # 		expr = {
# # 			sim_expr <- rep(sim_expr_scalar, nrow(dat))
# # 			expr_r2 <- cor(sim_expr, dat[[x]])^2
# # 			fit <- summary(lm(sim_pheno ~ get(x), data = dat))
# # 			out <- data.frame('cohort'=x, 'expr_R2'=expr_r2, 
# # 				'pheno_R2'=fit$r.squared * 100, 'pheno_pval'=fit$coefficients[2, 'Pr(>|t|)'])
# # 			return(out)
# # 		},
# # 		error = function(e){
# # 			print(e)
# # 			out <- data.frame('cohort'=x, 'expr_R2'=0, 'pheno_R2'=0, 'pheno_pval'=1)
# # 			return(out)
# # 		})
# # 	return(out_final)
# # }

# # get_fit('ROSMAP', sim_data, test_sampleid_sim_expr)

# # datasets

# # nrow(sim_data)
# # rep(test_sampleid_sim_expr, 220686)

# # head(rep(test_sampleid_sim_expr, 220686))
# # head(sim_data$ROSMAP)

# # ########## 
# # # SETUP FUNCTIONS

# # # function to reshape dataframes
# # shapedat <- function(x, colname=NULL){
# # 	if(is.null(colname)){ 
# # 		colname <- deparse(substitute(x))
# # 	}
# # 	if ('job' %in% colnames(x)){ 
# # 		idvars = c('TargetID', 'job')
# # 	} else {
# # 		idvars = 'TargetID'
# # 	}
# # 	return(melt(x, id.vars = idvars, variable.name = 'sample_id', value.name = colname))
# # }

# # # memory usage of all objects 
# # memuse <- function() {
# # 	for (object in ls(envir=sys.frame())){
# # 		print(paste0(object, ':    ', format(object.size(eval(parse(text=object))), standard='SI', unit='auto')))
# # 	}
# # }


# # ############################################################################################################


# # #!/usr/bin/env Rscript

# # # set options
# # Sys.setlocale('LC_ALL', 'C')
# # .libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
# # # options(stringsAsFactors=FALSE, digits=16)
# # options(stringsAsFactors=FALSE)

# # # load libraries
# # library(reshape2)
# # library(doFuture)
# # library(foreach)

# # # set up parallel environment
# # registerDoFuture()
# # pln <- plan(multicore, workers=availableCores(methods='SGE'))
# # print(paste0('Starting job with ', availableCores(methods='SGE'), ' cores'))

# # ######################################################
# # # directory for simulation files
# # dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'
# # expr_dir <- 'raw_dosage/'

# # # sampleid
# # test_sampleids_path <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_test_800_sampleid.txt')

# # # read in test sampleid data
# # test_sampleids <- read.table(test_sampleids_path, header = FALSE)$V1




# # # output file paths
# # suffix <- '_overlap1'
# # all_out_pred_path <- paste0(dir, 'power/', 'all_pred_results', suffix, '.txt')
# # all_out_power_path <- paste0(dir, 'power/', 'all_power_results', suffix, '.txt')
# # # all_out_pred_path <- paste0(dir, 'power/', 'all_pred_results_test.txt')
# # # all_out_power_path <- paste0(dir, 'power/', 'all_power_results_test.txt')

# # # output columns
# # out_pred_cols <- c('job','i','TargetID','ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2','cohort','expr_R2','pheno_R2','pheno_pval')
# # out_power_cols <- c('job','ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'Hp2','cohort','power')

# # # # set up output files
# # # cat(paste(out_pred_cols, collapse='\t', sep=''), file=all_out_pred_path, sep='\n')
# # # cat(paste(out_power_cols, collapse='\t', sep=''), file=all_out_power_path, sep='\n')

# # ######################################################
# # ## FUNCTIONS

# # # function to get fit and results safely in case of target failed in one of the datasets
# # get_fit <- function(x, dat){
# # 	out_final <- tryCatch(
# # 		expr = {
# # 			expr_r2 <- cor(dat$sim_expr, dat[[paste0(x, '_pred')]])^2
# # 			fit <- summary(lm(sim_pheno ~ get(paste0(x, '_pred')), data = dat))
# # 			out <- data.frame('cohort'=x, 'expr_R2'=expr_r2, 
# # 				'pheno_R2'=fit$r.squared * 100, 'pheno_pval'=fit$coefficients[2, 'Pr(>|t|)'])
# # 			return(out)
# # 		},
# # 		error = function(e){
# # 			# print(e)
# # 			out <- data.frame('cohort'=x, 'expr_R2'=0, 'pheno_R2'=0, 'pheno_pval'=1)
# # 			return(out)
# # 		})
# # 	return(out_final)
# # }

# # # function for calculating power for a cohort
# # pwr <- function(x, data){
# # 	dat <- subset(data, cohort == x)
# # 	power <- (length(which(dat[['pheno_pval']] < 2.5e-6)) / nrow(dat)) * 100
# # 	## long
# # 	out <- data.frame('cohort' = x, 'power' = power)
# # 	## wide
# # 	# out <- setNames(power, paste0(x, '_power'))
# # 	return(out)
# # }

# # # get power values for a dataframe for a single Hp2
# # hp2_pwr <- function(data) {
# # 	pheno <- unique(data[,'Hp2'])
# # 	return(cbind(Hp2=pheno, t(sapply(
# # 		c('ROSMAP','GTEx','Naive','SR'), pwr, 
# # 		data=data, USE.NAMES=FALSE))))
# # }

# # # get power for job
# # job_pwr <- function(data) {
# # 	res <- do.call(rbind.data.frame, by(data, data[,'Hp2'], hp2_pwr))
# # 	rownames(res) <- NULL
# # 	return(res)
# # }


# # oldjobstr <- function(job){
# # 	x <- strsplit(job, '_')[[1]]
# # 	if(length(x) == 3){
# # 		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
# # 		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
# # 		return(job)
# # 	} else {
# # 		return(paste(x[[1]], x[[5]], x[[4]], sep='_'))
# # 	}
# # }

# # ######################################################
# # # load expression variance
# # load(paste0(dir, 'power/data/expr_var', '_overlap1','.Rdata'))

# # # class(sim_data$GTEx_pred)
# # # sim_data$GTEx_pred <- as.numeric(sim_data$GTEx_pred)
# # # sim_data$ROSMAP_pred <- as.numeric(sim_data$ROSMAP_pred)
# # # sim_data$Naive_pred <- as.numeric(sim_data$Naive_pred)
# # # sim_data$SR_pred <- as.numeric(sim_data$SR_pred)


# # # sim_data[801,]

# # # get results for each job
# # for (job in jobs) {
# # 	print(job)

# # 	# get start time
# # 	start_time <- Sys.time()

# # 	# load job data
# # 	load(paste0(dir, 'power/data/sim_data_', job, '.Rdata'))

# # 	# get job info
# # 	# job_vec <- as.numeric(strsplit(job, '_')[[1]])
# # 	# causal_prop <- job_vec[1]
# # 	# GTEx_He2 <- job_vec[2]
# # 	# ROSMAP_He2 <- job_vec[3]
# # 	job_info <- as.numeric(strsplit(job, '_')[[1]])
# # 	ROSMAP_causal_prop  <- job_info[1]
# # 	GTEx_causal_prop  <- job_info[2]
# # 	GTEx_overlap  <- job_info[3]
# # 	ROSMAP_He2  <- job_info[4]
# # 	GTEx_He2  <- job_info[5]

# # 	# should be 1:1000; included for testing with fewer
# # 	i_pos <- length(strsplit(sim_data$TargetID[[1]], '_')[[1]])
# # 	i_targets <- sort(unique(as.integer(
# # 		do.call(rbind, strsplit(sim_data$TargetID, '_'))[,i_pos])))

# # 	# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
# # 	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %:%
# # 	foreach(i=i_targets, .combine=rbind) %dopar% {

# # 		# get target id
# # 		target_id <- paste0(job, '_', i)

# # 		# subset sim_data for target
# # 		target_data <- sim_data[sim_data$TargetID == target_id, -which(colnames(sim_data) == 'TargetID')]

# # 		# variables for simulating phenotype
# # 		##var_expr <- var(target_data$sim_expr)
# # 		gamma <- sqrt(pheno_h2 / expr_var[[target_id]])

# # 		error_term <- rnorm(nrow(target_data), mean=0, sd=sqrt(1 - pheno_h2))

# # 		# Simulate phenotype from test expression
# # 		target_data$sim_pheno <- gamma * target_data$sim_expr + error_term

# # 		# output for target
# # 		out_i <- data.frame(
# # 			'job'=job,
# # 			'i'=i,
# # 			'TargetID'=target_id, 
# # 			'ROSMAP_causal_prop'=ROSMAP_causal_prop,
# # 			'GTEx_causal_prop'=GTEx_causal_prop,
# # 			'GTEx_overlap'=GTEx_overlap,
# # 			'ROSMAP_He2'=ROSMAP_He2, 
# # 			'GTEx_He2'=GTEx_He2,
# # 			'Hp2'=pheno_h2,
# # 			rbind(
# # 				get_fit('ROSMAP', target_data),
# # 				get_fit('GTEx', target_data), 
# # 				get_fit('Naive', target_data),
# # 				get_fit('SR', target_data)))

# # 	} %seed% 1234567

# # 	# get power info
# # 	out_power <- data.frame(
# # 		'job'=job,
# # 		'ROSMAP_causal_prop'=ROSMAP_causal_prop,
# # 		'GTEx_causal_prop'=GTEx_causal_prop,
# # 		'GTEx_overlap'=GTEx_overlap,
# # 		'ROSMAP_He2'=ROSMAP_He2, 
# # 		'GTEx_He2'=GTEx_He2,
# # 		job_pwr(out_data))
# # 	out_power <- apply(out_power, 2, as.character)

# # 	# output pred results
# # 	out_data <- out_data[, out_pred_cols]
# # 	write.table(
# # 		out_data,
# # 		all_out_pred_path,
# # 		quote=FALSE,
# # 		append=TRUE,
# # 		row.names=FALSE,
# # 		col.names=FALSE,
# # 		sep='\t')

# # 	# output power results
# # 	out_power <- out_power[, out_power_cols]
# # 	write.table(
# # 		out_power,
# # 		all_out_power_path,
# # 		quote=FALSE,
# # 		row.names=FALSE,
# # 		col.names=FALSE,
# # 		append=TRUE,
# # 		sep='\t')

# # 	save(out_data, out_power, file=paste0(dir, 'power/data/out_power_pred_', job, '.Rdata'))
# # 	rm(sim_data, out_data, out_power)

# # 	# time elapsed
# # 	end_time <- Sys.time()
# # 	print(paste0('Computation time: ', end_time - start_time))
# # }
