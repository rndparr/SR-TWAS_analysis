#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
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

######################################################

## FUNCTION
# function to get fit and results safely in case of target failed in one of the datasets; just get pheno

get_fit2 <- function(x, dat, sig_level=2.5*10^-6){
	out_final <- tryCatch(
		expr = {
			out <- as.integer(as.numeric(summary(lm(sim_pheno ~ get(x), data = dat))$coefficients[2, 'Pr(>|t|)']) < sig_level)
			return(out)
		},
		error = function(e){
			# print(e)
			out <- NA
			return(out)
		})
	return(out_final)
}

get_fit3 <- function(x, dat, ...){
	out_final <- tryCatch(
		expr = {
			out <- as.numeric(summary(lm(sim_pheno ~ get(x), data = dat))$coefficients[2, 'Pr(>|t|)'])
			return(out)
		},
		error = function(e){
			# print(e)
			out <- NA
			return(out)
		})
	return(out_final)
}


######
## PRED FILES/COLS
pred_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/target_files/'
# pred_cols <- colnames(read.table(paste0(pred_dir, 'header.txt'), header=TRUE))
pred_cols <- c('sample_id', 'sim_expr', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')
pred_col_classes <- c('character', rep('numeric', 5))


######
## OUTPUT HEADER TO FILE
out_cols <- c('i', 'ROSMAP', 'GTEx_EN', 'Naive_EN', 'SR_EN')
cat(paste(out_cols, collapse='\t', sep=''), file=paste0(dir, 'results/header.txt'), sep='\n')


######
## set seed
set.seed(1234567)

## DO PER-TARGET CALCULATIONS
# for (target_id in 1:10^6) {



foreach(target_id=1:10^6, .combine=c) %dopar% {

	if ((target_id %% 1000) == 0) {print(target_id)}

	# load data
	pred <- read.table(paste0(pred_dir, target_id, '.txt'),
		header=FALSE, colClasses=pred_col_classes)
	colnames(pred) <- pred_cols

	# Simulate phenotype from N(0,1)
	pred$sim_pheno <- rnorm(nrow(pred), mean=0, sd=1)

	# output for target
	out_data <- data.frame(
		'i'=target_id, 
		'ROSMAP'=get_fit3('ROSMAP', pred, siglvl),
		'GTEx_EN'=get_fit3('GTEx_EN', pred, siglvl), 
		'Naive_EN'=get_fit3('Naive_EN', pred, siglvl),
		'SR_EN'=get_fit3('SR_EN', pred, siglvl))

	# output pred results
	write.table(
		out_data,
		paste0(dir, 'results/pred_sig_results.txt'),
		quote=FALSE,
		append=TRUE,
		row.names=FALSE,
		col.names=FALSE,
		sep='\t')

} %seed% 7654321





sig_level_dat <- data.frame(
	siglvl=c(10^-4, 10^-5, 10^-6, 2.5*10^-6),
	text=c('1e-4', '1e-5', '1e-6', '2.5e-6'),
	seed=c(8654321, 6654321, 5654321, 7654321))

# for (i in 1:3){

# 	siglvl <- sig_level_dat[i, 'siglvl']
# 	siglvltext <- sig_level_dat[i, 'text']
# 	seedi <- sig_level_dat[i, 'seed']

# 	foreach(target_id=1:10^6, .combine=c) %dopar% {

# 		if ((target_id %% 1000) == 0) {print(target_id)}

# 		# load data
# 		pred <- read.table(paste0(pred_dir, target_id, '.txt'),
# 			header=FALSE, colClasses=pred_col_classes)
# 		colnames(pred) <- pred_cols

# 		# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
# 		# pheno_h2 <- pheno_h2_list[[1]]

# 		# Simulate phenotype from N(0,1)
# 		pred$sim_pheno <- rnorm(nrow(pred), mean=0, sd=1)

# 		# output for target
# 		out_data <- data.frame(
# 			'i'=target_id, 
# 			'ROSMAP'=get_fit2('ROSMAP', pred, siglvl),
# 			'GTEx_EN'=get_fit2('GTEx_EN', pred, siglvl), 
# 			'Naive_EN'=get_fit2('Naive_EN', pred, siglvl),
# 			'SR_EN'=get_fit2('SR_EN', pred, siglvl))

# 		# output pred results
# 		write.table(
# 			out_data,
# 			paste0(dir, 'results/pred_sig_results_', siglvltext,'.txt'),
# 			quote=FALSE,
# 			append=TRUE,
# 			row.names=FALSE,
# 			col.names=FALSE,
# 			sep='\t')

# 	} %seed% 7654321
# }



# 	print(paste0('Computation time: ', end_time - start_time))




# ############

# ## SIM PHENO NON-PHENO FROM N(0,1)


##########
## SIM_EXPR

# # sampleid
# test_sampleids_path <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_test_800_sampleid.txt')

# # read in test sampleid data
# test_sampleids <- read.table(test_sampleids_path, header = FALSE)$V1


# # path to simulated expression
# sim_expr_path <- paste0(dir, 'ROSMAP_expr.txt')

# # read in number of columns, set column classes
# sim_expr_ncols <- length(colnames(read.table(sim_expr_path, header = TRUE, check.names = FALSE, nrow=1)))
# sim_expr_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', sim_expr_ncols - 5))

# # read-in
# sim_expr <- read.table(sim_expr_path, header=TRUE, check.names=FALSE, colClasses=sim_expr_colClasses, skipNul=TRUE)
# sampleid <- colnames(sim_expr)[-1]

# # center expression data
# sim_expr[, sampleid] <- sim_expr[, sampleid] - rowMeans(sim_expr[, sampleid], na.rm=TRUE)

# # get variance
# expr_var <- var(unlist(sim_expr[1, -1]))

# # remove unneeded variables
# rm(sim_expr_colClasses, sampleid, sim_expr)

# ## FUNCTION
# # function to get fit and results safely in case of target failed in one of the datasets
# get_fit <- function(x, dat){
# 	out_final <- tryCatch(
# 		expr = {
# 			expr_r2 <- cor(dat$sim_expr, dat[[x]])^2
# 			fit <- summary(lm(sim_pheno ~ get(x), data = dat))
# 			out <- data.frame('cohort'=x, 'expr_R2'=expr_r2, 
# 				'pheno_R2'=fit$r.squared * 100, 'pheno_pval'=fit$coefficients[2, 'Pr(>|t|)'])
# 			return(out)
# 		},
# 		error = function(e){
# 			# print(e)
# 			out <- data.frame('cohort'=x, 'expr_R2'=0, 'pheno_R2'=0, 'pheno_pval'=1)
# 			return(out)
# 		})
# 	return(out_final)
# }



# ## pre-calculate variables for simulating phenotype
# gamma <- setNames(sqrt(pheno_h2_list / expr_var), pheno_h2_list)
# error_sd <- setNames(sqrt(1 - pheno_h2_list), pheno_h2_list)

# #### 
# ## replaces:
# 	# # variables for simulating phenotype
# 	# gamma <- sqrt(pheno_h2 / expr_var)
# 	# # 
# 	# error_term <- rnorm(nrow(pred), mean=0, sd=sqrt(1 - pheno_h2))
# 	# # Simulate phenotype from test expression
# 	# pred$sim_pheno <- gamma * pred$sim_expr + error_term

# ## PHENO SIM
# for (target_id in 1:10^6) {

# 	target_seed <- 1234567 + target_id
# 	# load data
# 	pred <- read.table(paste0(pred_dir, target_id, '.txt'),
# 		header=FALSE, colClasses=pred_col_classes)
# 	colnames(pred) <- pred_cols

# 	# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
# 	# pheno_h2 <- pheno_h2_list[[1]]
# 	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %dopar% {

# 		# Simulate phenotype from test expression
# 		pred$sim_pheno <- gamma[[pheno_h2]] * pred$sim_expr + rnorm(nrow(pred), mean=0, sd=error_sd[[pheno_h2]])

# 		# output for target
# 		out_i <- data.frame(
# 			'TargetID'=target_id, 
# 			'Hp2'=pheno_h2,
# 			rbind(
# 				get_fit('ROSMAP', pred),
# 				get_fit('GTEx_EN', pred), 
# 				get_fit('Naive_EN', pred),
# 				get_fit('SR_EN', pred)))

# 	} %seed% target_seed

# 	# output pred results
# 	out_data <- out_data[, out_pred_cols]
# 	write.table(
# 		out_data,
# 		all_out_pred_path,
# 		quote=FALSE,
# 		append=TRUE,
# 		row.names=FALSE,
# 		col.names=FALSE,
# 		sep='\t')

# }


# ######################################################
# 	######################################################
# ######################################################


# # ######################################################
# # pred_dir <- '/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/split_files/'
# # datasets <- c('Naive_EN', 'SR_EN', 'GTEx_EN', 'ROSMAP')

# # read_pred <- function(dataset, pred_dir, test_sampleid){
# # 	pred_path <- paste0(pred_dir, dataset, '_', test_sampleid, '.txt')
# # 	pred_dataset <- read.table(pred_path, header=FALSE, colClasses=c('character', 'numeric'), fill=TRUE)
# # 	colnames(pred_dataset) <- c('TargetID', dataset)
# # 	return(pred_dataset)
# # }

# # test_sampleid <- test_sampleids[[1]]
# # dataset <- 'scrap'

# # pred <- data.frame()

# # for (test_sampleid in test_sampleids) {

# # 	sim_data <- Reduce(function(...) merge(..., all=TRUE), 
# # 		sapply(datasets, read_pred, pred_dir=pred_dir, test_sampleid=test_sampleid, simplify=FALSE))


# # 	i_targets <- sort(unique(as.integer(sim_data$TargetID)))

# # 	# for each pheno_h2, do i_targets of: getting expression prediction R2, phenotype prediction and power analysis simulations
# # 	pheno_h2 <- pheno_h2_list[[1]]
# # 	out_data <- foreach(pheno_h2=pheno_h2_list, .combine=rbind) %dopar% {
# # 	# %:%
# # 	# foreach(i=i_targets, .combine=rbind) %dopar% {

# # 		test_sampleid_sim_expr <- sim_expr[sim_expr$sample_id == test_sampleid, 'sim_expr']

# # 		# variables for simulating phenotype
# # 		gamma <- sqrt(pheno_h2 / expr_var)

# # 		error_term <- rnorm(nrow(sim_data), mean=0, sd=sqrt(1 - pheno_h2))

# # 		# Simulate phenotype from test expression
# # 		sim_data$sim_pheno <- gamma * test_sampleid_sim_expr + error_term



# # 	} %seed% 1234567

# # }


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
