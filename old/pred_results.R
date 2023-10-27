#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=F, digits=6, scipen=10)

# load libraries
library(doFuture)
library(foreach)
library(reshape2)
library(sqldf)

# set up parallel environment
doFuture::registerDoFuture()
ncores <- availableCores(methods = 'SGE')
# plan(future.batchtools::batchtools_sge)
plan(sequential)

# read in passed arguments
# args <- (commandArgs(TRUE))
# if(length(args)==0) {
#   stop("Error: No arguments supplied!")
# } else {
# 	chrm = as.character(args[[1]])
# }

### 
## directories, input depend on whith test data set
test_dat <- 'MAP'
# test_dat <- 'ROS'
test_sample_size <- list('MAP'=114,'ROS'=119)[[test_dat]]
ROSMAP_train_cohort <- list('MAP'='ROS', 'ROS'='MAP')[[test_dat]]

## directories
SR_dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/'
ROSMAP_dir <- '/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/'
pred_dir <- paste0(SR_dir, 'SR_TWAS_pred/Brain_', test_dat, '_test_data/')

## paths
sampleid_path <- paste0(ROSMAP_dir, 'sampleid/', test_dat, '_WGS_test_', test_sample_size, '_sampleid.txt')
expr_path <- paste0(ROSMAP_dir, 'expression/ROSMAP_WGS_b38_expr_TIGAR.txt')
out_path <- paste0(SR_dir, 'SR_TWAS_pred/pred_results_', test_dat, '_test_data.txt')
out_long_path <- paste0(SR_dir, 'SR_TWAS_pred/pred_results_', test_dat, '_test_data_long.txt')

## add header to out_path
cat(paste(c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','ROSMAP','GTEx','PredDB','Naive','Naive3','SR'), collapse='\t', sep=''), file=out_path, sep='\n')
cat(paste(c('CHROM','GeneStart','GeneEnd','TargetID','GeneName','cohort','R2'), collapse='\t', sep=''), file=out_long_path, sep='\n')

## columns
info_cols <- c('CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName')
sampleid <- read.table(sampleid_path, header = FALSE)$V1
keep_cols <- c(info_cols, sampleid)

## function to reshape the dataframes to combine
shapedat <- function(x, colname=NULL){
	if (is.null(colname)){
		colname <- deparse(substitute(x))
	}
	xinfo_cols <- intersect(info_cols, colnames(x))
	x$TargetID <- trimws(x$TargetID)
	x$GeneName <- trimws(x$GeneName)
	x$TargetID <- do.call(rbind, strsplit(x$TargetID, '.', fixed = TRUE))[,1]
	return(melt(x, id.vars=xinfo_cols, 
		variable.name='sample_id', value.name=colname))
}

## function to read in and shape prediction output
readshapepred <- function(path){
	colname <- gsub('_pred_path', '', deparse(substitute(path)))
	df <- read.table(path, header=TRUE, check.names=FALSE, fill=TRUE)[, keep_cols]
	df <- shapedat(df, colname)
	return(df)
}

## GET PREDICTION RESULTS FOR EACH CHROMOSOME

for (chrm in 1:22){
	print(paste('CHR:', chrm))
	## pred output
	ROSMAP_pred_path <- paste0(pred_dir, ROSMAP_train_cohort, '_CHR', chrm, '_Pred_GReX.txt')
	GTEx_pred_path <- paste0(pred_dir, 'GTEx_CHR', chrm, '_Pred_GReX.txt')
	PredDB_pred_path <- paste0(pred_dir, 'PredDB_CHR', chrm, '_Pred_GReX.txt')
	Naive_pred_path <- paste0(pred_dir, 'Naive_CHR', chrm, '_Pred_GReX.txt')
	Naive3_pred_path <- paste0(pred_dir, 'Naive3_CHR', chrm, '_Pred_GReX.txt')
	SR_pred_path <- paste0(pred_dir, 'SR_CHR', chrm, '_Pred_GReX.txt')

	## read in expression for chrom and fix colnames
	sql_query <- paste0('select ', paste('`', gsub('-', '.', keep_cols), '`', collapse=',', sep=''), ' from file where CHROM=', chrm)
	expr <- read.csv.sql(expr_path, sql=sql_query, header=TRUE, sep='\t', row.names=NULL); suppressWarnings(closeAllConnections())
	colnames(expr) <- gsub('.', '-', colnames(expr), fixed=TRUE)

	# center expression data
	expr[,sampleid] <- expr[,sampleid] - rowMeans(expr[,sampleid], na.rm=TRUE)

	# read in and shape prediction output
	ROSMAP_pred <- readshapepred(ROSMAP_pred_path)
	GTEx_pred <- readshapepred(GTEx_pred_path)
	PredDB_pred <- readshapepred(PredDB_pred_path)
	Naive_pred <- readshapepred(Naive_pred_path)
	Naive3_pred <- readshapepred(Naive3_pred_path)
	SR_pred <- readshapepred(SR_pred_path)

	## merge dataframes
	data <- merge(
		shapedat(expr),
		Reduce(function(...) merge(..., all=TRUE), list(ROSMAP_pred, GTEx_pred, PredDB_pred, Naive_pred, Naive3_pred, SR_pred)),
		all.x=TRUE)

	############################
	## GET RESULTS FOR CHRM

	## start time
	start_time <- Sys.time()

	## function to get fit and results safely in case of target failed in one of the datasets
	get_fit <- function(x, dat, wide=TRUE){
		out_final <- tryCatch(
			expr = {
				expr_r2 <- summary(lm(expr ~ get(x), data = dat))$r.squared
				if (wide){
					out <- setNames(expr_r2, x)
				} else {
					out <- data.frame('cohort' = x, 'expr_R2' = expr_r2)
				}
				return(out)
			},
			error = function(e){
				if (wide){
					out <- setNames(0, x)
				} else {
					out <- data.frame('cohort' = x, 'expr_R2' = 0)
				}
				return(out)
			})
		return(out_final)
	}

	# number of targets to compare
	targets <- unique(data$TargetID)
	N_targets <- length(targets)

	# do N_targets of: getting expression prediction R2
	out_data <- foreach(i = 1:N_targets, .combine = rbind) %dopar% {
		target_id <- targets[i]

		# subset data for target
		target_data <- data[data$TargetID == target_id, ]
		target_info <- target_data[1, info_cols]
		rownames(target_info) <- NULL

		out_i <- as.data.frame(c(target_info, 
				get_fit('ROSMAP', target_data),
				get_fit('GTEx', target_data),
				get_fit('PredDB', target_data),
				get_fit('Naive', target_data),
				get_fit('Naive3', target_data),
				get_fit('SR', target_data)))

	} %seed% 1234567

	# output R2 values for each target; append to file of all results
	write.table(
		out_data,
		out_path,
		quote=FALSE,
		row.names=FALSE,
		col.names=FALSE,
		append=TRUE,
		sep='\t')

	# output R2 values for each target; append to file of all results
	out_long_data <- melt(out_data, id.vars=info_cols, 
		variable.name='cohort', value.name='R2')
	write.table(
		out_long_data,
		out_long_path,
		quote=FALSE,
		row.names=FALSE,
		col.names=FALSE,
		append=TRUE,
		sep='\t')
}

pred_dat <- read.table(out_path, header=TRUE, check.names=FALSE, fill=TRUE)
pred_dat_long <- read.table(out_long_path, header=TRUE, check.names=FALSE, fill=TRUE)

colMeans(pred_dat[,6:11])
head(pred_dat)

####

datm <- melt(dat, measure.vars = c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'), variable.name='cohort', value.name='x')

datm$cohort <- factor(datm$cohort,
	levels=c('PredDB','GTEx',ROSMAP_train_cohort,'Naive'),
	labels = c(
		'PrediXcan_GTEx',
		'TIGAR_GTEx', 
		paste0('TIGAR_', ROSMAP_train_cohort), 
		'Naive')
	)


datm2 <- melt(pred_dat, measure.vars = c('ROSMAP','GTEx','PredDB','Naive','SR'), variable.name='cohort', value.name='x')


r2_cut <- 0.005


# only those greater greater than r2_cut
datm2_f1 <- datm2[datm2$x > r2_cut, ]
merge(aggregate(x ~ cohort, data=datm2_f1, median),
	aggregate(x ~ cohort, data=datm2_f1, mean),
	by='cohort')

# only targets where at least one of SR or other cohort is > r2_dat
# targets_r2 <- datm[datm$color != '4', 'TargetID']
datm2_f2 <- datm2[datm2$TargetID %in% targets_r2, ]
merge(aggregate(x ~ cohort, data=datm2_f2, median),
	aggregate(x ~ cohort, data=datm2_f2, mean),
	by='cohort')


save(pred_dat, pred_dat_long, file=paste0(SR_dir, 'SR_TWAS_pred/pred_results_', test_dat, '_test_data.Rdata'))

###############

# SR_pred <- read.table(SR_pred_path, header=TRUE, check.names=FALSE, fill=TRUE)[, keep_cols]
# Naive_pred <- read.table(Naive_pred_path, header=TRUE, check.names=FALSE, fill=TRUE)[, keep_cols]
# GTEx_pred <- read.table(GTEx_pred_path, header=TRUE, check.names=FALSE, fill=TRUE)[, keep_cols]
# ROSMAP_pred <- read.table(ROSMAP_pred_path, header=TRUE, check.names=FALSE, fill=TRUE)[, keep_cols]

# # output for target (long version)
# out_i <- cbind(target_info, 
# 	rbind(
# 		get_fit('GTEx', target_data, wide=FALSE),
# 		get_fit('ROSMAP', target_data, wide=FALSE),
# 		get_fit('Naive', target_data, wide=FALSE),
# 		get_fit('SR', target_data, wide=FALSE)))


# # write to chrmosome chrm, only 
# write.table(
# 	out_data,
# 	out_path,
# 	quote=FALSE,
# 	row.names=FALSE,
# 	col.names=TRUE,
# 	sep='\t')

###############
