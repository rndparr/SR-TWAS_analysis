#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)

# load libraries
library(reshape2)

# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/'

# sampleid
test_sampleid_path <- paste0('/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleid <- read.table(test_sampleid_path, header = FALSE)$V1

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

# memory usage of all objects 
memuse <- function() {
	for (object in ls(envir=sys.frame())){
		print(paste0(object, ':    ', format(object.size(eval(parse(text=object))), standard='SI', unit='auto')))
	}
}



##################
# GET SIM_EXPR, EXPR_VAR
pred_path <- paste0(dir, 'pred/ROSMAP_pred.txt')
# pred_path <- paste0(dir, 'pred/GTEx_EN_pred.txt')
# pred_path <- paste0(dir, 'pred/Naive_pred_EN.txt')
# pred_path <- paste0(dir, 'pred/SR_pred_EN.txt')

pred_cols <- colnames(read.table(pred_path, header=TRUE, check.names=FALSE, nrow=1))
pred_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', length(pred_cols) - 7))
pred <- shapedat(read.table(pred_path, header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE))
rm(pred)

pred <- shapedat(read.table(pred_path, header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE, nrows=50))

pred <- read.table(pred_path, header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE, nrows=50)
pred <- shapedat(pred)
# sr_pred <- shapedat(read.table(paste0(dir, 'pred/SR_pred_EN.txt'), header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE, nrow=20), 'SR')

# naive_pred <- shapedat(read.table(paste0(dir, 'pred/Naive_pred_EN.txt'), header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE, nrow=20), 'Naive')



# 	sim_data <- merge(
# 		sim_expr,
# 		Reduce(function(...) merge(..., all=TRUE), list(naive_pred, sr_pred)),
# 		all=TRUE)

# head(sim_data)[,c(1:7, 805)]

# read in prediction data, shape
get_pred_job <- function(dataset, job){
	# ret_pred <- pred[pred$dataset == dataset, ]

	# ret_pred$oldjob <- sub('_[^_]+$', '', ret_pred$TargetID)

	# job_df <- do.call(rbind, lapply(strsplit(ret_pred$oldjob, '_'), newjobvec))
	# ret_pred$job <- paste(job_df[,1], job_df[,2], job_df[,3], job_df[,4], job_df[,5], sep='_')

	# ret_pred <- ret_pred[ret_pred$job == job, ]

	# ret_pred <- ret_pred[, -1]
	# ret_pred <- ret_pred[, -1]

	ret_pred <- pred[which((pred$dataset == dataset) & (pred$job == job)), ]

	ret_pred <- ret_pred[, which(!(colnames(ret_pred)) %in% c('oldjob', 'dataset', 'suffix'))]
	# return(shapedat(ret_pred, paste0(dataset, '_pred')))
	return(shapedat(ret_pred, paste0(sub('_[^_]+$', '', dataset), '_pred')))
}


## for suffix jobs
## SIM_EXPR
sim_expr <- data.frame()

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

sim_expr$TargetID <- sub('0.1_0.1_0.5_0.1_0.1', '', sim_expr$TargetID)


sub('_(?!.*_)', '', sim_expr$TargetID, perl=TRUE)

all_sim_expr <- rbind(all_sim_expr, sim_expr)


save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr_noshape', out_suf, '.Rdata'))


# save unshaped dataset with only test sampleids
all_sim_expr_old <- all_sim_expr
all_sim_expr <- all_sim_expr[, c('TargetID', 'job', test_sampleid)]
save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr_noshape_testsampleid', out_suf, '.Rdata'))
all_sim_expr <- all_sim_expr_old; rm(all_sim_expr_old)

# reshape
all_sim_expr <- shapedat(all_sim_expr[,c('TargetID', 'job', test_sampleid)], 'sim_expr')

# expression variance (from ALL samples) for each simulation; used to get gamma for phenotype
expr_var <- by(all_sim_expr$sim_expr, all_sim_expr$TargetID, var, simplify=FALSE)

# filter sim_expr for only test_sampleids
all_sim_expr <- all_sim_expr[all_sim_expr$sample_id %in% test_sampleid, ]
all_sim_expr$sample_id <- droplevels(all_sim_expr$sample_id)
rownames(all_sim_expr) <- NULL

rm(list=c('expr_dir', 'test_sampleid_path', 'sim_expr_path', 'sim_expr_colClasses', 'sampleid')); # gc()

# save
save(expr_var, file=paste0(dir, 'power/data/expr_var', out_suf, '.Rdata'))
save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr', out_suf, '.Rdata'))

rm(expr_var, all_sim_expr)

#########
## EXPORT DATA FOR EACH JOB
for (job in jobs) {
	load(paste0(dir, 'power/data/sim_expr_noshape_testsampleid', out_suf,  '.Rdata'))
	sim_expr <- all_sim_expr; rm(all_sim_expr)
	sim_expr <- sim_expr[sim_expr$job == job, ]

	# reshape
	sim_expr <- shapedat(sim_expr)

	ROSMAP_pred <- get_pred_job('ROSMAP', job)
	GTEx_pred <- get_pred_job('GTEx', job)
	Naive_pred <- get_pred_job('Naive', job)
	SR_pred <- get_pred_job('SR', job)

	sim_data <- merge(
		sim_expr,
		Reduce(function(...) merge(..., all=TRUE), list(ROSMAP_pred, GTEx_pred, Naive_pred, SR_pred)),
		all=TRUE)
	rm(list=c('ROSMAP_pred', 'GTEx_pred', 'Naive_pred', 'SR_pred', 'sim_expr'))


	job_info <- strsplit(job, '_')[[1]]
	sim_data$ROSMAP_causal_prop <- job_info[1]
	sim_data$GTEx_causal_prop <- job_info[2]
	sim_data$GTEx_overlap <- job_info[3]
	sim_data$ROSMAP_He2 <- job_info[4]
	sim_data$GTEx_He2 <- job_info[5]

	sim_data <- sim_data[, c('ROSMAP_causal_prop', 'GTEx_causal_prop', 'GTEx_overlap', 'ROSMAP_He2', 'GTEx_He2', 'job', 'TargetID', 'sample_id', 'sim_expr', 'GTEx_pred', 'ROSMAP_pred', 'Naive_pred', 'SR_pred')]

	# job_info <- do.call(rbind, strsplit(sim_data$job, '_'))
	# sim_data$causal_prop <- job_info[, 1]
	# sim_data$GTEx_He2 <- job_info[, 2]
	# sim_data$ROSMAP_He2 <- job_info[, 3]

	# sim_data <- sim_data[, c('causal_prop', 'GTEx_He2', 'ROSMAP_He2', 'job', 'TargetID', 'sample_id', 'sim_expr', 'GTEx_pred', 'ROSMAP_pred', 'Naive_pred', 'SR_pred')]

	# save(sim_data, file=paste0(dir, 'power/data/sim_data_', job, suffix, '.Rdata'))
	save(sim_data, file=paste0(dir, 'power/data/sim_data_', job, '.Rdata'))
	rm(sim_data)
}


