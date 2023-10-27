#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
.libPaths('/home/rparrish/R/x86_64-redhat-linux-gnu-library/3.6')
options(stringsAsFactors=FALSE)

# load libraries
library(reshape2)

# load grguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else {
	suffix = as.character(args[[1]])
	jobs = as.character(args[[2]])
}



# directory for simulation files
dir <- '/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/'
expr_dir <- 'raw_dosage/'

# sampleid
test_sampleid_path <- paste0(dir, 'sampleid/', 'ROSMAP_test_800_sampleid.txt')

# read in test sampleid data
test_sampleid <- read.table(test_sampleid_path, header = FALSE)$V1


# do_jobs <- c('0.001_0.1_0.1', '0.001_0.2_0.2', '0.001_0.5_0.5', '0.01_0.1_0.1', 
# '0.01_0.2_0.2', '0.01_0.5_0.5', '0.05_0.1_0.1', '0.05_0.2_0.2', 
# '0.05_0.5_0.5', '0.1_0.1_0.1', '0.1_0.2_0.2', '0.1_0.5_0.5', 
# '0.001_0.05_0.1', '0.001_0.1_0.2', '0.001_0.25_0.5', '0.01_0.05_0.1', 
# '0.01_0.1_0.2', '0.01_0.25_0.5', '0.05_0.05_0.1', '0.05_0.1_0.2', 
# '0.05_0.25_0.5', '0.1_0.05_0.1', '0.1_0.1_0.2')


# do_jobs <- c('0.1_0.25_0.5')

# do_jobs <- c( '0.1_0.1_0.5_0.5_0.5', '0.1_0.1_0.5_0.2_0.2')

do_jobs <- c('0.001_0.001_0.5_0.1_0.1', '0.001_0.001_0.5_0.2_0.2', '0.001_0.001_0.5_0.5_0.5', 
'0.01_0.01_0.5_0.1_0.1', '0.01_0.01_0.5_0.2_0.2', '0.01_0.01_0.5_0.5_0.5', 
'0.05_0.05_0.5_0.1_0.1', '0.05_0.05_0.5_0.2_0.2', '0.05_0.05_0.5_0.5_0.5', 
'0.1_0.1_0.5_0.1_0.1')

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

newjobvec <- function(x){
	if(length(x) == 3){
		# (causal_prop for both)_(GTEx_He2)_(ROSMAP_He2) ->
		# (ROSMAP_causal_prop)_(GTEx_causal_prop)_(GTEx_overlap)_(ROSMAP_He2)_(GTEx_He2)
		return(c(x[[1]], x[[1]], '1', x[[3]], x[[2]]))
	} else {
		return(x)
	}
}


##################
# GET SIM_EXPR, EXPR_VAR
# path <- paste0(dir, 'pred/', dataset, '_pred.txt')
# path <- paste0(dir, 'pred/', dataset, '_pred_all.txt')
pred_path <- paste0(dir, 'pred/all_pred_results.txt')
pred_cols <- colnames(read.table(pred_path, header=TRUE, check.names=FALSE, nrow=1))
pred_colClasses <- c(rep('character', 2), rep('NULL', 3), 'character', 'NULL', rep('numeric', length(pred_cols) - 7))
pred <- read.table(pred_path, header=TRUE, check.names=FALSE, colClasses=pred_colClasses, skipNul=TRUE)
pred$oldjob <- sub('_[^_]+$', '', pred$TargetID)
job_df <- do.call(rbind, lapply(strsplit(pred$oldjob, '_'), newjobvec))
pred$job <- paste(job_df[,1], job_df[,2], job_df[,3], job_df[,4], job_df[,5], sep='_')

# pred[which(pred$suffix %in% c('_EN', '_EN_123121')), 'dataset'] <- paste0(pred[which(pred$suffix %in% c('_EN', '_EN_123121')), 'dataset'], '_EN')
# pred[which(pred$suffix %in% c('_EN_120621','_EN_121721')), 'dataset'] <- paste0(pred[which(pred$suffix %in% c('_EN_120621','_EN_121721')), 'dataset'], '_EN')

pred[which(pred$suffix %in% c('_EN','_EN_123121','_EN_120621','_EN_121721','_EN_overlap')), 'dataset'] <- paste0(pred[which(pred$suffix %in% c('_EN','_EN_123121','_EN_120621','_EN_121721','_EN_overlap')), 'dataset'], '_EN')

# unique(pred$dataset)

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
# suffix <- ''
# suffix <- '_120621'
# suffix <- '_121721'

# suffixes <- c('_120621', '_121721')


# suffixes <- c('_EN', '_EN_123121')
# suffix='_EN_123121'

# suffixes <- c('_EN_120621','_EN_121721')

suffixes <- '_EN_overlap'

## UPDATE JOB NAMING

out_EN <- 'EN_overlap'

## SIM_EXPR
jobs <- c()
all_sim_expr <- data.frame()
for (suffix in suffixes){
	# path to simulated expression
	suffix2 <- gsub('_EN', '', suffix)
	sim_expr_path <- paste0(dir, 'expression/raw_dosage/ROSMAP_expr', suffix2, '.txt')
	# read in number of columns, set column classes
	sim_expr_ncols <- length(colnames(read.table(sim_expr_path, header = TRUE, check.names = FALSE, nrow=1)))
	sim_expr_colClasses <- c(rep('NULL', 3), 'character', 'NULL', rep('numeric', sim_expr_ncols - 5))
	# read-in
	sim_expr <- read.table(sim_expr_path, header=TRUE, check.names=FALSE, colClasses=sim_expr_colClasses, skipNul=TRUE)
	sampleid <- colnames(sim_expr)[-1]
	# center expression data
	sim_expr[,sampleid] <- sim_expr[,sampleid] - rowMeans(sim_expr[,sampleid], na.rm=TRUE)

	sim_expr$oldjob <- sub('_[^_]+$', '', sim_expr$TargetID)

	# FILTER OUT UNFINISHED JOB
	sim_expr <- sim_expr[which(sim_expr$oldjob %in% do_jobs), ]

	# get job name
	job_df <- do.call(rbind, lapply(strsplit(sim_expr$oldjob, '_'), newjobvec))
	sim_expr$job <- paste(job_df[,1], job_df[,2], job_df[,3], job_df[,4], job_df[,5], sep='_')

	## get job info columns
	sim_expr$ROSMAP_causal_prop <- job_df[,1]
	sim_expr$GTEx_causal_prop <- job_df[,2]
	sim_expr$GTEx_overlap <- job_df[,3]
	sim_expr$ROSMAP_He2 <- job_df[,4]
	sim_expr$GTEx_He2 <- job_df[,5]

	jobs <- c(jobs, unique(sim_expr$job))

	all_sim_expr <- rbind(all_sim_expr, sim_expr)
	rm(sim_expr)
}
save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr_noshape', out_EN,'.Rdata'))

# save unshaped dataset with only test sampleids
all_sim_expr_old <- all_sim_expr
all_sim_expr <- all_sim_expr[, c('TargetID', 'job', test_sampleid)]
save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr_noshape_testsampleid', out_EN,'.Rdata'))
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
save(expr_var, file=paste0(dir, 'power/data/expr_var', out_EN,'.Rdata'))
save(all_sim_expr, file=paste0(dir, 'power/data/sim_expr', out_EN,'.Rdata'))

rm(expr_var, all_sim_expr)

#########
## EXPORT DATA FOR EACH JOB
for (job in jobs) {
	load(paste0(dir, 'power/data/sim_expr_noshape_testsampleid', out_EN, '.Rdata'))
	sim_expr <- all_sim_expr; rm(all_sim_expr)
	sim_expr <- sim_expr[sim_expr$job == job, ]

	# reshape
	sim_expr <- shapedat(sim_expr)

	ROSMAP_pred <- get_pred_job('ROSMAP', job)
	GTEx_pred <- get_pred_job('GTEx_EN', job)
	Naive_pred <- get_pred_job('Naive_EN', job)
	SR_pred <- get_pred_job('SR_EN', job)

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
	save(sim_data, file=paste0(dir, 'power/data/sim_data_', job, '_EN.Rdata'))
	rm(sim_data)
}

