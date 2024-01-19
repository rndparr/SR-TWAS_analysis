#!/usr/bin/env Rscript

# set options
Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=FALSE)

# load arguments
args=(commandArgs(TRUE))
cat(paste0('\nargs:\n'))
print(args)
cat()
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==1) {
	sim_dir = as.character(args[[1]])
	scenario_sim_i_suffix '0.1_0.1_0.5_0.1_0.1_667'
	suffix = ''
} else {
	sim_dir = as.character(args[[1]])
	scenario_sim_i_suffix as.character(args[[2]])
	suffix = as.character(args[[3]])
	
}

args=(commandArgs(TRUE))
cat(paste0('\nargs:\n'))
print(args)
cat()
if(length(args)==0) {
	stop("Error: No arguments supplied!")
} else if(length(args)==1) {
	sim_dir = as.character(args[[1]])
	scenario_sim_i_suffix '0.1_0.1_0.5_0.1_0.1_667'
} else {
	sim_dir = as.character(args[[1]])
	scenario_sim_i_suffix as.character(args[[2]])
}


# load libraries
library(reshape2)
library(doFuture)
library(foreach)

# set up parallel environment
registerDoFuture()

pln <- plan(multicore, 
	workers=availableCores(methods='SGE'))
print(paste0('Starting job with ', availableCores(methods='SGE'), ' cores'))

######################################################
# directory for simulation files
out_dir <- paste0(sim_dir, '/type1_error/', scenario_sim_isuffix, '/')

######################################################
# read in test sampleid data
test_sampleids <- read.table(paste0(sim_dir, 'sampleid/ROSMAP_test_800_sampleid.txt'), header=FALSE)$V1


## FUNCTION
# function to get fit and results safely in case of target failed in one of the datasets; just get pheno
get_fit <- function(x, dat, ...){
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
pred_dir <- paste0(out_dir, 'pred/target_files/')
pred_cols <- c('sample_id', 'sim_expr', 'TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg')
pred_col_classes <- c('character', rep('numeric', 7))


######
## OUTPUT HEADER TO FILE
out_cols <- c('i', 'TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg')
cat(paste(out_cols, collapse='\t', sep=''), file=paste0(out_dir, 'results/header.txt'), sep='\n')


######
## set seed
set.seed(1234567)

## DO PER-TARGET CALCULATIONS
foreach(target_id=1:10^6, .combine=c) %dopar% {

	# for log file to keep track of ~what job the script is on
	if ((target_id %% 1000) == 0) {print(target_id)}

	# load data
	pred <- read.table(paste0(pred_dir, target_id, '.txt'),
		header=FALSE, colClasses=pred_col_classes, fill=TRUE)
	# deal with bad rows
	pred <- pred[pred$V1 %in% test_sampleids, ]
	colnames(pred) <- pred_cols

	# Simulate phenotype from N(0,1)
	pred$sim_pheno <- rnorm(nrow(pred), mean=0, sd=1)

	# output for target
	out_data <- data.frame(
		'i'=target_id, 
		'TIGAR-ROSMAP'=get_fit('TIGAR-ROSMAP', pred, siglvl),
		'PrediXcan-GTEx'=get_fit('PrediXcan-GTEx', pred, siglvl), 
		'TIGAR-ROSMAP_valid'=get_fit('TIGAR-ROSMAP_valid', pred, siglvl),
		'Naive'=get_fit('Naive', pred, siglvl),
		'SR'=get_fit('SR', pred, siglvl),
		'Avg'=get_fit('Avg', pred, siglvl), check.names=FALSE)
	# 'TIGAR-ROSMAP', 'PrediXcan-GTEx', 'TIGAR-ROSMAP_valid', 'Naive', 'SR', 'Avg'
	# print(out_data)

	# # output pred results
	write.table(
		out_data,
		paste0(out_dir, 'results/pred_sig_results.txt'),
		quote=FALSE,
		append=TRUE,
		row.names=FALSE,
		col.names=FALSE,
		sep='\t')

} %seed% 7654321


# sig_level_dat <- data.frame(
# 	siglvl=c(10^-4, 10^-5, 10^-6, 2.5*10^-6),
# 	text=c('1e-4', '1e-5', '1e-6', '2.5e-6'))
